###########################
### Integration methods ###
###########################

library(Seurat)
library(tidyverse)
library(liger)
library(SingleCellExperiment)
library(MultiAssayExperiment)
library(magrittr)
source("~/multiOmic_benchmark/preprocess/selectFeatures.R")

### Wrapper function ###
run_integration <- function(sce.list, method, n_features, reference="RNA", query="ATAC"){
  ## Select integration features
  integrate_features <- HVG_Seurat(sce.list[[reference]], nfeatures = n_features) %>% {.[. %in% rownames(sce.list[[query]])]}
  
  if (method == "CCA") {
    integrate <- integrate_seuratCCA
  } else if (method == "liger") {
    integrate <- integrate_liger
  } else {
    stop("invalid integration method. Please select one of CCA, liger")
  }
  
  int_output <- integrate(sce.list, integrate_features, reference=reference, query=query)
  int_output[["integrate_features"]] <- integrate_features
  return(int_output)  
  }


### Integration models ###

#' Seurat CCA integration
#' @param sce.list list of SingleCellExperiment objects for RNA and ATAC seq datasets
#' @param integrate_features selected features to perform CCA on
#' @param reference reference dataset for FindTransferAnchors
#' @param query query dataset for FindTransferAnchors
#' 
#' @return list of integration output (see details)
#' 
#' @details Function outputs list of 
#' 1) MultiAssayExperiment object containing input dataset + integrated (imputed) dataset for RNA, 
#' 2) misc= CCA model (Anchor object from Seurat)
integrate_seuratCCA <- function(sce.list, integrate_features, reference="RNA", query="ATAC"){
  seurat.list <- imap(sce.list, ~ as.Seurat(.x, assay=.y))
  ## Scale data
  seurat.list <- map(seurat.list, ~ ScaleData(.x))
  ## Calculate CCA anchors
  transfer.anchors <- FindTransferAnchors(reference = seurat.list[[reference]], query = seurat.list[[query]], 
                                          features = integrate_features, 
                                          reduction = "cca")
  ## Impute transcriptome for ATAC-seq cells
  refdata <- GetAssayData(seurat.list[[reference]], slot = "scale.data")
  imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, 
                             weight.reduction = 'cca')
  
  ## Merge imputed RNA data with real RNA data
  imputed.sce <- SingleCellExperiment(list(normcounts=imputation@data), 
                                      colData=colData(sce.list[[query]]))
  sce.rna.bind <- SingleCellExperiment(list(normcounts=refdata), 
                                       colData=colData(sce.list[[reference]]))
  merged.sce <- cbind(sce.rna.bind, imputed.sce)
  colData(merged.sce)[["tech"]] <- ifelse(colnames(merged.sce) %in% colnames(sce.rna.bind), 
                                          reference, 
                                          query)
  ## Prepare output object
  sce.list[[str_c("integrated.", reference)]] <- merged.sce
  intMAE <- MultiAssayExperiment(sce.list, colData = colData(merged.sce))
  misc <- list(transfer.anchors = transfer.anchors)
  return(list(intOut=intMAE, misc=misc))
  }

#' LIGER NMF integration
#' @param sce.list list of SingleCellExperiment objects for RNA and ATAC seq datasets
#' @param integrate_features selected features to perform CCA on
#' 
#' @return list of integration output (see details)
#' 
#' @details Function outputs list of 
#' 1) MultiAssayExperiment object containing input dataset + integrated (imputed) dataset for RNA, 
#' 2) misc
integrate_liger <- function(sce.list, integrate_features, reference="RNA", query="ATAC"){
  data.list <- map(sce.list, ~ assay(.x, "counts"))
  liger.obj <- createLiger(data.list)
  liger.obj@norm.data <- map(sce.list, ~ assay(.x, "cpm"))
  liger.obj@var.genes <- integrate_features
  ## Scale data
  liger.obj <- scaleNotCenter(liger.obj)
  ## Perform matrix factorization
  liger.obj <- optimizeALS(liger.obj, k=20, lambda = 5.0)
  ## Quantile normalization step
  liger.obj <- quantileAlignSNF(liger.obj)
  ## Imputation: reconstruct RNA data for ATAC cells using weight matrix from NMF model
  Y_ref <- liger.obj@scale.data[[reference]]
  V_ref <- liger.obj@V[[reference]]
  H_query_cells <- liger.obj@H[[query]]
  H_ref_cells <- liger.obj@H.norm
  W <- liger.obj@W
  
  # Y_ref_cells <- H_ref_cells %*% (V_ref + W)
  # Y_ref_cells <- H_ref_cells %*% (W)
  # Y_query_cells <- H_query_cells %*% (V_ref + W)
  Y_impute <- H_ref_cells %*% (W + V_ref) ## Using the normalized factor matrix
    
  # Y_all <- rbind(Y_ref, Y_query_cells)
  # 
  # merged.sce <- SingleCellExperiment(list(normcounts=Y_all), 
  #                                   colData=colData(sce.list$ATAC))
  # Y_smp <- Y_all[, sample(1:3198, 100)]
  # 
  # Y_smp %>% melt(varnames=c("cell", "gene")) %>%
  #   mutate(dataset=ifelse(str_detect(cell, "-"), "ATAC", "RNA")) %>%
  #   ggplot(aes(dataset,value, color=dataset )) + geom_quasirandom()
  
  ## Merge imputed RNA data with real RNA data
  imputed.sce <- SingleCellExperiment(list(normcounts = t(Y_impute)), 
                                      colData=rbind(colData(sce.list[[reference]]),colData(sce.list[[query]])))
  merged.sce <- imputed.sce
  colData(merged.sce)[["tech"]] <- ifelse(colnames(merged.sce) %in% colnames(sce.rna.bind), 
                                          reference, 
                                          query)
  ## Prepare output object
  sce.list[[str_c("integrated.", reference)]] <- merged.sce
  intMAE <- MultiAssayExperiment(sce.list, colData = colData(merged.sce))
  misc <- list(H = liger.obj@H, H.norm = liger.obj@H.norm, V= liger.obj@V, W=liger.obj@W, liger.obj@parameters)
  return(list(intOut = intMAE, misc = misc))
  }
# 
# sce.list <- readRDS("~/my_data/integrated_thymus/F74_SCElist_20191017.RDS")
# f74.list <- sce.list
# sce.list %<>% map(~ .x[,1:500])
# integrate_features <- HVG_Seurat(sce.list$RNA, nfeatures = 4000) %>% {.[. %in% rownames(sce.list$ATAC)]}
# 
# int_seurat <- integrate_seuratCCA(sce.list, integrate_features)
# int_liger <- integrate_liger(sce.list, integrate_features)

# 
# liger.obj <- liger::runTSNE(liger.obj, use.raw = T)
# p1 <- plotByDatasetAndCluster(liger.obj, return.plots = T)
# print(p1[[1]])
# 
# liger.obj <- liger::runTSNE(liger.obj)
# p1 <- plotByDatasetAndCluster(liger.obj, return.plots = T)
# print(p1[[1]])
# 
# pheatmap::pheatmap(liger.obj@H.norm)
# 
# CD8A = plotGene(liger.obj,gene="CD4",return.plots=T)
# plot_grid(CD8A[[1]] + ylim(-30, 30), CD8A[[2]] + ylim(-30, 30))
# 
# 
# p1[[1]]










