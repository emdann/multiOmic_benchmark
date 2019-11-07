###########################
### Integration methods ###
###########################

suppressPackageStartupMessages({
library(Seurat)
library(tidyverse)
library(liger)
library(conos)
library(pdist)
library(SeuratWrappers)
library(SingleCellExperiment)
library(MultiAssayExperiment)
library(magrittr)
library(FastKNN)
source("~/multiOmic_benchmark/preprocess/selectFeatures.R")
source("~/multiOmic_benchmark/utils.R")
})

### Wrapper function ###
run_integration <- function(sce.list, method, n_features, feature.selection = "union_hvg", reference="RNA", query="ATAC"){
  ## Select integration features
  if (feature.selection=="reference_hvg") {
    integrate_features <- HVG_Seurat(sce.list[[reference]], nfeatures = n_features) %>% {.[. %in% rownames(sce.list[[query]])]}
  } else if (feature.selection == "union_hvg"){
    integrate_features_ref <- HVG_Seurat(sce.list[[reference]], nfeatures = n_features) %>% {.[. %in% rownames(sce.list[[query]])]}
    integrate_features_query <- HVG_Seurat(sce.list[[query]], nfeatures = n_features) %>% {.[. %in% rownames(sce.list[[reference]])]}
    integrate_features <- union(integrate_features_ref, integrate_features_query)
  } else {
      stop("Invalid feature selection method: please specify one of 'reference_hvg' or 'union_hvg' ")
    }
  if (method == "CCA") {
    integrate <- integrate_seuratCCA
  } else if (method == "liger") {
    integrate <- integrate_liger
  } else {
    stop("invalid integration method. Please specify one of 'CCA', 'liger'")
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


### Transfer Labels ###
#' Label transfer from Seurat CCA
#' @param seurat.list list of Seurat objects used as input for CCA integration
#' @param transfer_anchor output of Seurat::FindTransferAnchors for the seurat objects in input
#' @param reference reference dataset for FindTransferAnchors
#' @param query query dataset for FindTransferAnchors
#' 
#' @return list of Seurat objects with predicted annotations in the query dataset metadata
#' 
labelTransfer_seuratCCA <- function(seurat.list, transfer.anchors, reference="RNA", query="ATAC"){
  ## Transfer cell type labels
  celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = seurat.list[[reference]]$annotation, weight.reduction = "cca")
  seurat.list[[query]] <- AddMetaData(seurat.list[[query]], metadata = celltype.predictions)
  return(seurat.list)
}

#' Label transfer for LIGER
#' @param liger.obj liger object (output of run_liger)
#' @param sce.list list of SingleCellExperiment objects used as input for run_liger
#' @param annotation column of meta 
#' @param reference reference dataset for label transfer
#' @param query query dataset for label transfer
#' 
#' @return list of Seurat objects with predicted annotations (with score) in the query dataset metadata
#' 
#' @details Label transfer is performed by computing a cross-dataset k-nearest neighbor graph in the aligned factor space. 
#' Then we set the value of each missing label is set to the most abundant label among its k nearest neighbors in the reference dataset. 
#' Ties are dealt with by taking the first of the maximum values. As prediction score I take the number of neighbors harboring the predicted.id 
#' divided by k.
labelTransfer_liger <- function(liger.obj, sce.list, annotation.col="annotation", k=30, reference="RNA", query="ATAC"){
  seurat.list <- imap(sce.list, ~ as.Seurat(.x, assay=.y))
  seurat.list <- imap(seurat.list, ~ RenameCells(.x, add.cell.id=.y))
  annotation <- seurat.list[[reference]][[annotation.col]]
  
  ## Compute distance matrix between query and dataset 
  small.seu.liger <- ligerToSeurat(liger.obj, renormalize = F)
  ref.cells <- rownames(small.seu.liger@meta.data[which(small.seu.liger$orig.ident==reference),])
  query.cells <- rownames(small.seu.liger@meta.data[which(small.seu.liger$orig.ident==query),])
  nmf.mat <- small.seu.liger@reductions$inmf@cell.embeddings
  
  dist.mat <- as.matrix(pdist(nmf.mat[query.cells,], nmf.mat[ref.cells,]))
  rownames(dist.mat) <- query.cells
  colnames(dist.mat) <- ref.cells
  cross.dataset.NN <- matrix(0, nrow = nrow(dist.mat), ncol = ncol(dist.mat))
  for (i in 1:nrow(dist.mat)) {
    nn = k.nearest.neighbors(i, dist.mat, k=k)
    cross.dataset.NN[i,nn] <- 1
  }
  rownames(cross.dataset.NN) <- query.cells
  colnames(cross.dataset.NN) <- ref.cells
  # small.seu.liger@graphs <- cross.dataset.NN
  nn.list <- getNNlist(cross.dataset.NN, is.seurat = F)
  predicted.labels <- 
    map_dfr(nn.list, ~ propagateNNannotation(.x, annotation)) %>%
    mutate(cell=names(nn.list))
  query.metadata <- predicted.labels %>%
    filter(cell %in% colnames(small.seu.liger)[which(small.seu.liger$orig.ident==query)]) %>%
    column_to_rownames("cell") 
  seurat.list[[query]] <- AddMetaData(seurat.list[[query]], metadata = query.metadata)
  return(seurat.list)
}

labelTransfer_conos <- function(conos.out, sce.list, annotation.col="annotation", reference="RNA", query="ATAC"){
  seurat.list <- imap(sce.list, ~ as.Seurat(.x, assay=.y))
  ## Extract cell type annotation
  annotation <- seurat.list[[reference]][[annotation.col]]
  annotation <- setNames(annotation[,1], rownames(annotation))
  ## propagate labels
  new.label.probabilities <- conos.out$propagateLabels(labels = annotation, verbose = T)
  new.annot <- setNames(colnames(new.label.probabilities)[apply(new.label.probabilities,1,which.max)], rownames(new.label.probabilities))
  score <- apply(new.label.probabilities, 1, max)
  seurat.list[[query]] <- AddMetaData(seurat.list[[query]], metadata = data.frame(predicted.id = new.annot[colnames(seurat.list[[query]])],
                                                                                  score = score[colnames(seurat.list[[query]])]))
  return(seurat.list)
}

### Run models ### 
#' Seurat CCA integration model
#' @param sce.list list of SingleCellExperiment objects for RNA and ATAC seq datasets
#' @param integrate_features selected features to perform CCA on
#' @param reference reference dataset for FindTransferAnchors
#' @param query query dataset for FindTransferAnchors
#' 
#' @return list of integration output (see details)
#' 
#' @details Function outputs list of 
#' 1) model: CCA model (Anchor object from Seurat)
#' 2) input: Seurat object list containing input datasets
run_SeuratCCA <- function(sce.list, integrate_features, reference="RNA", query="ATAC"){
  seurat.list <- imap(sce.list, ~ as.Seurat(.x, assay=.y))
  seurat.list <- imap(seurat.list, ~ RenameCells(.x, add.cell.id=.y))
  ## Scale data
  seurat.list <- map(seurat.list, ~ ScaleData(.x))
  ## Calculate CCA anchors
  transfer.anchors <- FindTransferAnchors(reference = seurat.list[[reference]], query = seurat.list[[query]], 
                                          features = integrate_features, 
                                          reduction = "cca")
  return(list(model=transfer.anchors, input=seurat.list))
  }
#' LIGER NMF model

#' @param sce.list list of SingleCellExperiment objects for RNA and ATAC seq datasets
#' @param integrate_features selected features to perform CCA on
#' 
#' @return list of integration output (see details)
#' 
#' @details Function outputs list of 
#' 1) model: 
#' 2) misc
run_liger <- function(sce.list, integrate_features, reference="RNA", query="ATAC"){
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
  ## Run UMAP otherwise Seurat conversion fails
  liger.obj <- runUMAP(liger.obj)
  return(list(model = liger.obj, input=sce.list))
}

#' Run CONOS
#' 
run_conos <- function(sce.list, integrate_features, reference="RNA", query="ATAC"){
  data.processed <- map(sce.list, ~ as.Seurat(.x)) 
  VariableFeatures(data.processed[[reference]]) <- integrate_features
  VariableFeatures(data.processed[[query]]) <- integrate_features
  data.processed <- map(data.processed, ~ ScaleData(.x) %>% RunPCA(dims=1:30))
  l.con <- Conos$new(data.processed,n.cores=30)
  l.con$buildGraph(k=15,k.self=5,k.self.weigh=0.01,ncomps=30,n.odgenes=5e3,space='PCA') 

  l.con$findCommunities(resolution=1.5)
  l.con$embedGraph(alpha=1/2)
  return(list(model=l.con, input=sce.list))
}

### Utils ###

propagateNNannotation <- function(nn.vector, annotation){
  nn.annotations <- table(annotation[nn.vector,])
  if (length(which(nn.annotations==max(nn.annotations))) > 1) {
    warning("Tie in most abundant annotation. Selecting the first one.")
  }
  if (max(nn.annotations)==0) {
    predicted.id <- NA
  } else {
    predicted.id <- names(which.max(nn.annotations))
  }
  label.fraction <- max(nn.annotations)[1]/sum(nn.annotations)
  list(predicted.id=predicted.id, score=label.fraction)
}

## Running benchmark
# sce.list <- readRDS(file = "my_data/integrated_thymus/F74_SCElist_20191101.RDS")
# small.sce.list <- map(sce.list, ~ .x[,1:500])
# 
# reference='RNA'
# query="ATAC"
# n_features = 2000
# 
# ## Feature selection
# integrate_features_ref <- HVG_Seurat(sce.list[[reference]], nfeatures = n_features) %>% {.[. %in% rownames(sce.list[[query]])]}
# integrate_features_query <- HVG_Seurat(sce.list[[query]], nfeatures = n_features) %>% {.[. %in% rownames(sce.list[[reference]])]}
# integrate_features <- union(integrate_features_ref, integrate_features_query)
# 
# # integrate_features <- HVG_Seurat(sce.list[[reference]], nfeatures = n_features) %>% {.[. %in% rownames(sce.list[[query]])]}
# 
# ## Run models
# small.cca <- run_SeuratCCA(small.sce.list, integrate_features)
# small.liger <- run_liger(small.sce.list, integrate_features)
# small.conos <- run_conos(small.sce.list, integrate_features)
# 
# ## Transfer labels
# small.seu.cca <- labelTransfer_seuratCCA(seurat.list = small.cca$input,
#                                          transfer.anchors = small.cca$model)
# small.seu.liger <- labelTransfer_liger(liger.obj = small.liger$model,
#                                        sce.list = small.liger$input, 
#                                        k = 50)
# small.seu.conos <- labelTransfer_conos(conos.out = small.conos$model,sce.list = small.conos$input)
# 
# small.seu.cca
# small.seu.conos
# small.seu.liger
# # 
# # small.liger <- run_liger(small.sce.list, integrate_features)
# # small.seurat.liger <- labelTransfer_liger(liger.obj = small.liger$model, sce.list =  small.liger$input)
# # small.conos <- run_conos(small.sce.list, integrate_features)
# # seurat.list.conos <- labelTransfer_conos(small.conos$model, small.conos$input)
# 
# 
# 
