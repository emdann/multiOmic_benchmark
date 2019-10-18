###########################
### Integration methods ###
###########################

library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(MultiAssayExperiment)
source("multiOmic_benchmark/selectFeatures.R")

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
  ## Calculate CCA anchors
  transfer.anchors <- FindTransferAnchors(reference = seurat.list[[reference]], query = seurat.list[[query]], 
                                          features = integrate_features, 
                                          reduction = "cca")
  ## Impute transcriptome for ATAC-seq cells
  refdata <- GetAssayData(seurat.list[[reference]], slot = "counts")
  imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, 
                             weight.reduction = 'cca')
  
  ## Merge imputed RNA data with real RNA data
  imputed.sce <- SingleCellExperiment(list(logcounts=imputation@data), 
                                      colData=colData(sce.list$ATAC))
  sce.rna.bind <- SingleCellExperiment(list(logcounts=assay(sce.list$RNA, "logcounts")), 
                                       colData=colData(sce.list$RNA))
  merged.sce <- cbind(sce.rna.bind, imputed.sce)
  colData(merged.sce)[["tech"]] <- ifelse(str_detect(colnames(merged.sce), "_"), 
                                          "RNA", 
                                          "ATAC")
  ## Prepare output object
  sce.list[["integrated.RNA"]] <- merged.sce
  intMAE <- MultiAssayExperiment(sce.list, colData = colData(merged.sce))
  misc <- list(transfer.anchors = transfer.anchors)
  return(list(intOut=intMAE, misc=misc))
  }

