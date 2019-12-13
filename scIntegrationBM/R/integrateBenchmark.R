###########################
### Integration methods ###
###########################

# suppressPackageStartupMessages({
# library(Seurat)
# library(tidyverse)
# library(liger)
# library(conos)
# library(pdist)
# library(SeuratWrappers)
# library(SingleCellExperiment)
# library(MultiAssayExperiment)
# library(magrittr)
# library(FastKNN)
# source("~/multiOmic_benchmark/preprocess/selectFeatures.R")
# source("~/multiOmic_benchmark/utils.R")
# })

### Transfer Labels ###
#' Label transfer from Seurat CCA
#' @param seurat.list list of Seurat objects used as input for CCA integration
#' @param transfer_anchor output of Seurat::FindTransferAnchors for the seurat objects in input
#' @param reference reference dataset for FindTransferAnchors
#' @param query query dataset for FindTransferAnchors
#'
#' @return list of Seurat objects with predicted annotations in the query dataset metadata
#'
#' @import Seurat
#'
#' @export
labelTransfer_seuratCCA <- function(transfer.anchors, seurat.list, annotation.col="annotation", reference="RNA", query="ATAC", ... ){
  labels <- seurat.list[[reference]]@meta.data[,annotation.col]
  names(labels) <- rownames(seurat.list[[reference]]@meta.data)
  ## Transfer cell type labels
  celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = labels,
                                       ... )
  seurat.list[[query]] <- AddMetaData(seurat.list[[query]], metadata = celltype.predictions)
  pred.labels <- getPredictedLabels(seurat.list, int.name = "CCA", id.col = 'predicted.id', score.col = "prediction.score.max" )
  return(pred.labels)
}

#' Label transfer for LIGER
#' @param liger.obj liger object (output of run_liger)
#' @param seurat.list list of Seurat objects used as input for run_liger
#' @param annotation.col column of meta
#' @param reference reference dataset for label transfer
#' @param query query dataset for label transfer
#'
#' @return list of Seurat objects with predicted annotations (with score) in the query dataset metadata
#'
#' @details Label transfer is performed by computing a cross-dataset k-nearest neighbor graph in the aligned factor space.
#' Then we set the value of each missing label is set to the most abundant label among its k nearest neighbors in the reference dataset.
#' Ties are dealt with by taking the first of the maximum values. As prediction score I take the number of neighbors harboring the predicted.id
#' divided by k.
#'
#' @import liger
#' @import purrr
#' @import dplyr
#' @importFrom FastKNN k.nearest.neighbors
#' @importFrom Seurat as.Seurat
#' @importFrom pdist pdist
#'
#' @export
labelTransfer_liger <- function(liger.obj, seurat.list, annotation.col="annotation", k=50, reference="RNA", query="ATAC"){
  # seurat.list <- imap(sce.list, ~ as.Seurat(.x, assay=.y))
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

  nn.list <- getNNlist(cross.dataset.NN, is.seurat = F)
  predicted.labels <-
    map_dfr(nn.list, ~ propagateNNannotation(.x, annotation)) %>%
    dplyr::mutate(cell=names(nn.list))
  query.metadata <- predicted.labels %>%
    dplyr::filter(cell %in% colnames(small.seu.liger)[which(small.seu.liger$orig.ident==query)]) %>%
    column_to_rownames("cell")
  seurat.list[[query]] <- AddMetaData(seurat.list[[query]], metadata = query.metadata)
  pred.labels <- getPredictedLabels(seurat.list, int.name = "Liger", id.col = 'predicted.id', score.col="score")
  return(pred.labels)
}

### Transfer Labels ###
#' Label transfer from Conos
#' @param conos.out model output of `run_conos`
#' @param sce.list list of SingleCellExperiment objects used as input of `run_conos` (saved as `input` in output list of `run_conos`)
#' @param annotation.col column of meta
#' @param reference reference dataset for FindTransferAnchors
#' @param query query dataset for FindTransferAnchors
#'
#' @return list of Seurat objects with predicted annotations in the query dataset metadata
#'
#' @import Seurat
#'
#' @export
labelTransfer_conos <- function(conos.out, sce.list, annotation.col="annotation", reference="RNA", query="ATAC"){
  # seurat.list <- imap(sce.list, ~ as.Seurat(.x, assay=.y))
  ## Extract cell type annotation
  annotation <- seurat.list[[reference]][[annotation.col]]
  annotation <- setNames(annotation[,1], rownames(annotation))
  ## propagate labels
  new.label.probabilities <- conos.out$propagateLabels(labels = annotation, verbose = T, fixed.initial.labels=T)
  new.annot <- setNames(colnames(new.label.probabilities)[apply(new.label.probabilities, 1, function(x) which.max(x)[1]) ], rownames(new.label.probabilities))
  score <- apply(new.label.probabilities, 1, max)
  seurat.list[[query]] <- AddMetaData(seurat.list[[query]], metadata = data.frame(predicted.id = new.annot[colnames(seurat.list[[query]])],
                                                                                  score = score[colnames(seurat.list[[query]])]))
  pred.labels <- getPredictedLabels(seurat.list, int.name = "Conos", id.col = 'predicted.id', score.col="score")
  return(pred.labels)
}

### Run models ###

#' Seurat CCA integration model
#' @param sce.list list of SingleCellExperiment objects for multi-omic datasets
#' @param integrate_features selected features to perform CCA on
#' @param reference reference dataset for FindTransferAnchors
#' @param query query dataset for FindTransferAnchors
#'
#' @return list of integration output (see details)
#'
#' @details Function outputs list of
#' 1) model: CCA model (Anchor object from Seurat)
#' 2) input: Seurat object list containing input datasets
#'
#' @import Seurat
#'
#' @export
run_SeuratCCA <- function(seurat.list, integrate_features, reference="RNA", query="ATAC"){
    # seurat.list <- imap(sce.list, ~ as.Seurat(.x, assay=.y))
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
#' @param sce.list list of SingleCellExperiment objects for multi-omic datasets
#' @param integrate_features selected features to perform CCA on
#' @param reference reference dataset for FindTransferAnchors
#' @param query query dataset for FindTransferAnchors
#'
#' @return list of integration output (see details)
#'
#' @details Function outputs list of
#' 1) model:
#' 2) misc
#'
#' @import liger
#' @importFrom purrr map
#'
#' @export
run_liger <- function(seurat.list, integrate_features, reference="RNA", query="ATAC"){
  liger.obj <- seuratToLiger(seurat.list[c(reference,query)], names = c(reference,query), renormalize = F)
  liger.obj@norm.data <- map(seurat.list, ~ Matrix(.x@assays[[1]]@data))
  # data.list <- map(seurat.list, ~ assay(.x, "counts"))
  # liger.obj <- createLiger(data.list)
  # liger.obj@norm.data <- map(sce.list, ~ assay(.x, "logcounts"))
  # liger.obj@norm.data[[query]] <- Matrix(logcounts(sce.list[[query]]))
  ## Select genes
  integrate_features <- integrate_features[which(integrate_features %in% rownames(liger.obj@norm.data[[reference]]) &
                                                   integrate_features %in% rownames(liger.obj@norm.data[[query]]))]
  liger.obj@var.genes <- integrate_features
  ## Scale data
  liger.obj <- scaleNotCenter(liger.obj)
  ## Perform matrix factorization
  liger.obj <- optimizeALS(liger.obj, k=20, lambda = 5.0)
  ## Quantile normalization step
  liger.obj <- quantileAlignSNF(liger.obj)
  ## Run UMAP otherwise Seurat conversion fails
  liger.obj <- runUMAP(liger.obj)
  return(list(model = liger.obj, input=seurat.list))
}

#' Run CONOS
#' @param sce.list list of SingleCellExperiment objects for multi-omic datasets
#' @param integrate_features selected features to perform CCA on
#' @param reference reference dataset for FindTransferAnchors
#' @param query query dataset for FindTransferAnchors
#'
#' @return list of integration output (see details)
#'
#' @details Function outputs list of
#' 1) model:
#' 2) misc
#'
#' @import conos
#' @import Seurat
#'
#' @export
run_conos <- function(seurat.list, integrate_features, reference="RNA", query="ATAC"){
  data.processed <- seurat.list
  VariableFeatures(data.processed[[reference]]) <- integrate_features
  VariableFeatures(data.processed[[query]]) <- integrate_features
  data.processed <- map(data.processed, ~ ScaleData(.x) %>% RunPCA(dims=1:30))
  l.con <- Conos$new(data.processed,n.cores=30)
  l.con$buildGraph(k=15,k.self=5,k.self.weigh=0.01,ncomps=30,n.odgenes=5e3,space='PCA')

  l.con$findCommunities(resolution=1.5)
  l.con$embedGraph(alpha=1/2)
  return(list(model=l.con, input=seurat.list))
}

### Utils ###

#' Utility function for labelTransfer with Liger
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

#' Get cell type prediction from benchmark output
#' @param seu.int list of Seurat objects of label transfer output
#' @param int.name label for integration method
#' @param id.col column of predicted labels in obj metadata (default: "predicted.id")
#' @param score.col column of prediction score in obj metadata (default: "score")
#' @param filter_score numeric indicating the minimum prediction score to keep assigned label (default: 0)
#'
#' @return returns data.frame of label predictions and prediction scores
#'
#' @import dplyr
#' @import stringr
#'
#' @export
getPredictedLabels <- function(seu.int, int.name, id.col="predicted.id", score.col="score", filter_score=0, query="ATAC"){
  pred.df <- seu.int[[query]]@meta.data[,c(id.col, score.col), drop=F]
  colnames(pred.df) <- c('predicted.id', "score")
  pred.df <- pred.df %>%
    rownames_to_column("cell") %>%
    mutate(predicted.id = ifelse(score < filter_score, NA, as.character(predicted.id))) %>%
    mutate(method = int.name) %>%
    column_to_rownames("cell")
  rownames(pred.df) <- str_remove(rownames(pred.df), paste0("^", query,"_"))
  # colnames(pred.df) <- c(str_c("predicted.id", "_", int.name), str_c("score", "_", int.name))
  pred.df
}









# ----------------- Deprecated functions ----------------- #

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

  # integrate_features <- integrate_features[which(integrate_features %in% rownames(seurat.list[[reference]]) &
                                                   # integrate_features %in% rownames(seurat.list[[query]]))]
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
  # liger.obj@norm.data <- map(sce.list, ~ assay(.x, "logcounts"))
  # liger.obj@raw.data[[query]] <- Matrix(logcounts(sce.list[[query]]))
  liger.obj <- normalize(liger.obj)
  liger.obj@norm.data[[query]] <- Matrix(logcounts(sce.list[[query]]))
  ## Select genes
  integrate_features <- integrate_features[which(integrate_features %in% rownames(liger.obj@norm.data[[reference]]) &
                                                   integrate_features %in% rownames(liger.obj@norm.data[[query]]))]
  liger.obj@var.genes <- integrate_features
  ## Scale data
  liger.obj <- scaleNotCenter(liger.obj, remove.missing = F)
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
