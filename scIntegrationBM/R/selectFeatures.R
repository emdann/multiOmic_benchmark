#########################
### Feature selection ###
#########################

#' Select highly variable features
#' Wrapper around Seurat function
#' @param sce SingleCellExperiment object (must have counts and logcounts for log-normalized data)
#' 
select_highlyVariable <- function(sce){
  HVGs <- VariableFeatures(FindVariableFeatures(as.Seurat(sce, counts = "counts",data = "logcounts"), nfeatures = nfeatures,
                                                selection.method="mvp", dispersion.cutoff=c(0.5, 10), mean.cutoff=c(0.0125, 3)
                                                ))
  HVGs
  }

#' Select highly covered features
#' 
#' @param sce SingleCellExperiment object (must have counts and logcounts for log-normalized data)
#' @param frac_cells fraction of cells in which a feature has to be different from 0 to be included
#' 
select_highlyCovered <- function(sce, frac_cells=0.1){
  HCFs <-
    names(which(Matrix::rowSums(logcounts(sce)) > frac_cells*ncol(sce)))
  HCFs
    }


