#########################
### Feature selection ###
#########################

library(Seurat)

select_highlyVariable <- function(sce, nfeatures=2000){
  HVGs <- VariableFeatures(FindVariableFeatures(as.Seurat(sce, counts = "counts",data = "logcounts"), nfeatures = nfeatures
                                                selection.method="mvp", dispersion.cutoff=c(0.5, 10), mean.cutoff=c(0.0125, 6)
                                                ))
  HVGs
  }


select_highlyCovered <- function(sce, frac_cells=0.1){
  HCFs <-
    names(which(Matrix::rowSums(logcounts(sce)) > frac_cells*ncol(sce)))
  HCFs
    }


