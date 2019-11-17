#########################
### Feature selection ###
#########################

library(Seurat)

HVG_Seurat <- function(sce, nfeatures=2000){
  HVGs <- VariableFeatures(FindVariableFeatures(as.Seurat(sce), nfeatures = nfeatures, 
                                                selection.method="mvp", dispersion.cutoff=c(0.5, 100000), mean.cutoff=c(0.0125, 6)))
}





