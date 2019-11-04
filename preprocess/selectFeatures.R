#########################
### Feature selection ###
#########################

library(Seurat)

HVG_Seurat <- function(sce, nfeatures=2000){
  HVGs <- VariableFeatures(FindVariableFeatures(as.Seurat(sce), nfeatures = nfeatures))
}





