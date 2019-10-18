#########################
### Feature selection ###
#########################

HVG_Seurat <- function(sce, nfeatures=2000){
  HVGs <- VariableFeatures(FindVariableFeatures(as.Seurat(sce), nfeatures = nfeatures))
}



