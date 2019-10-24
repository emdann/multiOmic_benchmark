#############
### Utils ###
#############

#' Convert output of integration method to Seurat object
#' @param intOut output of integration
#' @param int.features features used for integration (HVGs in reference dataset)
#' @param reference name of reference assay
intOutput2seurat <- function(intOut, int.features, reference="RNA"){
  int.seu <- as.Seurat(intOut[[paste0("integrated.", reference)]], counts=NULL, data="normcounts")
  int.seu@assays$RNA@var.features <- int.features
  return(int.seu)
  }