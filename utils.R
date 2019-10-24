#############
### Utils ###
#############

library(RColorBrewer)

#' Convert output of integration method to Seurat object
#' @param intOut output of integration
#' @param int.features features used for integration (HVGs in reference dataset)
#' @param reference name of reference assay
intOutput2seurat <- function(intOut, int.features, reference="RNA"){
  int.seu <- as.Seurat(intOut[[paste0("integrated.", reference)]], counts=NULL, data="normcounts")
  int.seu@assays$RNA@var.features <- int.features
  return(int.seu)
}


### Plotting utils ###

brewer_palette_4_values <- function(vec, palette, seed=42){
  set.seed(seed)
  pal=suppressWarnings(colorRampPalette(brewer.pal(9, palette))(length(vec)))
  names(pal) <- sample(vec)
  return(pal)
}
