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

#' Get list of NN for each cell from NN graph in seurat object
getNNlist <- function(obj, is.seurat=TRUE){
  if (is.seurat) {
    nn.graph <- obj@graphs[[1]]
  } else {
      nn.graph <- obj
    }
  nn.graph %>%
    apply(1, function(x) names(which(x==1))) %>%
    asplit(2)
}


eucl_distance <- function(p,q){
  sqrt(sum((p - q)^2))
}


### Plotting utils ###

brewer_palette_4_values <- function(vec, palette, seed=42){
  set.seed(seed)
  pal=suppressWarnings(colorRampPalette(brewer.pal(9, palette))(length(vec)))
  # names(pal) <- sample(vec)
  return(pal)
}

DimPlotCluster <- function(atac.seu, annotation_col, cluster, label=NULL){
  if (is.null(label)) {
    label <- cluster
  }
  highlight = which(atac.seu@meta.data[,annotation_col]==cluster)
  DimPlot(atac.seu, reduction = "umap.snap",cells.highlight = highlight, cols.highlight = "red", pt.size = 0.02, sizes.highlight = 0.1) +
    guides(color="none") +
    ggtitle(label = label)
}
