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
  pal=suppressWarnings(colorRampPalette(brewer.pal(12, palette))(length(vec)))
  # names(pal) <- sample(vec)
  return(pal)
}

DimPlotCluster <- function(atac.seu, annotation_col, cluster, label=NULL, reduct="umap.snap"){
  if (is.null(label)) {
    label <- cluster
  }
  highlight = which(atac.seu@meta.data[,annotation_col]==cluster)
  DimPlot(atac.seu, reduction = reduct, cells.highlight = highlight, cols.highlight = "red", pt.size = 0.02, sizes.highlight = 0.1) +
    guides(color="none") +
    ggtitle(label = label) +
    theme(axis.ticks = element_blank(), axis.text = element_blank())
}


FeaturePlotCluster <- function(atac.seu, annotation_col, feature_col, cluster, label=NULL, reduct="umap.snap"){
  if (is.null(label)) {
    label <- cluster
  }
  atac.seu.copy <- atac.seu
  highlight = which(atac.seu.copy@meta.data[,annotation_col]==cluster)
  grey_cells = which(atac.seu.copy@meta.data[,annotation_col]!=cluster)
  atac.seu.copy@meta.data[grey_cells,feature_col] <- NA
  FeaturePlot(atac.seu.copy, reduction = reduct, feature=feature_col, order=FALSE, pt.size = 0.05, cells=highlight) +
    scale_color_viridis_c(na.value = "grey80") +
    # guides(color="none") +
    ggtitle(label = label) 
    # theme(axis.ticks = element_blank(), axis.text = element_blank())
}
