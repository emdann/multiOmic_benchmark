##########################
### Data preprocessing ###
##########################

library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(Matrix)


#' Make list of SingleCellExperiment objects of assays to integrate
makeSCElist <- function(matrix.list){
  sce.list <- map(matrix.list, ~ SingleCellExperiment(assays=list(counts=.x)))
  return(sce.list)
}

#' Filter low quality cells
#' (empty cells, dead cells)
filterCells <- function(sce, nCounts=1000, nFeatures=600, fracMito=0.10){
  non.empty.cells <- assay(sce, "counts") %>% colSums() %>% {which(. > nCounts)} 
  highQ.cells <- diff(assay(sce, "counts")@p) %>% {which(. > nFeatures)}
  sce <- sce[,intersect(highQ.cells, non.empty.cells)]
  
  mito.genes <- rownames(assay(sce, "counts")) %>% {which(str_detect(., pattern = "MT"))}
  frac.mito.genes <- colSums(assay(sce, "counts")[mito.genes,])/(colSums(assay(sce, "counts")))
  # hist(frac.mito.genes, breaks=100)
  sce <- sce[,which(frac.mito.genes < fracMito)]
  return(sce)
  }

#' Normalize per cell
normalizePerCell <- function(sce){
  counts <- assay(sce, "counts")
  cpm <- apply(counts, 2, function(x) (x/sum(x))*1e-06)
  cpm(sce) <- as(object = cpm, Class = "dgCMatrix")
  return(sce)
}

## From SC for pedestrians
colVars_spm <- function( spm ) {
  stopifnot( is( spm, "dgCMatrix" ) )
  ans <- sapply( seq.int(spm@Dim[2]), function(j) {
    mean <- sum( spm@x[ (spm@p[j]+1):spm@p[j+1] ] ) / spm@Dim[1]
    sum( ( spm@x[ (spm@p[j]+1):spm@p[j+1] ] - mean )^2 ) +
      mean^2 * ( spm@Dim[1] - ( spm@p[j+1] - spm@p[j] ) ) } ) / ( spm@Dim[1] - 1 )
  names(ans) <- spm@Dimnames[[2]]
  ans
}

rowVars_spm <- function( spm ) {
  colVars_spm( t(spm) )
}

# counts <- assay(sce.list$RNA, "counts")
# counts <- counts[ rowSums(counts) > 0,  ]
# size_fact <- colSums(counts)
# nrm_counts <- t( t(counts) / colSums(counts) )
# 
# gene_means <- rowMeans( counts )
# gene_vars <- rowVars_spm( counts )
# 
# # plot( gene_means, gene_vars / gene_means,
# #       log = "xy", cex = .3, col = adjustcolor("black", alpha=.3), 
# #       xlab = "mean", ylab = "variance / mean" )
# 
# poisson_vmr <- mean( 1 / colSums( nrm_counts ) )
# 
# informative_genes <- names(which( 
#   gene_vars / gene_means  >  2 * poisson_vmr ))


# calculateQCMetrics(sce, 
#                    feature_controls = list(mito = mito.genes),
#                    cell_controls = list(empty = empty.cells, damaged = lowQ.cells))





















