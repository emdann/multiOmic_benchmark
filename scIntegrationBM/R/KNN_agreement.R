### KNN agreement score for label transfer ###

library(Seurat)
# library(tidyverse)
library(parallel)
library(purrr)

#' Get NN list from original data reduction
#' @param query.seu Seurat object of query dataset
#' @param k numeric indicating number of nearest neighbors to build graph
#' @param reduction.name character indicating reduction to use for NN graph (must be an element of `query.seu@reductions`)
#'
#' @return list of nearest neighbors per cell
#'
getQueryKNNgraph <- function(query.seu, k, reduction.name){
  query.seu <- FindNeighbors(query.seu, reduction = reduction.name, k.param = k)
  nn.list <- getNNlist(query.seu)
  nn.list
}


#' Calculate KNN agreement score with NN list
#' @param nn.list named list of k nearest neighbors for each cell
#' @param pred.labels named vector of predicted labels per cell
#'
#' @return vector of KNN agreement score for each cell
#'
#' @import purrr
#'
#' @export
KNN_score <- function(nn.list, pred.labels ){
                      # pred.id.col = "predicted.id"
                      # ){
  k <- unique(map_dbl(nn.list, length))
  nn.list <- nn.list[which(!is.na(pred.labels))]
  pred.labels <- pred.labels[which(!is.na(pred.labels))]

  if (length(k) > 1) {
    stop("KNN list is malformed. Different k for different cells.")
  }
  knn.score <- imap_dbl(nn.list, ~ sum(pred.labels[.x] == pred.labels[.y], na.rm = TRUE)/k)
  # data.frame(knn.score) %>% rownames_to_column("cell")
  knn.score
  }

#' KS-test for true KNN ECDF against ECDF of null distribution
#' @param nn.list list of nearest neighbors
#' @param pred.labels data.frame of predicted labels
#' @param n_perm number of permutations
#'
#' @importFrom parallel detectCores
#'
#' @export
KNN_purity_test <- function(nn.list, pred.labels, n_perm = 100, n_cores=detectCores()){
  true <- KNN_score(nn.list, pred.labels)
  null <- null_KNN_score(nn.list, pred.labels, n_perm = n_perm, n_cores=n_cores)
  ks <- ks.test(x = true, y=null, alternative = "less")
  list(fracKNN = true, null = null, KNN_purity=ks$statistic[1], p_val=format.pval(ks$p.value))
}

#' Calculate KNN agreement score with NN list
#' @param query.seu Seurat object of query dataset
#' @param pred.labels named vector of predicted labels per cell
#'
#' @return distribution of KNN purity per differen values of K
#'
#' @import purrr
#'
#' @export
KNN_purity_dist <- function(query.seu, pred.labels, reduction.name){
  ## Select Ks to range from 1 to 10% of total cells in dataset
  tot.cells <- ncol(query.seu)
  k.vector <- tot.cells*seq(0.01,0.1, by=0.005)
  ## Calculate KNN purity for each K
  knn_purity.df <- data.frame(k=k.vector, knn_purity=NA)
  for (k in k.vector) {
    nn.list <- getQueryKNNgraph(query.seu, k=k, reduction.name = reduction.name)
    KNN_purity <- KNN_purity_test(nn.list, pred.labels)
    df[df$k==k, "knn_purity"] <- KNN_purity$KNN_purity
  }
  knn_purity.df
}

#' Get list of NN for each cell from NN graph in seurat object
#' @param obj Seurat object (with graph) or neighborhood graph
#' @param is.seurat logical indicating whether `obj` is a Seurat object (default: TRUE)
#'
#' @return list of NN for each cell
#'
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

#' Calculate null KNN score for predicted labels
#' @details Implements a permutation test to avoid that a high KNN score is assigned to a label transfer method that just
#' gives the same label to everything
#'
#' @importFrom parallel mclapply
#'
null_KNN_score <- function(nn.list, pred.labels, n_perm = 10, n_cores=10){
  null <- mclapply(X = list(1:n_perm), function(x) KNN_score(nn.list, setNames(sample(pred.labels), nm = names(pred.labels))), mc.cores = n_cores)
  Reduce( c, null)
}
