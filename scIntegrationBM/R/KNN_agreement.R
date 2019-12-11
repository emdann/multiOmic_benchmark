### KNN agreement score for label transfer ###

library(Seurat)
# library(tidyverse)
library(parallel)
library(purrr)

#' Calculate KNN agreement score with NN list
#' @param nn.list named list of k nearest neighbors for each cell
#' @param pred.labels named vector of predicted labels per cell
#' 
#' @return vector of KNN agreement score for each cell
#' 
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

#' Calculate null KNN score for predicted labels
#' @details A permutation test is implemented to avoid that a high KNN score is assigned to a label transfer method that just 
#' gives the same label to everything
#' 
#' 
null_KNN_score <- function(nn.list, pred.labels, n_perm = 10, n_cores=10){
  null <- mclapply(X = list(1:n_perm), function(x) KNN_score(nn.list, setNames(sample(pred.labels), nm = names(pred.labels))), mc.cores = n_cores)
  Reduce( c, null)
}

#' KS-test for true KNN ECDF against ECDF of null distribution
test.knn <- function(nn.list, pred.labels,  n_perm = 100, n_cores=detectCores()){
  true <- KNN_score(nn.list, pred.labels)
  null <- null_KNN_score(nn.list, pred.labels, n_perm = n_perm, n_cores=n_cores)
  ks <- ks.test(x = true, y=null, alternative = "less")
  list(KNN_score = true, null = null, D=ks$statistic[1], p.val=format.pval(ks$p.value))
}
