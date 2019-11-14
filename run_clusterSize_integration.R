library(Seurat)
library(glue)
source("~/multiOmic_benchmark/integrateBenchmark.R")

model.cca <- readRDS("~/models/modelCCA_union_hvg_PBMC_SCElist_20191105.RDS")
model.liger <- readRDS("~/models/modelLiger_union_hvg_PBMC_SCElist_20191105.RDS")
model.conos <- readRDS("~/models/modelConos_union_hvg_PBMC_SCElist_20191105.RDS")

integrate_features <- scan("~/models/intFeatures_union_hvg_2000_F74_SCElist_20191101.txt", what = "")

reclusterReference <- function(RNA.input.seu, res){
  RNA.input.seu <- FindClusters(RNA.input.seu, resolution = res, algorithm = 4) ## clustering w leiden algorithm
  annotation <- RNA.input.seu@meta.data$seurat_clusters
  return(annotation)
  }

getPredictedLabels <- function(seu.int, int.name, id.col="predicted.id", score.col="score"){
  pred.df <- seu.int$ATAC@meta.data[,c(id.col, score.col), drop=F] 
  rownames(pred.df) <- str_remove(rownames(pred.df), "^ATAC_")
  colnames(pred.df) <- c(str_c("predicted.id"), str_c("score"))
  pred.df
  }

labelTransfer_clusters <- function(model, res, annotation, method, reference="RNA", query="ATAC"){
  if (method=="CCA") {
    labelTransfer <- labelTransfer_seuratCCA
    score.col="prediction.score.max"
  } else if (method=="conos") {
    labelTransfer <- labelTransfer_conos
    score.col = "score"
  } else if (method=="liger") {
    labelTransfer <- labelTransfer_liger
    score.col = "score"
  } else {
    stop("Method must be one of CCA, liger, conos")
  }
  new.input <- model$input
  new.input[[reference]]$annotation <- annotation
  lt.seu <- labelTransfer(model$model, new.input)
  pred.labels.df <- getPredictedLabels(seu.int = lt.seu, int.name = glue("{method}_res{res}"), score.col = score.col)  
  return(pred.labels.df)
  }

## Convert predicted labels to long data.frame
long.pred.labels <- function(pred.labels, method, resolutions){
  map(seq_along(resolutions), ~ pred.labels[[.x]] %>% 
        rownames_to_column('cell') %>%
        mutate(res=resolutions[.x], method=method)) %>%
    purrr::reduce(bind_rows)
}


## Find NN graph on reference dataset
RNA.input <- model.cca$input$RNA
if ( class(RNA.input)!="Seurat" ) { 
  RNA.input.seu <- as.Seurat(RNA.input)  
} else {
  RNA.input.seu <- RNA.input
}

RNA.input.seu <- FindVariableFeatures(RNA.input.seu)
RNA.input.seu <- ScaleData(RNA.input.seu)
RNA.input.seu <- RunPCA(RNA.input.seu)
RNA.input.seu <- FindNeighbors(RNA.input.seu, reduction = "pca", dims=1:30)

## Clustering w different resolutions
resolutions <- seq(0.5,1.5, by = 0.1)
annotation_list <- map(resolutions, ~ reclusterReference(RNA.input.seu, res = .x))

pred.labels.liger <- map2(annotation_list, resolutions, ~ 
                            labelTransfer_clusters(model.liger, annotation=.x, res= .y, method = "liger"))
pred.labels.cca <- map2(annotation_list, resolutions, ~ 
                            labelTransfer_clusters(model.cca, annotation=.x, res= .y, method = "CCA"))
pred.labels.conos <- map2(annotation_list, resolutions, ~ 
                            labelTransfer_clusters(model.conos, annotation=.x, res= .y, method = "conos"))

cluster.size.df <- imap(list(CCA=pred.labels.cca, conos=pred.labels.conos, liger=pred.labels.liger), ~ 
      long.pred.labels(.x, method = .y, resolutions = resolutions)) %>%
      purrr::reduce(bind_rows)

saveRDS(cluster.size.df, "models/clusterSize_labelTransfer_PBMC.RDS")
saveRDS(annotation_list, "~/models/clusterSize_trueAnnotations_PBMC.RDS")














