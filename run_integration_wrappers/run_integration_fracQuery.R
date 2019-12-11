#### RUN INTEGRATION / QUERY DATASET FRACTION ####

library(argparse)
library(stringr)

parser <- ArgumentParser()
parser$add_argument("sce_list_path", type="character",
                    help = "Path to .RDS file of SingleCellExperiment objects list to integrate")
parser$add_argument("int_method", type="character",
                    help = "Integration method to use (one of CCA, liger, Conos)")
parser$add_argument("percent_query", type="double",
                    help = "Fraction of cells in query dataset to keep")
parser$add_argument("seed", type="double",
                    help = "Random seed")
parser$add_argument("--feature_selection_method", default='union_hvg', type='character',
                    help="How to select features to use for integration? (reference_hvg: HVGs from reference dataset, union_hvg: union of HVGs from reference and query)")
parser$add_argument("--n_features", default=2000, type='double',
                    help="Number of highly variable genes from the reference dataset to use for integration (default=2000)")
parser$add_argument("--reference", default="RNA", type="character",
                    help="name of reference dataset")
parser$add_argument("--query", default="ATAC", type="character",
                    help="name of query dataset")
args <- parser$parse_args()

source("~/multiOmic_benchmark/integrateBenchmark.R")

sce.list.path <- args$sce_list_path
seed <- args$seed
reference <- args$reference
query <- args$query
feature.selection <- args$feature_selection_method
perc <- args$percent_query
n_features <- args$n_features
if (args$int_method=="CCA") {
  run <- run_SeuratCCA
} else if (args$int_method=="liger") {
  run <- run_liger
} else if (args$int_method=="conos") {
  run <- run_conos
} else {
  stop("Wrong integration method: choose one of CCA, liger, conos")
}

sce.list <- readRDS(sce.list.path)
# int_output <- run_integration(sce.list, method, n_features, reference=reference, query=query)

## Subset query dataset 
# perc <- 0.75
set.seed(seed)
keep_query <- sample(colnames(sce.list[[query]]), size = ncol(sce.list[[query]])*perc)
sce.list[[query]] <- sce.list[[query]][,keep_query]

## Feature selection
message(" --- Feature selection ---", appendLF = T)
if (feature.selection=="reference_hvg") {
  integrate_features <- HVG_Seurat(sce.list[[reference]], nfeatures = n_features) %>% {.[. %in% rownames(sce.list[[query]])]}
} else if (feature.selection == "union_hvg"){
  integrate_features_ref <- HVG_Seurat(sce.list[[reference]], nfeatures = n_features) %>% {.[. %in% rownames(sce.list[[query]])]}
  integrate_features_query <- HVG_Seurat(sce.list[[query]], nfeatures = n_features) %>% {.[. %in% rownames(sce.list[[reference]])]}
  integrate_features <- union(integrate_features_ref, integrate_features_query)
} else {
  stop("Invalid feature selection method: please specify one of 'reference_hvg' or 'union_hvg' ")
}

outdir <- "~/models/cell_frac/"

start_time <- Sys.time()

## Run models
message(" --- Running model ---", appendLF = T)
model <- run(sce.list, integrate_features)

# message(" --- Saving models ---", appendLF = T)
# model_outfiles <- 
#   map(c("CCA", "Liger", "Conos"), ~ str_c("model", .x, "_",feature.selection,"_",str_remove(sce.list.path, ".+/")))
# map2(list(CCA=cca.model, Liger=liger.model, Conos=conos.model), model_outfiles, ~ saveRDS(.x, file = str_c(outdir, .y)))

## Transfer labels
message(" --- Label transfer ---", appendLF = T)
if (args$int_method=="CCA") {
  seu <- labelTransfer_seuratCCA(seurat.list = model$input,
                                     transfer.anchors = model$model)
  colnames(seu$ATAC@meta.data) <- str_replace(colnames(seu$ATAC@meta.data), "prediction.score.max", "score")
} else if (args$int_method=="liger"){
  seu <- labelTransfer_liger(liger.obj = model$model,
                                   sce.list = model$input, 
                                   k = 50)
    } else if (args$int_method=="conos") {
  seu <- labelTransfer_conos(conos.out = model$model,sce.list = model$input)
  }

end_time <- Sys.time()

time_diff <- end_time - start_time

outtab <- seu$ATAC@meta.data[,c("predicted.id", "score")]
outtab["time"] <- time_diff

outfile <- str_c("fracQuery", args$int_method, "_frac", perc, "_seed", seed, "_", str_remove_all(sce.list.path, ".+/|.RDS"), ".csv")
write.csv(x = outtab, file = str_c(outdir, outfile))
# 
# message(" --- Save label transfer ---", appendLF = T)
# labeltransfer_outfiles <- 
#   map(c("CCA", "Liger", "Conos"), ~ str_c("labelTransfer", .x, "_",feature.selection,"_", str_remove(sce.list.path, ".+/")))
# map2(list(CCA=seu.cca, Liger=seu.liger, Conos=seu.conos), labeltransfer_outfiles, ~ saveRDS(.x, file = str_c(outdir, .y)))

