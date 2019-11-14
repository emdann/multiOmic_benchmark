#### RUN INTEGRATION ####

library(argparse)
library(stringr)

parser <- ArgumentParser()
parser$add_argument("sce_list_path", type="character",
                    help = "Path to .RDS file of SingleCellExperiment objects list to integrate")
# parser$add_argument("int_method", type="character",
#                     help = "Integration method to use (one of CCA, liger)")
parser$add_argument("--feature_selection_method", default='union_hvg', type='character',
                    help="How to select features to use for integration? (reference_hvg: HVGs from reference dataset, union_hvg: union of HVGs from reference and query)")
parser$add_argument("--n_features", default=2000, type='double',
                    help="Number of highly variable genes from the reference dataset to use for integration (default=2000)")
parser$add_argument("--reference", default="RNA", type="character",
                    help="name of reference dataset (default: RNA)")
parser$add_argument("--query", default="ATAC", type="character",
                    help="name of query dataset (default: ATAC)")
args <- parser$parse_args()

source("~/multiOmic_benchmark/integrateBenchmark.R")

sce.list.path <- args$sce_list_path
reference <- args$reference
query <- args$query
feature.selection <- args$feature_selection_method
n_features <- args$n_features

sce.list <- readRDS(sce.list.path)
# int_output <- run_integration(sce.list, method, n_features, reference=reference, query=query)

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


outdir <- "~/models/"

## Run models
message(" --- Running models ---", appendLF = T)
cca.model <- run_SeuratCCA(sce.list, integrate_features)
liger.model <- run_liger(sce.list, integrate_features)
conos.model <- run_conos(sce.list, integrate_features)

message(" --- Saving models ---", appendLF = T)
model_outfiles <- 
  map(c("CCA", "Liger", "Conos"), ~ str_c("model", .x, "_",feature.selection,"_",str_remove(sce.list.path, ".+/")))
map2(list(CCA=cca.model, Liger=liger.model, Conos=conos.model), model_outfiles, ~ saveRDS(.x, file = str_c(outdir, .y)))

## Transfer labels
message(" --- Label transfer ---", appendLF = T)
seu.cca <- labelTransfer_seuratCCA(seurat.list = cca.model$input,
                                         transfer.anchors = cca.model$model)
seu.liger <- labelTransfer_liger(liger.obj = liger.model$model,
                                       sce.list = liger.model$input, 
                                       k = 50)
seu.conos <- labelTransfer_conos(conos.out = conos.model$model,sce.list = conos.model$input)

message(" --- Save label transfer ---", appendLF = T)
labeltransfer_outfiles <- 
  map(c("CCA", "Liger", "Conos"), ~ str_c("labelTransfer", .x, "_",feature.selection,"_", str_remove(sce.list.path, ".+/")))
map2(list(CCA=seu.cca, Liger=seu.liger, Conos=seu.conos), labeltransfer_outfiles, ~ saveRDS(.x, file = str_c(outdir, .y)))

## Save features 
message(" --- Save features ---", appendLF = T)
featfile <- str_c("intFeatures","_",feature.selection,"_", n_features, "_", str_remove(str_remove(sce.list.path, ".+/"), ".RDS"), ".txt")
write(integrate_features, featfile)
