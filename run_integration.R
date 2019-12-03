#### RUN INTEGRATION ####

library(argparse)
library(stringr)

parser <- ArgumentParser()
parser$add_argument("sce_reference", type="character",
                    help = "Path to .RDS file of SingleCellExperiment object of reference dataset")
parser$add_argument("sce_query", type="character",
                    help = "Path to .RDS file of SingleCellExperiment object of query dataset")
# parser$add_argument("int_method", type="character",
#                     help = "Integration method to use (one of CCA, liger)")
parser$add_argument("integration_features", type='character',
                    help="Path to txt file of common features (e.g. genes) selected for integration")
parser$add_argument("study_id", type="character",
                    help="name of study")
parser$add_argument("--reference", default="RNA", type="character",
                    help="name of reference dataset (default: RNA)")
parser$add_argument("--query", default="ATAC", type="character",
                    help="name of query dataset (default: ATAC)")
parser$add_argument("--outdir", default="~/models2/", type="character",
                    help="Path to directory where to save integration outputs")
args <- parser$parse_args()

source("~/multiOmic_benchmark/integrateBenchmark.R")

sce.ref.path <- args$sce_reference
sce.query.path <- args$sce_query
reference <- args$reference
query <- args$query
integrate_features.path <- args$integration_features
study.id <- args$study_id
outdir <- args$outdir

feats.id <- str_remove_all(integrate_features.path, ".+/intFeatures_|_.+")

## Read input datasets
sce.ref <- readRDS(sce.ref.path)
sce.query <- readRDS(sce.query.path)
sce.list <- list(sce.ref, sce.query) %>% setNames(nm=c(reference, query))

## Read integration features
integrate_features <- scan(integrate_features.path, what="")

## Run models
message(" --- Running models ---", appendLF = T)

cca.model <- run_SeuratCCA(sce.list, integrate_features)
liger.model <- run_liger(sce.list, integrate_features)
conos.model <- run_conos(sce.list, integrate_features)

# cca.model <- readRDS("~/models2/modelCCA__unionHVGnHCGFeatures.RDS")
# liger.model <- readRDS("~/models2/modelLiger__unionHVGnHCGFeatures.RDS")
# conos.model <- readRDS("~/models2/modelConos__unionHVGnHCGFeatures.RDS")

message(" --- Saving models ---", appendLF = T)
model_outfiles <- 
  map(c("CCA", "Liger", "Conos"), ~ str_c("model", .x, "_",study.id,"_",feats.id,"Features.RDS"))
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
  map(c("CCA", "Liger", "Conos"), ~ str_c("labelTransfer", .x, "_",study.id,"_",feats.id,"Features.RDS"))
map2(list(CCA=seu.cca, Liger=seu.liger, Conos=seu.conos), labeltransfer_outfiles, ~ saveRDS(.x, file = str_c(outdir, .y)))
