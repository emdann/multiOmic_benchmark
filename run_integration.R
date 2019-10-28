#### RUN INTEGRATION ####

library(argparse)
library(stringr)

parser <- ArgumentParser()
parser$add_argument("sce_list_path", type="character",
                    help = "Path to .RDS file of SingleCellExperiment objects list to integrate")
parser$add_argument("int_method", type="character",
                    help = "Integration method to use (one of CCA, liger)")
parser$add_argument("--n_features", default=4000, type='double',
                    help="Number of highly variable genes from the reference dataset to use for integration")
parser$add_argument("--reference", default="RNA", type="character",
                    help="name of reference dataset")
parser$add_argument("--query", default="ATAC", type="character",
                    help="name of query dataset")
args <- parser$parse_args()

source("~/multiOmic_benchmark/integrateBenchmark.R")

sce.list.path <- args$sce_list_path
reference <- args$reference
query <- args$query
method <- args$int_method
n_features <- args$n_features

sce.list <- readRDS(sce.list.path)
int_output <- run_integration(sce.list, method, n_features, reference=reference, query=query)

outdir <- "~/models/"
outfile <- str_c("integrate", method, "_", str_remove(sce.list.path, ".+/"))
featfile <- str_c("intFeatures","_", n_features, "_", str_remove(str_remove(sce.list.path, ".+/"), ".RDS"), ".txt")
saveRDS(int_output[1:2], str_c(outdir, outfile))
saveRDS(int_output[[3]], str_c(outdir, featfile))