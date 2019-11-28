### Preprocess Thymus dataset
library(magrittr)
source("~/multiOmic_benchmark/preprocess/preprocess.R")

thymus.rna_1 <- Read10X(data.dir = "~/my_data/Human_colon_16S7985396/")
thymus.rna_2 <- Read10X(data.dir = "~/my_data/Human_colon_16S7985397/")
colnames(thymus.rna_1) <- str_c(colnames(thymus.rna_1), "_1")
colnames(thymus.rna_2) <- str_c(colnames(thymus.rna_2), "_2")
thymus.rna <- cbind(thymus.rna_1, thymus.rna_2)

rna.sce <- SingleCellExperiment(assays=list(counts=thymus.rna))
# sce.list <- makeSCElist(list(RNA=thymus.rna, ATAC=thymus.atac.act))
rna.sce <- filterCells(rna.sce)
# sce.list$ATAC <- filterCells(sce.list$ATAC, fracMito = 1)

rna.sce <- normalizePerCell(rna.sce)

logcounts(rna.sce) <- log1p(cpm(rna.sce))

## Compatibility w Seurat
rownames(rna.sce) %<>% str_replace_all("_", "-")

## Add cell type annotation
annotation.df <- read.csv("~/my_data/F74_RNA_obs_v3.csv")
# annotation.df <- annotation.df %>%
#   dplyr::mutate(cell=str_remove(as.character(X), "F74_1_") %>% str_c(ifelse(batch==0,'_1', "_2"))) 

coldata <- annotation.df[,c("X", "anno_v3")] %>%
  dplyr::rename(annotation = anno_v3) %>%
  column_to_rownames("X") 
colData(rna.sce) <- DataFrame(coldata[colnames(rna.sce),, drop=F])

## Save SingleCellExperiment list
saveRDS(object = rna.sce, file = "~/my_data/F74_RNA_seurat_processed.RDS")
