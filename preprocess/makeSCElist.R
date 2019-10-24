### Preprocess Thymus dataset
library(magrittr)
source("./multiOmic_benchmark/preprocess.R")

thymus.rna_1 <- Read10X(data.dir = "my_data/Human_colon_16S7985396/")
thymus.rna_2 <- Read10X(data.dir = "my_data/Human_colon_16S7985397/")
colnames(thymus.rna_1) <- str_c(colnames(thymus.rna_1), "_1")
colnames(thymus.rna_2) <- str_c(colnames(thymus.rna_2), "_2")
thymus.rna <- cbind(thymus.rna_1, thymus.rna_2)

thymus.atac.act <- readRDS("./my_data/cellranger-atac110_count_30439_WSSS8038360_GRCh38-1_1_0.geneActivity.RDS")

sce.list <- makeSCElist(list(RNA=thymus.rna, ATAC=thymus.atac.act))
sce.list$RNA <- filterCells(sce.list$RNA)
sce.list$ATAC <- filterCells(sce.list$ATAC, fracMito = 1)

sce.list <- map(sce.list, ~ normalizePerCell(.x))

logcounts(sce.list$RNA) <- log1p(cpm(sce.list$RNA))
logcounts(sce.list$ATAC) <- log1p(cpm(sce.list$ATAC))

## Compatibility w Seurat
rownames(sce.list$RNA) %<>% str_replace_all("_", "-")
rownames(sce.list$ATAC) %<>% str_replace_all("_", "-")

## Save SingleCellExperiment list
saveRDS(object = sce.list, file = "my_data/integrated_thymus/F74_SCElist_20191017.RDS")
