### Preprocess Thymus dataset
library(magrittr)
source("./multiOmic_benchmark/preprocess/preprocess.R")

thymus.rna_1 <- Read10X(data.dir = "my_data/Human_colon_16S7985396/")
thymus.rna_2 <- Read10X(data.dir = "my_data/Human_colon_16S7985397/")
colnames(thymus.rna_1) <- str_c(colnames(thymus.rna_1), "_1")
colnames(thymus.rna_2) <- str_c(colnames(thymus.rna_2), "_2")
thymus.rna <- cbind(thymus.rna_1, thymus.rna_2)

# thymus.atac.act <- readRDS("./my_data/cellranger-atac110_count_30439_WSSS8038360_GRCh38-1_1_0.geneActivity.RDS")
x.sp <- readRDS("~/my_data/cellranger-atac110_count_30439_WSSS8038360_GRCh38-1_1_0.snapATAC.RDS")
thymus.atac.act <- t(x.sp@gmat)
thymus.atac.act <- thymus.atac.act[names(which(table(rownames(thymus.atac.act)) == 1)),] # Remove duplicate rows

sce.list <- makeSCElist(list(RNA=thymus.rna, ATAC=thymus.atac.act))
sce.list$RNA <- filterCells(sce.list$RNA)
sce.list$ATAC <- filterCells(sce.list$ATAC, fracMito = 1)

sce.list <- map(sce.list, ~ normalizePerCell(.x))

logcounts(sce.list$RNA) <- log1p(cpm(sce.list$RNA))
logcounts(sce.list$ATAC) <- log1p(cpm(sce.list$ATAC))

## Compatibility w Seurat
rownames(sce.list$RNA) %<>% str_replace_all("_", "-")
rownames(sce.list$ATAC) %<>% str_replace_all("_", "-")

## Add cell type annotation to RNA data
annotation.df <- read.csv("~/my_data/F74_RNA_obs_v2.csv")
annotation.df <- annotation.df %>%
  dplyr::mutate(cell=str_remove(as.character(X), "F74_1_") %>% str_c(ifelse(batch==0,'_1', "_2"))) 

coldata <- annotation.df[,c("cell", "anno_v2")] %>%
  dplyr::rename(annotation = anno_v2) %>%
  column_to_rownames("cell") 
colData(sce.list$RNA) <- DataFrame(coldata[colnames(sce.list$RNA),, drop=F])

## Save SingleCellExperiment list
saveRDS(object = sce.list, file = "my_data/integrated_thymus/F74_SCElist_20191113.RDS")
