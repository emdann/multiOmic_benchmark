### Make ATAC counts per genes (as suggested by LIGER tutorial, but custom code)
library(liger)

genes.bc <- read.table(file = "~/my_data/cellranger-atac110_count_30439_WSSS8038360_GRCh38-1_1_0.genes_bc.bed", header = FALSE)
promoters.bc <- read.table(file = "~/my_data/cellranger-atac110_count_30439_WSSS8038360_GRCh38-1_1_0.promoters_bc.bed", header = FALSE)

## Convert barcode column into list
genes.bc <- mutate(genes.bc, V2=strsplit(as.character(V2), ','))
promoters.bc <- mutate(promoters.bc, V2=strsplit(as.character(V2), ','))

## Make gene x cell count matrix
genes_count_mat <- unnest(genes.bc, V2) %>%
  group_by(V1, V2) %>%
  summarise(count=n()) %>%
  spread(V2, count) %>%
  column_to_rownames("V1") %>%
  as.matrix() %>%
  {ifelse(is.na(.), 0,.)}

promoters_count_mat <- unnest(promoters.bc, V2) %>%
  group_by(V1, V2) %>%
  summarise(count=n()) %>%
  spread(V2, count) %>%
  column_to_rownames("V1") %>%
  as.matrix() %>%
  {ifelse(is.na(.), 0,.)}

rownames(promoters_count_mat) %in% rownames(genes_count_mat)
promoters_count_mat[rownames(genes_count_mat),colnames(genes_count_mat)] + genes_count_mat

#The makeFeatureMatrix function requires a list of barcodes as a vector
barcodes = names(atac_clusts)
gene.counts <- makeFeatureMatrix(genes.bc, barcodes)
promoter.counts <- makeFeatureMatrix(promoters.bc, barcodes)

gene.counts <- gene.counts[order(rownames(gene.counts)),]
promoter.counts <- promoter.counts[order(rownames(promoter.counts)),]
pbmc.atac <- gene.counts + promoter.counts