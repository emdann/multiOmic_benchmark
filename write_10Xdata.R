### RDS to AnnData (PBMC data)
library(DropletUtils)

## Read data
rna_clusts = readRDS("~/10X_data/rna_cluster_assignments.RDS")
atac_clusts = readRDS("~/10X_data/atac_cluster_assignments.RDS")
pbmc.atac <- readRDS('~/10X_data/pbmc.atac.expression.mat.RDS')
pbmc.rna <- readRDS('~/10X_data/pbmc.rna.expression.mat.RDS')

## Select common genes 
my.genes <- intersect(rownames(pbmc.atac), rownames(pbmc.rna))

## Concatenate matrices
int.pbmc <- cbind(pbmc.rna[my.genes,], pbmc.atac[my.genes,])

cell.ids <- colnames(int.pbmc)

ngenes <- nrow(int.pbmc)
gene.symb <- rownames(int.pbmc)

# Writing this to file
write10xCounts(path = "./10X_data/intPBMC_10x.hdf5",type = "HDF5", int.pbmc,
               gene.symbol=gene.symb, barcodes=cell.ids)

as.data.frame(rna_clusts) %>%
  write.csv("./10X_data/rna_cluster_assignments.csv")

as.data.frame(atac_clusts) %>%
  write.csv("./10X_data/atac_cluster_assignments.csv")
