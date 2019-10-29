Benchmarking methods for alignment of scRNA-sea and scATAC-seq data
================

## Dataset details

Incl. stats about dataset quality (median no. of fragments per cell,
perc. of fragments mapping to peaks, median no. of fragments in peaks
per cell)

## Data preprocessing

### ATAC-seq

raw 10x sequencing reads were preprocessed to call peaks with Cellranger
(v …). For aligment with scRNA-seq data we need to reduce ATAC peaks to
gene-level features. We do this as in the Seurat pipeline, summing all
counts in peaks within gene bodies + 2kb upstream.

### RNA-seq

Both datasets where preprocessed with std steps: removal of empty/lowQ
cells, normalization per cell coverage,

#### Feature selection

  - Highly-variable genes in RNA
  - genes that are expressed only in a few cells

## Tested integration methods

  - **Seurat CCA**:
      - select K?
  - **LIGER:**
      - K factors
      - Imputation strategy: projecting ATAC cells on NMF factor space.
        Testing if this is a valid imputation strategy using scRNA-seq
        data only
  - **SnapATAC pipeline:** it just wraps CCA alignment
  - **scGen:** requires cell type annotation also on the ATAC dataset
  - **totalVI:** not applicable as it assumes matching between cells
    (cite-Seq and RNA-seq from the same single-cells)
  - **BBKNN:** (how do you do an imputation step?)

## Uniform output for all methods

  - Impute transcriptomic data

## Metrics for comparison of integration models

Problems: finding optimal distance metric for gene accessibility and
expression (correlation doesn’t work, too sparse). (Or more general: how
to relate features from different datasets) Nearest neighbors in PCA?
How to make a nearest neighbor graph between modalities - MNN cosine
normalizes each batch and then calculates distances

  - PCA on integrated space and find genes with high eigen values
  - How to denoise the ATAC data??

<!-- end list -->

1)  Robustness to different methods of feature selection: HVGs in the
    RNA
2)  Robustness to different fractions of cells in ATAC dataset
3)  Leave-one-out approach for imputed data
4)  **Fraction of unassigned cells** (but how to distinguish unassigned
    and badly assigned?)
5)  **Joint clustering: purity of cell type annotation inside a cluster,
    mixing within the same cluster between different technologies**
6)  **Robustness to parameter picking (e.g. no. of factors)**
7)  Agreement meric defined by Welch et al. 2019: compare KNN graph of
    single datasets with KNN graph in integrated space, then calculate
    how many of each cell’s NNs in the single dataset graph are also NNs
    in the integrated graph. Welch et al. compare NN graphs built with
    different factor models to compare CCA and LIGER performance. I
    build for all methods KNN graphs from the PCA projection of imputed
    values.
8)  **Expression not at random of markers after integration:** is the
    structure that is clear in the RNA maintained after embedding w the
    ATAC? From idea of JP collaborator
9)  **pySCENIC on RNA only or integration w ATAC seq**

## Ideas for less “agnostic” integration

1)  Select only a certain lineage of cells (e.g. that you can align in
    pseudotime)
2)  Annotation of cell types also in ATAC-seq data (e.g. to use scGen)
3)  Considering enhancer accessibility (matching them to genes??)

## Does adding the ATAC information improve the inference of gene regulatory networks?

  - running SCENIC on full integrated data VS just on RNA

## Other random things

  - How does ATAC improve the RNA? Can we detect cell
    clusters/populations that are highly homogeneous in the RNA but
    display significant variability at the accessibility level? But what
    is variability in super noisy ATAC seq?
  - Studying pioneering TFs: temporal relationship between TF expression
    and accessibility (are they really opening chromatin?) Could be done
    on the time series data
