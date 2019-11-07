Benchmarking methods for alignment of scRNA-sea and scATAC-seq data
================

## Dataset details

Incl. stats about dataset quality (median no. of fragments per cell,
perc. of fragments mapping to peaks, median no. of fragments in peaks
per cell)

## Data preprocessing

### ATAC-seq

ATAC-seq data is inherently sparse and noisy, more so than RNA-seq data.
In addition, methods for alignment with scRNA-seq data mostly require to
reduce accessibility signal to gene-level features. Different papers use
different strategies to preprocess ATAC-seq data and featurization
(benchmarked by Chen et al.
([2019](#ref-chenAssessmentComputationalMethods2019a))):

  - **Seurat pipeline**: raw 10x sequencing reads were preprocessed to
    call peaks with Cellranger (v …). Then all counts in peaks are
    summed up all counts in peaks within gene bodies + 2kb upstream
    (Welch et al. [2019](#ref-welchSingleCellMultiomicIntegration2019a);
    Stuart et al.
    [2019](#ref-stuartComprehensiveIntegrationSingleCell2019a)).
  - **SNAP-ATAC pipeline:** the genome is binned, then a cell x bin
    binary matrix is constructed. From this binary matrix a gene-level
    matrix is made (“Fast and Accurate Clustering of Single Cell
    Epigenomes Reveals Cis-Regulatory Elements in Rare Cell Types,”
    [n.d.](#ref-FastAccurateClusteringa)).
  - **PCA based on bulk:** some studies (Buenrostro et al.
    [2018](#ref-buenrostroIntegratedSingleCellAnalysis2018)) find
    principal components of differentiation/cell type from bulk ATAC-seq
    datasets and then project the single-cells on the same space.
    Arguably, this will lead to miss out on accessibility dynamics of
    rare cell populations.
  - **LSI method:** this procedure starts by making a bin x cell matrix.
    Then normalization and rescaling are done using term
    frequency-inverse document frequency transformation (something from
    text mining), then SVD is performed to obtain a PC x cell matrix.
    This is used to cluster cells and peaks are called on the clusters.
    Then the peak x cell matrix is used to do PCA again (I am getting
    dizzy).

### RNA-seq

Both datasets where preprocessed with std steps: removal of empty/lowQ
cells, normalization per cell coverage,

#### Feature selection

  - Highly-variable genes in RNA
  - genes that are expressed only in a few cells
  - Also ATAC specific markers are informative, contain signal from
    known marker genes in thymic development. Probably beneficial for
    alignment to take the union and not only the HVGs from RNA.

## Tested integration methods

| Method         | Reference                                                                                                  | Included in benchmark |                Reason for Excluding                 |
| -------------- | ---------------------------------------------------------------------------------------------------------- | :-------------------: | :-------------------------------------------------: |
| Seurat CCA     | Stuart et al. ([2019](#ref-stuartComprehensiveIntegrationSingleCell2019a))                                 |          Yes          |                          /                          |
| LIGER          | Welch et al. ([2019](#ref-welchSingleCellMultiomicIntegration2019a))                                       |          Yes          |                          /                          |
| Conos          | (<span class="citeproc-not-found" data-reference-id="barkasJointAnalysisHeterogeneous2019">**???**</span>) |          Yes          |                          /                          |
| scGen          | Lotfollahi, Wolf, and Theis ([2019](#ref-lotfollahiScGenPredictsSinglecell2019))                           |          No           |   Requires cell type annotation in both datasets    |
| totalVI        | Gayoso et al. ([2019](#ref-gayosoJointModelRNA2019))                                                       |          No           | Requires multi-omic data from the same single-cells |
| BBKNN          | Polański et al. ([n.d.](#ref-polanskiBBKNNFastBatch))                                                      |          No           |            Bad alignment during testing             |
| Cusanovich2018 | Cusanovich et al. ([2018](#ref-cusanovichSingleCellAtlasVivo2018a))                                        |          No           |                  Code unavailable                   |

### Label transfer

One of the key tasks for integration methods is to be able to transfer
cell type annotations learnt from a reference to a query dataset. This
is especially useful if the query is a scATAC-seq dataset, where calling
of cell types based on prior knownledge on marker genes is often not
possible. Different models are adapted to transfer discrete cell state
labels derived from gene expression to cells measured with scATAC-seq.

##### Seurat CCA

Identified anchor pairs are weighted based on the query cell local
neighboorhood (the k nearest anchors) and the anchor score. The obtained
reference cells x query cells weight matrix is then multiplied by the
matrix of annotation x reference cells, to generate a query cell x
annotation matrix. This returns a prediction score for each class for
every cell in the query dataset, ranging from 0 to 1 (Stuart et al.
[2019](#ref-stuartComprehensiveIntegrationSingleCell2019a)).

##### LIGER

While the authors do not describe a method for transferring discrete
labels, I adapted their strategy for feature imputation. I build a
cross-dataset KNN graph in the aligned factor space, then I assign each
query cell to the most abundant label between the k nearest neighbors in
the reference dataset (k=30). The prediction score for label \(l\) is
computed as the fraction of nearest neighbors that have the predicted
label. \[
score = \frac{count(l)}{k}
\] \#\#\#\#\# Conos Label transfer is treated as a general problem of
information propagation between vertices of the common graph (detailed
in
(<span class="citeproc-not-found" data-reference-id="barkasJointAnalysisHeterogeneous2019">**???**</span>)).
The label score is the label probability updating during the diffusion
process.

<!-- ## Uniform output for all methods  -->

<!-- - Impute transcriptomic data -->

<!-- - **Transfer labels from RNA-seq to ATAC-seq** -->

## Metrics for comparison of integration models

Problems: finding optimal distance metric for gene accessibility and
expression (correlation doesn’t work, too sparse). (Or more general: how
to relate features from different datasets) Nearest neighbors in PCA?
How to make a nearest neighbor graph between modalities - MNN cosine
normalizes each batch and then calculates distances

  - PCA on integrated space and find genes with high eigen values
  - How to denoise the ATAC data??

<!-- end list -->

2)  Robustness to different fractions of cells in ATAC dataset
3)  Leave-one-out approach for imputed data
4)  **Fraction of unassigned cells** (but how to distinguish unassigned
    and badly assigned?): prediction score from label transfer
5)  **Joint clustering: purity of cell type annotation inside a cluster,
    mixing within the same cluster between different technologies**
6)  **Robustness to parameter picking (e.g. no. of factors)**
7)  Agreement metric defined by Welch et al. 2019: compare KNN graph of
    single datasets with KNN graph in integrated space, then calculate
    how many of each cell’s NNs in the single dataset graph are also NNs
    in the integrated graph. Welch et al. compare NN graphs built with
    different factor models to compare CCA and LIGER performance. I
    build for all methods KNN graphs from the PCA projection of imputed
    values.
8)  **Expression not at random of markers after integration:** is the
    structure that is clear in the RNA maintained after embedding w the
    ATAC? From idea of JP collaborator
      - Find marker genes from RNA only
      - Measure non-random expression in RNA only
9)  **pySCENIC on RNA only or integration w ATAC seq**

## Ideas for less “agnostic” integration

1)  Select only a certain lineage of cells (e.g. that you can align in
    pseudotime)
2)  Annotation of cell types also in ATAC-seq data (e.g. to use scGen)
3)  Considering enhancer accessibility (matching them to genes??)

## Biological interpretation of integration

In Cusanovich et al. 2018 they identify differentially accessible sites
& cluster specific accessibility patterns.

## Other random things

  - How does ATAC improve the RNA? Can we detect cell
    clusters/populations that are highly homogeneous in the RNA but
    display significant variability at the accessibility level? But what
    is variability in super noisy ATAC seq?
  - Studying pioneering TFs: temporal relationship between TF expression
    and accessibility (are they really opening chromatin?) Could be done
    on the time series data Svensson and Pachter
    ([2019](#ref-svenssonInterpretableFactorModels2019))
  - Clonality of epigenetic modifications ATAC + TCR tracing

## Bibliography

<div id="refs" class="references">

<div id="ref-buenrostroIntegratedSingleCellAnalysis2018">

Buenrostro, Jason D., M. Ryan Corces, Caleb A. Lareau, Beijing Wu,
Alicia N. Schep, Martin J. Aryee, Ravindra Majeti, Howard Y. Chang, and
William J. Greenleaf. 2018. “Integrated Single-Cell Analysis Maps the
Continuous Regulatory Landscape of Human Hematopoietic Differentiation.”
*Cell* 173 (6): 1535–1548.e16.
<https://doi.org/10.1016/j.cell.2018.03.074>.

</div>

<div id="ref-chenAssessmentComputationalMethods2019a">

Chen, Huidong, Caleb Lareau, Tommaso Andreani, Michael E. Vinyard, Sara
P. Garcia, Kendell Clement, Miguel A. Andrade-Navarro, Jason D.
Buenrostro, and Luca Pinello. 2019. “Assessment of Computational Methods
for the Analysis of Single-Cell ATAC-Seq Data.” *bioRxiv*, August,
739011. <https://doi.org/10.1101/739011>.

</div>

<div id="ref-cusanovichSingleCellAtlasVivo2018a">

Cusanovich, Darren A., Andrew J. Hill, Delasa Aghamirzaie, Riza M. Daza,
Hannah A. Pliner, Joel B. Berletch, Galina N. Filippova, et al. 2018. “A
Single-Cell Atlas of in Vivo Mammalian Chromatin Accessibility.” *Cell*
174 (5): 1309–1324.e18. <https://doi.org/10.1016/j.cell.2018.06.052>.

</div>

<div id="ref-FastAccurateClusteringa">

“Fast and Accurate Clustering of Single Cell Epigenomes Reveals
Cis-Regulatory Elements in Rare Cell Types.” n.d., 41.

</div>

<div id="ref-gayosoJointModelRNA2019">

Gayoso, Adam, Romain Lopez, Zoë Steier, Jeffrey Regier, Aaron Streets,
and Nir Yosef. 2019. “A Joint Model of RNA Expression and Surface
Protein Abundance in Single Cells.” *bioRxiv*, October, 791947.
<https://doi.org/10.1101/791947>.

</div>

<div id="ref-lotfollahiScGenPredictsSinglecell2019">

Lotfollahi, Mohammad, F. Alexander Wolf, and Fabian J. Theis. 2019.
“scGen Predicts Single-Cell Perturbation Responses.” *Nature Methods*
16 (8): 715. <https://doi.org/10.1038/s41592-019-0494-8>.

</div>

<div id="ref-polanskiBBKNNFastBatch">

Polański, Krzysztof, Matthew D. Young, Zhichao Miao, Kerstin B. Meyer,
Sarah A. Teichmann, and Jong-Eun Park. n.d. “BBKNN: Fast Batch Alignment
of Single Cell Transcriptomes.” *Bioinformatics*. Accessed October 3,
2019. <https://doi.org/10.1093/bioinformatics/btz625>.

</div>

<div id="ref-stuartComprehensiveIntegrationSingleCell2019a">

Stuart, Tim, Andrew Butler, Paul Hoffman, Christoph Hafemeister,
Efthymia Papalexi, William M. Mauck, Yuhan Hao, Marlon Stoeckius, Peter
Smibert, and Rahul Satija. 2019. “Comprehensive Integration of
Single-Cell Data.” *Cell* 177 (7): 1888–1902.e21.
<https://doi.org/10.1016/j.cell.2019.05.031>.

</div>

<div id="ref-svenssonInterpretableFactorModels2019">

Svensson, Valentine, and Lior Pachter. 2019. “Interpretable Factor
Models of Single-Cell RNA-Seq via Variational Autoencoders.” Preprint.
Bioinformatics. <https://doi.org/10.1101/737601>.

</div>

<div id="ref-welchSingleCellMultiomicIntegration2019a">

Welch, Joshua D., Velina Kozareva, Ashley Ferreira, Charles Vanderburg,
Carly Martin, and Evan Z. Macosko. 2019. “Single-Cell Multi-Omic
Integration Compares and Contrasts Features of Brain Cell Identity.”
*Cell* 177 (7): 1873–1887.e17.
<https://doi.org/10.1016/j.cell.2019.05.006>.

</div>

</div>
