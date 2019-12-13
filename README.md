# Optimizing integration of scRNA-seq and scATAC-seq datasets

### Content

##### Analysis notebooks

- `/PBMC_benchmark/` comparison of integration methods on PBMC dataset from 10X genomics
  - `20191111_fracQuery_EDA.Rmd` contains method run time and robustness analysis
  - `20191111_labelTransferEDA_PBMC.Rmd` contains comparison of label transfer outcomes (visualization, prediction score analysis, KNN purity analysis) 
  - (`20191113_clusterSize_PBMC.Rmd`)[https://github.com/EmmaDann/multiOmic_benchmark/blob/master/PBMC_benchmark/20191113_clusterSize_PBMC.Rmd] contains comparison of dependency of label transfer on cluster size
- `/Thymus_benchmark/`



### Notebooks to reproduce report figures

- Fig.1A-B: `20191111_fracQuery.Rmd`
- Fig.1C-D: `20191106_labelTransferEDA_PBMC.Rmd`
- Fig.2A-B,D,F: `20191125_gmat_EDA.Rmd` 
- Fig.2C: `20191128_labelTransferEDA_thymus_v2_countgmat.Rmd`
- Fig.2E: `20191128_labelTransferEDA_thymus_v2_bgmat.Rmd`
- Fig.3A-G: `20191127_tcellTrajectory.Rmd`
