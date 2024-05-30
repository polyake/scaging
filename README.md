# Supplementary data for the manuscript titled Profiling the transcriptomic age of single-cells in humans

## Abstract
Although aging clocks predicting the age of individual organisms have been extensively studied, the age of individual cells remained largely unexplored. Most recently single-cell omics clocks were developed for the mouse, however, profiling the age of human cells is still lacking. To fill this gap, here we use available scRNA-seq data of 1,058,909 blood cells of 508 healthy, human donors (between 19 and 75 years), for developing single-cell transcriptomic clocks and predicting the age of human blood cells. By the application of the proposed cell-type-specific single-cell clocks, our main observations are that (i) naive CD8-positive T cells are transcriptomically younger than CD8-positive memory T cells; (ii) the transcriptomic age of classical monocyte cells is decreased in COVID-19 subjects; and (iii) the transcriptomic age of cells decreased followed by an increase during embryonic development. In summary, here we demonstrate that single-cell transcriptomic clocks are useful tools to investigate the aging process at the single-cell level.

## Code and data for the application of the single-cell clocks
### Requirements
* Python 3.8.17
* requirements.txt

### Pipeline
1. Single-cell gene expression data has to be in a *./Data* folder in the form of an AnnData object, named *local.h5ad* <br />
It is assumed that the raw gene expressions are in the main (X) layer of the object and the age of individuals are contained in the *development_stage* column of the .obs layer. Small modification of the code is needed otherwise.
2. Single-cell clocks has to be copied to a *./models* folder
3. Data for missing value imputation has to be copied to a *./data_for_imputation* folder
4. Cell types clocks are applied to can be specified in the *apply_clocks_general.py* file
5. Apply the single-cell clocks to new data: <br />
*python apply_clocks_general.py*
