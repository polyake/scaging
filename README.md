# Human single-cell transcriptomics aging clocks

We developed cell-type-specific single-cell transcriptomic clocks predicting the age of human blood cells. Here, we provide the code and pipeline for the application of the clocks to new data. The clocks and other necessary files can be found on Zenodo: 10.5281/zenodo.10405106

## Code and data for the application of the single-cell clocks
### Requirements
* Python 3.8.17
* requirements.txt

### Data
The single-cell clocks and all relevant data can be found at https://zenodo.org/records/10405106

### Pipeline
1. Single-cell gene expression data has to be in a *./Data* folder in the form of an AnnData object, named *local.h5ad* <br />
It is assumed that the raw gene expressions are in the main (X) layer of the object and the age of individuals are contained in the *development_stage* column of the .obs layer. Small modification of the code is needed otherwise.
2. Single-cell clocks has to be copied to a *./models* folder
3. Data for missing value imputation has to be copied to a *./data_for_imputation* folder
4. Cell types clocks are applied to can be specified in the *apply_clocks_general.py* file
5. Apply the single-cell clocks to new data: <br />
*python apply_clocks_general.py*
