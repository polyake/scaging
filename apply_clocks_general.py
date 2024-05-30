print("It works!")
import pandas as pd
import scanpy as sc
import numpy as np
import os
import json
import anndata as ad

### reading the data that we want to apply the clock to
adata = ad.read('./data/local.h5ad', backed='r')


### this function reads data from a h5ad file only for the given cell type into a dataframe -- this is a memory saving way to read the needed data, because in this way we don't load the entire h5ad file into memory
### it assumes that raw expression counts are available (either through the main layer or the .raw layer)
### for more details see the comments inside the function
def get_df(data, cell_type):
    ###########
        # data: a loaded anndata file, but it does not have to be loaded into memory, loading with "backed='r'" option is enough
        # cell_type: cell type in a string format, that you would like to load into memory; it has to be character-wisely the same as the one in the "cell_type" column of the data.obs dataframe
    ###########
    
    k=data[data.obs["cell_type"]==cell_type]
    sh=(k.obs.shape[0],k.var.shape[0])
    cols=k.var_names
    indexes=k.obs_names
    
    ### If there is a .raw layer, then the next 3 rows should be commented out
    spcd=k.X.data
    spcind=k.X.indices
    spcptr=k.X.indptr
    ### with .raw layer, the next 3 rows should be used instead
    #spcd=k.raw.X.data
    #spcind=k.raw.X.indices
    #spcptr=k.raw.X.indptr
    
    sufnimtr=np.zeros(shape=sh,dtype=np.float32)
    for i in range(sh[0]):
        indmsk=spcind[spcptr[i]:spcptr[i+1]]
        data=spcd[spcptr[i]:spcptr[i+1]]
        sufnimtr[i][indmsk]=data
    df=pd.DataFrame(sufnimtr,index=indexes,columns=cols)
    print("base_df Done!")
    
    ### normalization step --> if the data is already normalized, the next 2 rows can be skipped
    df = df.div(df.sum(axis=1), axis=0)*10000
    df = df.applymap(np.log1p)
    
    ### adding age and donor id info to the dataframe
    ### here, age is given in the "development_stage" column in the following format: "45-year-old human stage --> we split this on the character "-" and keep the first item, i.e. the age as an integer number
    age = k.obs["development_stage"].str.split("-", expand=True)
    age_list = age[0].astype(int)
    donor_list = k.obs["donor_id"]
    df["age"] = age_list
    df["donor_id"] = donor_list
    print("df Done!")
    return df


### applying one clock to data from one cell type
### it is assumed that the name of the genes are the same in the training and application dataset, thus they can be mapped onto each other directly
### for more details see the comments inside the function
def apply_clock_celltype(data, cell_type, impute_data, model):
    ###########
        # data: dataframe containing the normalized expression counts and age and donor info --> the output of the get_df function
        # cell_type: cell type in a string format
        # impute_data: dataframe containing imputation values for genes in the training dataset --> it is used to impute missing values (expression value of genes that are not present in the application dataset)
        # model: clock to be applied, in a dataframe format containing the coefficients and the interception value
    ###########
    
    predictions = []
    ages = []
    donors = []
    cell_names = []
    
    ### taking the genes that are present in both the set of clock features and the application dataset
    genes_in_model = set(model.columns)&set(data.columns)
            
    ### keeping only the common genes in the application dataset + donor and age info
    cols_to_keep = list(genes_in_model)+["donor_id", "age"]
    df_celltype_age = data[cols_to_keep]
    
    ### coefficients and intercept of the clock
    model_coeff = model[model.columns[:-1]]
    model_intercept = model["intercept"]
    
    ### getting the set of genes that is present in the model, but is missing from the application dataset, and getting imputation values for these genes
    genes_rem = list(set(model_coeff.columns)-set(genes_in_model))
    data_for_imputation = impute_data[genes_rem]
    
    ### applying the clock to each cell, one-by-one
    for i in range(len(df_celltype_age)):
        cell_name = df_celltype_age.iloc[i].name
        ### missing value imputation
        row = pd.concat([df_celltype_age.iloc[i].to_frame().T.reset_index(drop=True), data_for_imputation], axis=1)
    
        donor = row["donor_id"]
        age = row["age"]
        row = row[list(set(row.columns)-set(["donor_id","age"]))]
        
        ### applying the clock --> dot product of the expression values and the model coefficients + interception
        prediction = row.dot(model_coeff.T) + model_intercept
    
        predictions.append(prediction[0][0])
        ages.append(age[0])
        donors.append(donor[0])
        cell_names.append(cell_name)
    
    print("Prediction Done!")
    
    ### output: list of predictions, ages, donor ids, cell names --> each element of the lists represents one cell
    return predictions, ages, donors, cell_names


### this is the main function
### the function "apply_clock_celltype" is applied with all 5 models of a given cell type specific clock (resulting from 5-fold cv)
def predict(data, cell_type, dataset_name):
    ###########
        # data: loaded anndata file; it will be the input of the function "get_df"
        # cell_type: cell type in a string format; here, it is assumed that the given cell types are present in the application dataset and the same names can be found in the name of the imputation and model files
        # dataset_name: technical argument for the name of the resulting file
    ###########
    
    ### reading imputation data and models
    ### these 2 rows can be changed according to the name of the files (e.g. if the cell type names are not the same in the training and application datasets)
    impute_data = pd.read_csv("./data_for_imputation/Impute_avg_"+cell_type+".csv")
    models_celltype = pd.read_csv("./models/"+cell_type+"_models5.csv")
    
    ### reading the expression data (application dataset) with the function "get_df"
    df_final = pd.DataFrame()
    df_base = get_df(data, cell_type)
    
    ### applying the 5 clocks one-by-one
    for j in range(len(models_celltype)):
        ### taking one clock from the dataframe of all 5 clocks
        current_model = models_celltype.iloc[j].to_frame().T.reset_index()
        current_model = current_model.drop("index",axis=1)
        
        ### taking only the non-zero coefficients from the clock
        m2 = (abs(current_model) > 0).any()
        a = m2.index[m2]
        model_to_use = current_model[list(a)]
        
        ### applying "apply_clock_celltype" with the current clock
        pred, ages, donors, cell_names = apply_clock_celltype(df_base, cell_type, impute_data, model_to_use)
        
        ### adding the predictions given by the current clock to the final dataframe, with donor, age and cell name info
        df = pd.DataFrame()
        df["donor"] = donors
        df["age"] = ages
        df["predicted_age"] = pred
        df["cell_name"] = cell_names
        df_final = pd.concat([df_final, df])
        
    ### writing the final datarframe into a csv file
    df_final.to_csv(dataset_name+"_predictions_selflognormed_"+cell_type+".csv", index=False)
    
### example for the applications
common_celltypes = ['CD14-low, CD16-positive monocyte',
 'CD16-negative, CD56-bright natural killer cell, human',
 'CD4-positive, alpha-beta T cell',
 'CD4-positive, alpha-beta cytotoxic T cell',
 'CD8-positive, alpha-beta T cell',
 'conventional dendritic cell',
 'dendritic cell',
 'effector memory CD4-positive, alpha-beta T cell',
 'erythrocyte',
 'gamma-delta T cell',
 'innate lymphoid cell',
 'memory B cell',
 'mucosal invariant T cell',
 'naive B cell',
 'naive thymus-derived CD8-positive, alpha-beta T cell',
 'natural killer cell',
 'plasmacytoid dendritic cell',
 'platelet',
 'regulatory T cell']
### running the prediction process --> here, "common_celltypes" is a list of cell types that we want to apply the speicific clocks to
for cell_type in common_celltypes:
    predict(adata, cell_type, "application")