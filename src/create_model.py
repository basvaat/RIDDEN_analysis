import pandas as pd
import numpy as np
import statsmodels.api as sm
from src.config import *
from src.model_creation import *

def read_in_ligand_receptor_info():
    print('Read in ligand receptor matrix')
    lr_matrix_path = LIG_REC_MATRIX
    return pd.read_csv(lr_matrix_path, index_col =0)

def read_in_data(pert_type, genes_filename):
    print('Read in signature and metadata')
    metadata = pd.read_csv(f'data/filtered_lincs_meta/filtered_{pert_type}_info_of_receptor_ligand_pert.csv', index_col = 0)
    signature = pd.read_csv(f'data/lincs_consensus/{genes_filename}_{pert_type}_pert_cell_liana.csv', index_col = 0)
    # NOTE : this is necessary because the signautre filtering can differ at the beginning (if publish: filter correct signature and delete these 2 lines)
    metadata['id'] = metadata['cmap_name'] + '_' + metadata['cell_iname']
    signature = signature[signature.index.isin(list(set(metadata['id'])))]
    metadata = metadata[metadata['id'].isin(list(set(signature.index)))]
    return metadata, signature

def read_in_filtered_compound_info():
    print('Read in filtered compound info')
    return pd.read_csv('data/filtered_lincs_meta/filtered_coumpound_info_to_receptor_perturbation_signatures_signed.csv', index_col =0)

def create_compound_model_matrix(metadata, signature, compound_info):
    print('Create compound model matix')
    compound_info = compound_info[['cmap_name', 'target', 'sign']].drop_duplicates()
    cp_pert_cell_meta = metadata[['cmap_name', 'cell_iname']].drop_duplicates().reset_index(drop = True)
    cp_targets = cp_pert_cell_meta.merge(compound_info[['cmap_name', 'target', 'sign']], on = 'cmap_name')
    cp_targets['pert'] = cp_targets['cmap_name'] + '_' + cp_targets['cell_iname']
    cp_targets.sign = cp_targets.sign.astype('int')
    assert(len(cp_targets.drop_duplicates(subset = ['pert', 'target'])) == len(cp_targets))
    model_matrix = cp_targets.pivot_table(index='pert', columns='target', values='sign', aggfunc='sum').fillna(0)
    return model_matrix

def create_receptor_model_matrix(metadata, signature):
    print('Creating receptor model matrix')
    receptors = metadata[metadata['receptor'] == 1].cmap_name.unique()
    if len(receptors) == 0:
        print('There is no receptor in dataset.')
        return pd.DataFrame() # return empty dataframe
    select_receptor_indices = list(filter(lambda x: x.split('_')[0] in receptors, signature.index))
    receptor_signatures = signature.loc[select_receptor_indices]
    model_matrix_receptor = receptor_signatures.reset_index()['index'].str.split('_', expand=True)[0].str.get_dummies().set_index(receptor_signatures.index)
    return model_matrix_receptor

def create_ligand_model_matrix(metadata, signature, lr_matrix):
    print('Creating ligand model matrix')
    # create ligand model matrx
    # select from dataframe ligands
    ligands = metadata[metadata['ligand'] == 1].cmap_name.unique()
    if len(ligands) == 0:
        print('There is no ligand in dataset.')
        return pd.DataFrame()
    select_ligand_indices = list(filter(lambda x: x.split('_')[0] in ligands, signature.index))
    ligand_signatures = signature.loc[select_ligand_indices]
    # select ligands from lr_matrix that are included in signatures
    lr_common_ligands = lr_matrix.loc[list(set(lr_matrix.index) & set(ligands))]
    # since I filter aout ligands, maybe there will be receptors, that are modulated by ligands that was filtered out -> remove receptors with 0 interactions
    cols = lr_common_ligands.columns[abs(lr_common_ligands).sum() != 0]
    lr_common_ligands = lr_common_ligands[cols]
    ligands = ligand_signatures.reset_index()['index'].str.split('_', expand = True)[0]
    ligands.index = ligand_signatures.index
    ligands.name = 'ligand'
    # create model matrix by merging filtered lr_matrix to index of ligand perturbations of signature
    model_matrix_ligand = pd.merge(ligands, lr_common_ligands, left_on = 'ligand', right_index = True, how = 'left').drop(columns = ['ligand'])
    return model_matrix_ligand

def create_model_matrix_from_ligand_and_receptor_model_matrices(signature, model_matrix_receptor, model_matrix_ligand):
    print('Concatenate model_matrixes to create a whole')
    # concatenate receptor and ligand model matrix and sort as original dataframe
    # model_matrices can be empty DataFrames, so they concatenate these to the other 
    # if there is no receptor or ligand in a dataset (eg. ligand perturbations - no receptor)
    model_matrix = pd.concat([model_matrix_receptor, model_matrix_ligand]).fillna(0)
    model_matrix = model_matrix.loc[signature.index]
    model_matrix = model_matrix.astype('int64')
    return model_matrix

def set_sign_modelmatrix(model_matrix, pert_type):
    return model_matrix*pert_sign[pert_type]


def save_model_matrix(model_matrix, pert_type):
    path = f'data/design_matrices/{pert_type}_pert_binary_liana.csv'
    print('Save binary matrix to '+  path)
    model_matrix.to_csv(path)
    print('Done')

def save_coeffitient_matrix(coeff_matrix, pert_type):
    path = f'data/coefficient_matrix/{pert_type}_pert_coef_liana.csv'
    print('Save coefficient matrix to '+  path)
    coeff_matrix.to_csv(path)
    print('Done')

def pipeline_calculating_coeffitient_matrix(pert_type, genes_filename):
    print('Calculating and saving coefficient matrix of LINCS signatures pert type: ', pert_type)
    lr_matrix = read_in_ligand_receptor_info()
    metadata, signature = read_in_data(pert_type, genes_filename)
    # if perturbation is compound create model matrix different way
    # metadata of compound perturbation contains only compound names 
    # and other dataframe contains corresponding receptor targets 
    if pert_type == 'cp':
        compound_info = read_in_filtered_compound_info()
        model_matrix = create_compound_model_matrix(metadata, signature, compound_info)
        signature = signature.loc[model_matrix.index]
    else:
        ligand_model_matrix = create_ligand_model_matrix(metadata, signature, lr_matrix)
        receptor_model_matrix = create_receptor_model_matrix(metadata, signature)
        model_matrix =  create_model_matrix_from_ligand_and_receptor_model_matrices(signature, receptor_model_matrix, ligand_model_matrix)

    assert all(model_matrix.apply(lambda s: pd.to_numeric(s, errors='coerce').notnull().all()))
    assert all(signature.index == model_matrix.index)

    model_matrix = set_sign_modelmatrix(model_matrix, pert_type)
    save_model_matrix(model_matrix, pert_type)
    # fit linear model using signature and model matrix (gene expression signature = model matrix x beta coeffitiens)
    coefficiens = fit_linear_model_and_get_coefficients_by_receptor(y = signature, X = model_matrix)
    save_coeffitient_matrix(coefficiens, pert_type)

# Calculate coefficient matrix for each perturbation type
pert_sign = {'xpr':-1, 'sh':-1, 'oe':1, 'lig':1, 'cp':1}
pert_types = ['sh', 'oe', 'lig', 'cp', 'xpr'] 
for pert_type in pert_types:
    pipeline_calculating_coeffitient_matrix(pert_type, 'lm')

