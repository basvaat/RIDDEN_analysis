import pandas as pd
import numpy as np
from cmapPy.pandasGEXpress.parse import parse
from scipy.stats import spearmanr
from src.config import LINCS_DATA_DIR

geneinfo = pd.read_table(f'{LINCS_DATA_DIR}geneinfo_beta.txt')

# landmark genes
landmark_genes = geneinfo[geneinfo['feature_space']=='landmark']
landmark_genes = landmark_genes.set_index('gene_id', drop = True)['gene_symbol']
landmark_genes.index = landmark_genes.index.astype(str)

lm_map = landmark_genes.copy()
lm_map.index = landmark_genes.index.astype(float)
lm_map = lm_map.to_dict()

# inferred genes
inf_genes = geneinfo[(geneinfo['feature_space']=='landmark') | (geneinfo['feature_space']=='best inferred')]
inf_genes = inf_genes.set_index('gene_id', drop = True)['gene_symbol']
inf_genes.index = inf_genes.index.astype(str)

inf_map = inf_genes.copy()
inf_map.index = inf_genes.index.astype(float)
inf_map = inf_map.to_dict()
# print('Number of inferred genes: ', len(inf_genes))

pert_types = ['lig', 'oe', 'sh', 'xpr', 'cp']

def calc_MODZ(data):
    """calculates MODZ weights based on the original CMAP/L1000 study
    use only lm genes for MODZ calculation!"""
    if data.shape[1]==1: #1 column
        weights = np.array([1.0])
    elif data.shape[1]==2: # 2 columns
        weights = np.array([[0.5], [0.5]])
    else:
        CM = spearmanr(data)[0] # correlations (r) array with shape col x col
        fil = CM<0
        CM[fil] = 0.01 # correlation less than 0 = 0.01
        weights = np.sum(CM, 1)-1 # sum by rows minus the correlation woth itself (r=1)
        weights = weights / np.sum(weights)
        weights = weights.reshape((-1, 1)) # reshape to column,each value in separate array [[1][2][3]]
    return weights

def parse_data(pert_type, genes_filename, mapping):

    signature = parse(f'data/filtered_lincs_data/signatures_{genes_filename}_{pert_type}_raw_liana.gctx').data_df
    signature.index = signature.index.astype('float').map(mapping)
    sig_info = pd.read_csv(f'data/filtered_lincs_meta/filtered_{pert_type}_info_of_receptor_ligand_pert.csv', index_col = 0)

    return (signature, sig_info)

def calculate_consensus(pert_type, signature, sig_info):
    print('Read in signature info')
    
    perts = sig_info.cmap_name.unique()
    cells = sig_info.cell_iname.unique()
    pert_signatures = {}
    i=0
    print('Calculate consensus of perturbations and cells')
    for pert in perts:
        # print('pert ', pert, end = ': ')
        if i%100 == 0:
            print(i, end = ',')
        cell_signatures = {}
        for cell in cells:
            sample_ids = list(sig_info[(sig_info['cmap_name'] == pert) & (sig_info['cell_iname'] == cell)].sig_id)
            data = signature.loc[:, sample_ids]
            if data.shape[1] > 0:
                weights = calc_MODZ(data)
                cell_signatures[cell] = pd.DataFrame(np.dot(data, weights), index=data.index, columns = [cell])    
                
        pert_signatures[pert] = cell_signatures
        i = i+1
    print()
    return pert_signatures

def from_dict_create_signatures(pert_signatures, signature):
    print('Concatenate signatures to dataframe')
    all_signatures = {}
    for pert in list(pert_signatures.keys()):
        pert_signature = pert_signatures[pert]

        cols = [pert +'_'+ cell for cell in list(pert_signatures[pert].keys())]
        cell_signature = pd.DataFrame(index = signature.index, columns = cols)
        for cell in pert_signature: # cell (keys of pert_signature)
            cell_signature.loc[signature.index, pert+'_'+cell] = pert_signatures[pert][cell].iloc[:, 0]
        all_signatures[pert] = cell_signature.T
    all_df = pd.concat(all_signatures.values(), join = 'inner')
    return all_df

def save_consensus_signatures(pert_type, genes_filename, mapping):

    signature, sig_info = parse_data(pert_type, genes_filename, mapping)
    pert_signatures = calculate_consensus(pert_type, signature, sig_info)
    all_df = from_dict_create_signatures(pert_signatures, signature)
    print(f'Save result to data/lincs_consensus/{pert_type}_pert_cell_liana.csv')
    all_df.to_csv(f'data/lincs_consensus/{genes_filename}_{pert_type}_pert_cell_liana.csv')
    print('Done')

# create consensus signatures for landmark genes
# pert_types = ['lig', 'sh', 'xpr', 'oe', 'cp']
# for i in pert_types:
#     print('Perturbation: '+ i)
#     save_consensus_signatures(i, genes_filename = 'lm', mapping = lm_map) # or if inferred genes: 'inferred' and inf_map
#     print('----- DONE -----')

# create consensus signatures for inferred genes
pert_types = ['lig', 'sh', 'xpr', 'oe']
for i in pert_types:
    print('Perturbation: '+ i)
    save_consensus_signatures(i, genes_filename = 'inf', mapping = inf_map) # or if inferred genes: 'inferred' and inf_map
    print('----- DONE -----')