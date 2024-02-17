
from tqdm import tqdm
import pandas as pd
import numpy as np
from cmapPy.pandasGEXpress.parse import parse
from cmapPy.pandasGEXpress.write_gctx import write
from glob import glob
from scipy.stats import spearmanr

DATA_DIR = # data path
LINCS_DATA_DIR = # LINCS path
geneinfo = pd.read_table(LINCS_DATA_DIR+'geneinfo_beta.txt')
inf_genes = geneinfo[(geneinfo['feature_space']=='landmark') | (geneinfo['feature_space']=='best inferred')]
inf_genes = inf_genes.set_index('gene_id', drop = True)['gene_symbol']
inf_genes.index = inf_genes.index.astype(str)

inf_map = inf_genes.copy()
inf_map.index = inf_genes.index.astype(float)
inf_map = inf_map.to_dict()

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

def calculate_consensus(pert_type, signature, sig_info):
    print('Read in signature info')
    cells = sig_info.cell_iname.unique()
    i=0
    print('Calculate consensus of perturbations and cells')
    cell_signatures = {}
    for cell in cells:
        sample_ids = list(sig_info[(sig_info['cell_iname'] == cell)].sig_id)
        data = signature.loc[:, sample_ids]
        if data.shape[1] > 0:
            weights = calc_MODZ(data)
            cell_signatures[cell] = pd.DataFrame(np.dot(data, weights), index=data.index, columns = [cell])                
        i = i+1
    print()
    return cell_signatures


def calc_pert_cons(pert_type, filename, genes, genes_filename, mapping):
    # read in filtered LINCS meta
    filepath = f'{DATA_DIR}/filtered_lincs_meta/filtered_{pert_type}_info_of_receptor_ligand_pert.csv'
    print(filepath)
    print(f'Filtering LINCS {pert_type} perturbation data for receptors and ligands')
    siginf = pd.read_csv(filepath, index_col = 0)
    siginf = siginf[(siginf['ligand'] == 1) | (siginf['receptor'] == 1)]
    print()
    perts = siginf.cmap_name.unique()
    for pert in perts:
        print('Drug consensus calculation:', pert)
        siginf_drug = siginf[(siginf['cmap_name'] == pert)]
        print(len(siginf_drug.cell_iname.unique()))
        print('Parsing from gctx...')
        lname = glob(LINCS_DATA_DIR+'level5_beta_'+filename+'*.gctx')[0]
        print(lname)
        signature = parse(lname, cid=siginf_drug.sig_id,rid=genes.index).data_df
        signature.index = signature.index.astype('float').map(mapping)
        print('Calculating consensus...')
        cell_signatures = calculate_consensus(pert_type, signature, siginf_drug)
        cell_signatures_df = pd.concat(cell_signatures.values(), axis=1)
        print('Save to file...')
        cell_signatures_df.to_csv(f'{DATA_DIR}/lincs_consensus/inferred_genes_signatures/signatures_'+genes_filename+'_'+pert_type+'_'+pert+'_consensus.csv')
        print('')


pert_types = {'lig':'trt_misc', 'oe':'trt_oe', 'sh':'trt_sh', 'xpr':'trt_xpr', 'cp':'trt_cp'}
calc_pert_cons('cp', filename = pert_types['cp'], genes = inf_genes, genes_filename = 'inf', mapping = inf_map)