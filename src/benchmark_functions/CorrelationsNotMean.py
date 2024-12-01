from pandas import DataFrame, isna
from numpy import mean
import random


def general_correlation(genes_to_iterate: list, activities: DataFrame, input_data: DataFrame, method = 'pearson'):
    """Compute correlation of gene values and activity values, by iterating through the given list of gene names
    
    Parameters
    ----------
    genes_to_iterate : list
        list of gene names in activities DataFrame to iterate over
    activities : DataFrame
        estimated activities, rows: sampels, columns: receptors
    input_data:  DataFrame
        TPM or z-score transformed values of samples' gene expression
    method : str
        DataFrame.corr correlation method (default is pearsonr)
    
    Returns
    -------
    dict
       Correlation values of the activity and gene values across samples, keys: genes, values: r (correlation)
    """    
    correlations = {}
    for gene in genes_to_iterate:
        if gene not in input_data.columns:
            continue
        correlations[gene] = activities.loc[:, gene].astype('float').corr(input_data.loc[:, gene].astype('float'), method = method)
    return correlations

def activity_and_gene_correlation(input_data:DataFrame, activities:DataFrame): 
    """Correlation of receptor/ligand pathway activity and receptors/ligand TPM/zscore
    
    Parameters
    ----------
    input_data : DataFrame
        TPM or z-score transformed values of samples' gene expression
    activities:  DataFrame
        estimated receptor/ligand activities, rows: sampels, columns: receptors
  
    Returns
    -------
    dict
       Correlation values of the receptor/ligand pathway activity and receptors/ligand TPM/zscore, keys: receptors/ligands, values: r (correlation)
    """
    genes_of_receptors = list(set(input_data.columns) & set(activities.columns))
    correlations = general_correlation(genes_of_receptors, activities, input_data, method = 'pearson')
    return correlations


def ligand_activity_receptor_correlation(lig_rec: DataFrame, input_data: DataFrame, ligand_activities: DataFrame):
    """Correlation of ligand pathway activity (eg. CytoSig prediction) and mean of their receptor TPM/zscore 
    
    Parameters
    ----------
    lig_rec : DataFrame
        prior knowledge of ligand - receptor interactions
    input : DataFrame
        TPM or z-score transformed values of samples' gene expression
    ligand_activities:  DataFrame
        estimated ligand activities, rows: sampels, columns: receptors
    
    Returns
    -------
    dict
       Correlation values of the ligand activity and receptor TPM/zscore, keys: ligands, values: r (correlation)
    """
    corr_dict = {}

    ligand_receptor_corr_tg = {}
    # target : receptor, source: ligand
    for ligand in ligand_activities.columns:
        ligtarget_df = lig_rec[lig_rec['source_genesymbol'] == ligand]
        if len(ligtarget_df) == 0: 
            continue
        corr_list = []

        for idx in ligtarget_df.index:
            targetgene = ligtarget_df.target_genesymbol.loc[idx]

            if targetgene not in input_data.columns:
                continue
            if all(input_data[targetgene] == 0):
                continue    
            method = 'pearson'
            gene_values = input_data.loc[:, targetgene].astype('float')
            correlation = ligand_activities.loc[:, ligand].astype('float').corr(gene_values, method = method)
            ligand_receptor_corr_tg[ligand + '_'+ targetgene] = correlation
    return ligand_receptor_corr_tg


def receptor_activity_ligand_correlation(lig_rec: DataFrame, input_data: DataFrame, receptor_activities: DataFrame):
    """Correlation of ligand pathway activity (eg. CytoSig prediction) and mean of their receptor TPM/zscore 
    
    Parameters
    ----------
    lig_rec : DataFrame
        prior knowledge of ligand - receptor interactions
    input : DataFrame
        TPM or z-score transformed values of samples' gene expression
    receptor_activities:  DataFrame
        estimated ligand activities, rows: sampels, columns: receptors
    
    Returns
    -------
    dict
       Correlation values of the ligand activity and receptor TPM/zscore, keys: ligands, values: r (correlation)
    """
    corr_dict = {}

    receptor_ligand_corr_tg = {}
    # target : receptor, source: ligand
    for receptor in receptor_activities.columns:
        recsource_df = lig_rec[lig_rec['target_genesymbol'] == receptor]
        if len(recsource_df) == 0: 
            continue
        # corr_list = []

        for idx in recsource_df.index:
            sourcegene = recsource_df.source_genesymbol.loc[idx]

            if sourcegene not in input_data.columns:
                continue
            if all(input_data[sourcegene] == 0):
                continue    
            method = 'pearson'
            gene_values = input_data.loc[:, sourcegene].astype('float')
            correlation = receptor_activities.loc[:, receptor].astype('float').corr(gene_values, method = method)
            receptor_ligand_corr_tg[receptor + '_'+ sourcegene] = correlation
    return receptor_ligand_corr_tg


def receptor_activity_and_receptorxligand_correlation(lig_rec: DataFrame, input:DataFrame, receptor_activities:DataFrame, zscore = True): # return_no_target_list = False
    """Correlation of receptor activity and mean of ligand x its receptors TPM/ranks
    
    Parameters
    ----------
    lig_rec : DataFrame
        prior knowledge of ligand - receptor interactions
    input : DataFrame
        TPM or z-score transformed values of samples' gene expression
    receptor_activities:  DataFrame
        estimated receptor activities, rows: sampels, columns: receptors
    zscore : bool
        if the input DataFrame zscore, if True, calculate ranks to avoid multiplication with 2 negative numbers
    
    Returns
    -------
    dict
       Correlation values of the receptor activity and ligand x receptor TPM/ranks, keys: receptors, values: r (correlation)
    """
    genes_of_receptors = list(set(input.columns) & set(receptor_activities.columns))
    gene_receptor_corr_tg = {}
    # no_no_target = []
    for gene in genes_of_receptors:
        ligtarget_df = lig_rec[lig_rec['target_genesymbol'] == gene]
        if len(ligtarget_df) == 0: 
        #     no_no_target.append(gene)
            continue

        for idx in ligtarget_df.index:
            targetgene = ligtarget_df.source_genesymbol.loc[idx]

            if targetgene not in input.columns:
                continue
            if all(input[targetgene] == 0):
                continue
            if zscore:
                method = 'spearman'
                # ranks
                gene_values = input.loc[:, gene].rank(method = 'min') * input.loc[:, targetgene].rank(method = 'min')
            else: # tpm
                method = 'pearson'
                gene_values = input.loc[:, gene].astype('float') * input.loc[:, targetgene].astype('float')
            correlation = receptor_activities.loc[:, gene].astype('float')\
                .corr(gene_values, method = method)
            gene_receptor_corr_tg[gene + '_' + targetgene] = correlation # gene: receptor, targetgene: ligand
    # if return_no_target_list: return no_no_target, gene_receptor_corr_tg
    return gene_receptor_corr_tg

def ligand_activity_and_receptorxligand_correlation(lig_rec: DataFrame, input_data:DataFrame, ligand_activities:DataFrame, zscore = True): # return_no_target_list = False
    """Correlation of receptor activity and mean of ligand x its receptors TPM/ranks
    
    Parameters
    ----------
    lig_rec : DataFrame
        prior knowledge of ligand - receptor interactions
    input : DataFrame
        TPM or z-score transformed values of samples' gene expression
    ligand_activities:  DataFrame
        estimated ligand activities, rows: sampels, columns: ligands
    zscore : bool
        if the input DataFrame zscore, if True, calculate ranks to avoid multiplication with 2 negative numbers
    
    Returns
    -------
    dict
       Correlation values of the receptor activity and ligand x receptor TPM/ranks, keys: receptors, values: r (correlation)
    """
    genes_of_ligands = list(set(input_data.columns) & set(ligand_activities.columns))
    gene_ligand_corr_tg = {}
    for gene in genes_of_ligands:
        ligtarget_df = lig_rec[lig_rec['source_genesymbol'] == gene] # get receptors of ligand
        if len(ligtarget_df) == 0: 
            continue
        corr_list = []
        for targetreceptorgene in ligtarget_df.target_genesymbol: # iterate over receptors
            if targetreceptorgene not in input_data.columns:
                continue
            if all(input_data[targetreceptorgene] == 0): # if all TPM == 0
                continue
            if zscore:
                method = 'spearman'
                gene_values = input_data.loc[:, gene].rank(method = 'min') * input_data.loc[:, targetreceptorgene].rank(method = 'min')
            else: # tpm
                method = 'pearson'
                gene_values = input_data.loc[:, gene].astype('float') * input.loc[:, targetreceptorgene].astype('float') # ligand x receptor
            correlation = ligand_activities.loc[:, gene].astype('float').corr(gene_values, method = method)
            gene_ligand_corr_tg[gene + targetreceptorgene] = correlation # gene: ligand , targetreceptorgene: receptor
    return gene_ligand_corr_tg

def perform_permutations(input_data: DataFrame):
    permuted_dataframe = input_data.copy() 
    random.seed(29723)
    permuted_genenames = random.sample(list(input_data.index), len(input_data.index))
    permuted_dataframe.index = permuted_genenames
    permuted_dataframe = permuted_dataframe.loc[input_data.index]
    return permuted_dataframe

def background_correlation(input_data: DataFrame, activities: DataFrame):
    """Calculation of background correlation via gene name permutations
    
    Parameters
    ----------
    input_data : DataFrame
        TPM or z-score transformed values of samples' gene expression
    receptor_activities : DataFrame
        estimated receptor activities, rows: sampels, columns: receptors
    
    Returns
    -------
    dict
       Correlation values of the receptor activity and receptor gene (TPM or zscore); keys: receptors, values: r (correlation)
    """
    z_perm = perform_permutations(input_data)
    correlation = activity_and_gene_correlation(z_perm, activities)
    return correlation

def background_correlation_repeated(input_data, activities, repetition = 10):
    """Repeat background correlation and calculate the mean of all repeated correlations per receptor
    
    Parameters
    ----------
    input_data : DataFrame
        TPM or z-score transformed values of samples' gene expression
    activities : DataFrame
        estimated activities, rows: sampels, columns: receptors/ligands
    x : int
        number of repetition, default = 10
    
    Returns
    -------
    dict
       Mean correlation values of activtities and genes where gene name permutation repeated 
    """
    correlation_perm = {}
    for i in range(0, repetition):
        print(i, end = ' ')
        corr = (background_correlation(input_data, activities))
        correlation_perm[i] = corr
    print('...', end = ' ')
    correlation_perm = DataFrame.from_dict(correlation_perm)
    correlation_perm = correlation_perm.mean(1)
    correlation_perm = dict(correlation_perm)
    # delete key-values where the correlation is nan due to the original data sum()=0
    correlation_perm = {k: v for k, v in correlation_perm.items() if v is not None and not isna(v)}

    print('done')

    return correlation_perm


def background_correlation_receptor(lig_rec: DataFrame, input_data: DataFrame, activities: DataFrame):
    """Calculation of background correlation via gene name permutations of ligand and its receptors
    
    Parameters
    ----------
    lig_rec: 
        prior knowledge of ligand - receptor interactions
    input_data : DataFrame
        TPM or z-score transformed values of samples' gene expression
    receptor_activities : DataFrame
        estimated receptor activities, rows: sampels, columns: receptors
    
    Returns
    -------
    dict
       Correlation values of the ligand activity and receptor gene (TPM or zscore); keys: ligands, values: r (correlation)
    """
    z_perm = perform_permutations(input_data)
    correlation = activity_and_gene_correlation(z_perm, activities)
    return correlation

def background_correlation_receptor_repeated(lig_rec: DataFrame, input_data: DataFrame, activities: DataFrame, repetition = 10):
    """Repeat background correlation and calculate the mean of all repeated correlations per receptor
    
    Parameters
    ----------
    lig_rec : DataFrame
        prior knowledge of ligand - receptor interactions
    input_data : DataFrame
        TPM or z-score transformed values of samples' gene expression
    activities : DataFrame
        estimated activities, rows: sampels, columns: receptors/ligands
    x : int
        number of repetition, default = 10
    
    Returns
    -------
    dict
       Mean correlation values of ligand activtities and receptor where gene name permutation repeated 
    """
    correlation_perm = {}
    for i in range(0, repetition):
        print(i, end = ' ')
        corr = (background_correlation_receptor(lig_rec, input_data, activities))
        correlation_perm[i] = corr
    print('...', end = ' ')
    correlation_perm = DataFrame.from_dict(correlation_perm)
    correlation_perm = correlation_perm.mean(1)
    correlation_perm = dict(correlation_perm)
    # delete key-values where the correlation is nan due to the original data sum()=0
    correlation_perm = {k: v for k, v in correlation_perm.items() if v is not None and not isna(v)}
    print('done')

    return correlation_perm

def receptor_activity_and_ligand_correlation(lig_rec: DataFrame, input:DataFrame, receptor_activities:DataFrame, zscore = False): # return_no_target_list = False
    """Correlation of receptor activity and mean of ligand x its receptors TPM/ranks
    
    Parameters
    ----------
    lig_rec : DataFrame
        prior knowledge of ligand - receptor interactions
    input : DataFrame
        TPM or z-score transformed values of samples' gene expression
    receptor_activities:  DataFrame
        estimated receptor activities, rows: sampels, columns: receptors
    zscore : bool
        if the input DataFrame zscore, if True, calculate ranks to avoid multiplication with 2 negative numbers
    
    Returns
    -------
    dict
       Correlation values of the receptor activity and ligand x receptor TPM/ranks, keys: receptors, values: r (correlation)
    """
    genes_of_receptors = list(set(input.columns) & set(receptor_activities.columns))
    gene_receptor_corr_tg = {}
    # no_no_target = []
    for gene in genes_of_receptors:
        ligtarget_df = lig_rec[lig_rec['target_genesymbol'] == gene]
        if len(ligtarget_df) == 0: 
        #     no_no_target.append(gene)
            continue

        for idx in ligtarget_df.index:
            targetgene = ligtarget_df.source_genesymbol.loc[idx]

            if targetgene not in input.columns:
                continue
            if all(input[targetgene] == 0):
                continue
            if zscore:
                method = 'spearman'
                # ranks
                gene_values = input.loc[:, targetgene].rank(method = 'min')
            else: # tpm
                method = 'pearson'
                gene_values = input.loc[:, targetgene].astype('float')
            correlation = receptor_activities.loc[:, gene].astype('float')\
                .corr(gene_values, method = method)
            gene_receptor_corr_tg[gene + '_' + targetgene] = correlation # gene: receptor, targetgene: ligand
    # if return_no_target_list: return no_no_target, gene_receptor_corr_tg
    return gene_receptor_corr_tg


# NEW
def background_correlation_receptor_ligand_repeated(lig_rec: DataFrame, input_data: DataFrame, receptor_activities: DataFrame, repetition = 10):
    """Repeat background correlation and calculate the mean of all repeated correlations per receptor
    
    Parameters
    ----------
    lig_rec : DataFrame
        prior knowledge of ligand - receptor interactions
    input_data : DataFrame
        TPM or z-score transformed values of samples' gene expression
    receptor_activities : DataFrame
        estimated activities, rows: sampels, columns: receptors/ligands
    repetition : int
        number of repetition, default = 10
    
    Returns
    -------
    dict
       Mean correlation values of ligand activtities and receptor where gene name permutation repeated 
    """
    correlation_perm = {}
    for i in range(0, repetition):
        print(i, end = ' ')
        corr = (background_correlation_receptor_ligand(lig_rec, input_data, receptor_activities))
        correlation_perm[i] = corr
    print('...', end = ' ')
    correlation_perm = DataFrame.from_dict(correlation_perm)
    correlation_perm = correlation_perm.mean(1)
    correlation_perm = dict(correlation_perm)
    # delete key-values where the correlation is nan due to the original data sum()=0
    correlation_perm = {k: v for k, v in correlation_perm.items() if v is not None and not isna(v)}
    print('done')

    return correlation_perm


def background_correlation_ligand_ligandxreceptor(lig_rec: DataFrame, input_data: DataFrame, ligand_activities: DataFrame):
    """Calcualtion background correlation via gene name permutations of ligand activity and ligandxreceptor activity
    
    Parameters
    ----------
    lig_rec: DataFrame
        prior knowledge of ligand - receptor interactions
    input_data : DataFrame
        TPM or z-score transformed values of samples' gene expression
    ligand_activities : DataFrame
        estimated receptor activities, rows: sampels, columns: receptors
    
    Returns
    -------
    dict
       Correlation values of the receptor activity and receptor gene (TPM or zscore); keys: receptors, values: r (correlation)
    """
    z_perm = perform_permutations(input_data)
    correlations = ligand_activity_and_receptorxligand_correlation(lig_rec, z_perm, ligand_activities)
    return correlations



def background_correlation_ligand_ligandxreceptor_repeated(lig_rec: DataFrame, input_data: DataFrame, activities: DataFrame, repetition = 10):
    """Repeat background correlation and calculate the mean of all repeated correlations per ligand
    
    Parameters
    ----------
    lig_rec : DataFrame
        prior knowledge of ligand - receptor interactions
    input_data : DataFrame
        TPM or z-score transformed values of samples' gene expression
    activities : DataFrame
        estimated activities, rows: sampels, columns: receptors/ligands
    x : int
        number of repetition, default = 10
    
    Returns
    -------
    dict
       Mean correlation values of ligand activtities and ligand x receptor where gene name permutation repeated 
    """
    correlation_perm = {}
    for i in range(0, repetition):
        print(i, end = ' ')
        corr = (background_correlation_ligand_ligandxreceptor(lig_rec, input_data, activities))
        correlation_perm[i] = corr
    print('...', end = ' ')
    correlation_perm = DataFrame.from_dict(correlation_perm)
    correlation_perm = correlation_perm.mean(1)
    correlation_perm = dict(correlation_perm)
    # delete key-values where the correlation is nan due to the original data sum()=0
    correlation_perm = {k: v for k, v in correlation_perm.items() if v is not None and not isna(v)}
    print('done')

    return correlation_perm


def background_correlation_receptor_ligandxreceptor(lig_rec: DataFrame, input_data: DataFrame, ligand_activities: DataFrame):
    """Calcualtion background correlation via gene name permutations of ligand activity and ligandxreceptor activity
    
    Parameters
    ----------
    lig_rec: DataFrame
        prior knowledge of ligand - receptor interactions
    input_data : DataFrame
        TPM or z-score transformed values of samples' gene expression
    ligand_activities : DataFrame
        estimated receptor activities, rows: sampels, columns: receptors
    
    Returns
    -------
    dict
       Correlation values of the receptor activity and receptor gene (TPM or zscore); keys: receptors, values: r (correlation)
    """
    z_perm = perform_permutations(input_data)
    correlations = receptor_activity_and_receptorxligand_correlation(lig_rec, z_perm, ligand_activities)
    return correlations

def background_correlation_receptor_ligand(lig_rec: DataFrame, input_data: DataFrame, receptor_activities: DataFrame):
    """Calcualtion background correlation via gene name permutations of ligand activity and ligandxreceptor activity
    
    Parameters
    ----------
    lig_rec: DataFrame
        prior knowledge of ligand - receptor interactions
    input_data : DataFrame
        TPM or z-score transformed values of samples' gene expression
    receptor_activities : DataFrame
        estimated receptor activities, rows: sampels, columns: receptors
    
    Returns
    -------
    dict
       Correlation values of the receptor activity and receptor gene (TPM or zscore); keys: receptors, values: r (correlation)
    """
    z_perm = perform_permutations(input_data)
    correlations = receptor_activity_and_ligand_correlation(lig_rec, z_perm, receptor_activities)
    return correlations



def background_correlation_receptor_ligandxreceptor_repeated(lig_rec: DataFrame, input_data: DataFrame, activities: DataFrame, repetition = 10):
    """Repeat background correlation and calculate the mean of all repeated correlations per ligand
    
    Parameters
    ----------
    lig_rec : DataFrame
        prior knowledge of ligand - receptor interactions
    input_data : DataFrame
        TPM or z-score transformed values of samples' gene expression
    activities : DataFrame
        estimated activities, rows: sampels, columns: receptors/ligands
    x : int
        number of repetition, default = 10
    
    Returns
    -------
    dict
       Mean correlation values of ligand activtities and ligand x receptor where gene name permutation repeated 
    """
    correlation_perm = {}
    for i in range(0, repetition):
        print(i, end = ' ')
        corr = (background_correlation_receptor_ligandxreceptor(lig_rec, input_data, activities))
        correlation_perm[i] = corr
    print('...', end = ' ')
    correlation_perm = DataFrame.from_dict(correlation_perm)
    correlation_perm = correlation_perm.mean(1)
    correlation_perm = dict(correlation_perm)
    # delete key-values where the correlation is nan due to the original data sum()=0
    correlation_perm = {k: v for k, v in correlation_perm.items() if v is not None and not isna(v)}
    print('done')

    return correlation_perm