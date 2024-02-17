
from pandas import DataFrame, concat
import random
from tqdm import tqdm
from numpy import array, mean, std


def permute(genename_list:list, seed:int):
    random.seed(seed)
    random.shuffle(genename_list)
    return genename_list

def permute_genenames(input_data:DataFrame, number_of_permutation:int):
    permuted_genenames_dict = {}
    permuted_genenames_dict[0] = list(input_data.columns)
    for s in range(1, number_of_permutation):
        # if s%100 == 0:
        #     print(s, end = ', ')
        genenames = list(input_data.columns)
        permuted_genenames_dict[s] = permute(genenames, s)
    return permuted_genenames_dict


def change_gene_name_order_based_on_permutations(chunk_df:DataFrame, permuted_genename_list:list):
    chunk_df_permuted = chunk_df.copy()
    chunk_df_permuted.columns = permuted_genename_list
    return chunk_df_permuted

def estimation_dataframe(lincs_model:DataFrame, chunk_df: DataFrame):
    genes = lincs_model.index.intersection(chunk_df.columns)
    # if the column names are duplicated - error
    return chunk_df.loc[:, genes] @ lincs_model.loc[genes]

def estimation_with_zscore_calculation_in_chunks(input_data:DataFrame, lincs_model:DataFrame, number_of_permutation: int, chunk_size = 100):
    # permute genenames and store as list of indexes
    permuted_genename_lists = permute_genenames(input_data,number_of_permutation)

    num_samples = len(input_data)
    print('Number of samples:', num_samples)
    chunks = [input_data.index[i:i+chunk_size] for i in range(0, num_samples, chunk_size)]
    print('Number of chunks:', len(chunks))
    receptors = list(lincs_model.columns)
    permutations = range(0, len(permuted_genename_lists.keys())) # number of permutation
    print('Number of permutations:', len(permutations))


    zscore_receptor_activities = {}
    i = 0
    for chunk in tqdm(chunks, total=len(chunks)):
        chunk_df = input_data.loc[chunk].copy() # subset dataframe

        chunk_estimated_values_dict = {}

        for permutation in permutations:
            chunk_df_permuted = change_gene_name_order_based_on_permutations(chunk_df, permuted_genename_lists[permutation])
            # estimated values of a sample and a receptor
            estimated_values_df = estimation_dataframe(lincs_model, chunk_df_permuted)
            chunk_estimated_values_dict[permutation] = estimated_values_df

            
        array_3d = array([df.values for df in chunk_estimated_values_dict.values()])
        mean_1d = mean(array_3d, axis=0)
        std_1d = std(array_3d, axis=0)
        chunk_receptor_activities_zscore = (chunk_estimated_values_dict[0] - mean_1d) / std_1d

        zscore_receptor_activities[i] = chunk_receptor_activities_zscore
        i = i+1


    zscore_receptor_activities_sample_receptor= concat(zscore_receptor_activities, axis=0)
    zscore_receptor_activities_sample_receptor.index = zscore_receptor_activities_sample_receptor.index.get_level_values(1)

    return zscore_receptor_activities_sample_receptor