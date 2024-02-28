from importlib.resources import path
import sys
import os
import scanpy as sc
import argparse

path_model = '../../models'

sys.path.append(os.path.join(path_model, 'scripts'))
from functions import *


parser = argparse.ArgumentParser(description='Calculate standard ditances.')

parser.add_argument('-input_adata', type=str, required=True, help='Path of the input adata.')
parser.add_argument('-name_peripheral_model', help='Name of the peripheral model.')
parser.add_argument('-k', type=int, required=True, help='Number of cells used to calculate standardized distances.')
parser.add_argument('-name_output_adata', type=str, required=True, help='Name of the output adata saved in models/results/name_output_adata/name_output_adata.h5ad')
parser.add_argument('-atlas_test', type=str, required=False, help='is input a dataset which has been preprocessed as the samples in culture e.g. reference datasets')

args = parser.parse_args()

input_adata = args.input_adata
name_peripheral_model = args.name_peripheral_model
k = args.k
name_output_adata = args.name_output_adata
atlas_test = args.atlas_test
print(args)


def main():

    adata = sc.read_h5ad(input_adata)

    d_file_models = {name_peripheral_model: os.path.join(path_model, 'results', name_peripheral_model)}    
    adata_p = sc.read_h5ad(os.path.join(d_file_models[name_peripheral_model], 'adata_peripheral_retina.h5ad'))

    path_results = '../results'

    for model in sorted(d_file_models.keys()):
        
        print_var(model, '')

        path_results_1 = os.path.join(path_results, name_output_adata, 'tmp')

        if not os.path.exists(path_results_1):
            os.makedirs(path_results_1)

        input_adata_model = os.path.join(d_file_models[model])

        name_output_adata_1 = '{}_{}_{}'.format(name_output_adata, 'sd_from', model)
        
        adata_sd = sd_calculation(adata, name_output_adata_1, path_results_1, k, path_model=input_adata_model, model_name=model, atlas_test=atlas_test)

main()

print('End.')
