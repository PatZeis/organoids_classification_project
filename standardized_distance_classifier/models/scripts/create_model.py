from importlib.resources import path
import sys
import os
import scanpy as sc

from functions import *

import argparse

parser = argparse.ArgumentParser(description='Create model.')

parser.add_argument('-input_adata', type=str, required=True, help='Path of the input adata.')
parser.add_argument('-k', type=int, required=True, help='Number of cells used to calculate standardized distances.')
parser.add_argument('-n_markers_ct_vs_rest', type=int, required=False, default=10, help='Number of markers to keep after the cell type group differential expression analysis: each cell type group vs rest.\
                                                                                         Default = 10')

parser.add_argument('-n_markers_ct_vs_ct', type=int, required=False, default=3, help='Number of markers to keep after the cell type differential expression analysis.\
                                                                                      This analysis is restricted to amacrine, bipolar and horizontal cell type group and it compares\
                                                                                      each cell type to each other inside the same cell type group, e.g. horzizontal: HC_01 vs HC_02 and HC_02 vs HC_01.\
                                                                                      Default = 3')
parser.add_argument('-name_marker_file', type=str, required=False ,help='Name of the marker file', default=None)

parser.add_argument('-name_output_adata', type=str, required=True, help='Name of the output adata saved in models/results/name_output_adata/name_output_adata.h5ad')

args = parser.parse_args()
print_var(args, 'Parameter values:')


name_output_adata = args.name_output_adata
input_adata = args.input_adata
path_results = os.path.join('..', 'results', name_output_adata)
k = args.k
n_markers_ct_vs_rest = args.n_markers_ct_vs_rest
n_markers_ct_vs_ct = args.n_markers_ct_vs_ct
marker_file = args.name_marker_file
def main():

    adata = sc.read_h5ad(input_adata)
    adata = prepare_model_adata(adata)
    if marker_file is None:
        l_markers = get_markers(adata, 'cell_type_group', n_markers_ct_vs_rest, path_results)
        l_markers = get_markers_two_groups(adata, n_markers_ct_vs_ct, path_results, l_markers, name_output_adata)
 
    else:
        print("load marker file: "+marker_file) 
        markers = pd.read_csv(marker_file, index_col=0)
        l_markers = set(list(markers.iloc[:,0]))

    adata.var.loc[:, 'markers'] = adata.var.index.isin(l_markers)

    if not os.path.exists(path_results):
        os.makedirs(path_results)

    adata_sd = sd_calculation(adata, name_output_adata, path_results, k)

main()

print('End.')

