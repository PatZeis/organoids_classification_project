
import os
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import pickle
from scipy.spatial import distance
from scipy import stats
import sys
from matplotlib import pyplot as plt

rng = np.random.default_rng(0)



def print_var(var_p, name_p, all_p=False):
    
    print(name_p)
    print()
    
    if type(var_p) == list:
        
        if all_p == False:
            print(var_p[0:10])
        else:
            print(var_p)
            
        print()
        print('Length:')
        print(len(var_p))
        
    else:
        print(var_p)

        try:
            print(var_p.shape)
        except:
            print()
    
    print()
    print()



def prepare_model_adata(adata):

    """
    Filter out some cell types (CM, immune and vascular).
    Create a new cell type group annotation (cell_type_group_mod)
    where amacrine, bipolar and horizontal cell type group are modified,
    for example: horizontal --> HC_01 and HC_02.

    Parameters:
        adata: Adata.

    Return:
        adata_filt: Filtered adata

    """

    lab_ct = 'cell_type'
    lab_ct_group = 'cell_type_group'

    adata_filt = adata.copy()

    adata_filt.obs.loc[:, lab_ct_group] = adata_filt.obs.loc[:, lab_ct_group].astype('str')
    mask_rpe = adata_filt.obs.loc[:,"cell_type"] == "RPE"
    mask_cm = adata_filt.obs.loc[:,"cell_type"] == "CM"
    celltype_cat = list(adata_filt.obs.loc[:,"cell_type"])
    cell_types_strip = [e.split("_")[0] for e in celltype_cat]
    mask_fb =  list(np.array(cell_types_strip) == "FB")
    mask_end = list(np.array(cell_types_strip) == "END")
    mask_per =  list(np.array(cell_types_strip) == "PER")
    mask_t = list(np.array(cell_types_strip) == "TCell")
    mask_nk =  list(np.array(cell_types_strip) == "NK")
    mask_mo = list(np.array(cell_types_strip) == "MO")
    mask_uG = list(np.array(cell_types_strip) == "uG")
    mask_mast =  list(np.array(cell_types_strip) == "MAST")  

    adata_filt.obs.loc[mask_rpe, lab_ct_group] = "RPE"
    adata_filt.obs.loc[mask_cm, lab_ct_group] = "CM"
    adata_filt.obs.loc[mask_fb, lab_ct_group] = "FB"
    adata_filt.obs.loc[mask_end, lab_ct_group] = "END"
    adata_filt.obs.loc[mask_per, lab_ct_group] = "PER"
    adata_filt.obs.loc[mask_t, lab_ct_group] = "TCell"
    adata_filt.obs.loc[mask_nk, lab_ct_group] = "NK"
    adata_filt.obs.loc[mask_mo, lab_ct_group] = "MO"
    adata_filt.obs.loc[mask_uG, lab_ct_group] = "uG"
    adata_filt.obs.loc[mask_mast, lab_ct_group] = "MAST"

    adata_filt.obs.loc[:, lab_ct_group] = adata_filt.obs.loc[:, lab_ct_group].astype('category')


    adata_filt.obs.loc[:, 'cell_type_group_mod'] = adata_filt.obs.loc[:, 'cell_type_group'].copy(deep=True).astype('str')
    
    for ct in ['amacrine', 'bipolar', 'horizontal']:
        mask = adata_filt.obs.loc[:, 'cell_type_group'] == ct

        s = adata_filt.obs.groupby(by=['cell_type_group', 'cell_type', 'subject']).size()

        if ct == 'amacrine':
            s = adata_filt.obs.loc[mask, 'cell_type'].str.split('_').str[0:2]

            s = s.str.join('_')

        if ct == 'bipolar':
            s = adata_filt.obs.loc[mask, 'cell_type'].str.split('_').str[0]

        if ct == 'horizontal':
            s = adata_filt.obs.loc[mask, 'cell_type']



        adata_filt.obs.loc[mask, 'cell_type_group_mod'] = s

    adata_filt.obs.loc[:, 'cell_type_group_mod'] = adata_filt.obs.loc[:, 'cell_type_group_mod'].astype('category')



    return adata_filt



def get_markers(adata, lab_ct_group, n_genes, path_results, l_test_ct=None):

    """
    Calculate differential expression genes.

    Parameters:
        adata: Adata
        lab_ct_group: Name of the cell type annotation to use as group for the differential expression analysis.
        n_genes: Number of genes to select.
        path_results: Path where to save the results.
        l_test_ct: If l_test_ct=None compare a group vs rest.
                   If l_test_ct=list, e.g. [A, B], compare A vs B.

        Return:
            l_markers: List of markers.
    """
    
    l_markers = []
    
    if l_test_ct is None:
    
        sc.tl.rank_genes_groups(adata, groupby=lab_ct_group, method='wilcoxon', use_raw=False, n_genes=n_genes)

        l_ct = sorted(adata.obs.loc[:, lab_ct_group].unique().tolist())

        for ct in l_ct:

            ct_diff_exp = sc.get.rank_genes_groups_df(adata, group=ct, key='rank_genes_groups')

            if ct_diff_exp.shape[0] < n_genes:
                sys.exit('Error: n_genes < {} in {}'.format(n_genes, ct))

            ct_markers = ct_diff_exp.loc[:, 'names'].to_list()
            l_markers = l_markers + ct_markers
        
    else:

        for ct_1 in l_test_ct:

            for ct_2 in l_test_ct:

                if ct_1 != ct_2:
                    
                    mask = (adata.obs.loc[:, lab_ct_group] == ct_1) | (adata.obs.loc[:, lab_ct_group] == ct_2)
                    adata_filt = adata[mask, :].copy()
                    
                    sc.tl.rank_genes_groups(adata_filt, groupby=lab_ct_group, method='wilcoxon', use_raw=False, n_genes=n_genes, groups=[ct_1], rest=ct_2)
                    
                    test_name = '{}_vs_{}'.format(ct_1, ct_2)

                    ct_diff_exp = sc.get.rank_genes_groups_df(adata_filt, group=ct_1, key='rank_genes_groups')
                    
                    if ct_diff_exp.shape[0] < n_genes:
                        sys.exit('Error: n_genes < {} in {}'.format(n_genes, ct))
                    
                    ct_markers = ct_diff_exp.loc[:, 'names'].to_list()
                    l_markers = l_markers + ct_markers

    l_markers = sorted(set(l_markers))
    
    return l_markers



def get_markers_two_groups(adata, n_genes, path_results, l_markers, dataset_name):

    """
    Calculate differential expression genes comparing two groups instead of one group vs rest.

    Parameters:
        adata: Adata.
        n_genes: Number of genes.
        path_results: Path where to save the results.
        l_markers: List of markers to exapand with new markers.
        dataset_name: Name to add to the marker file name.
        

    Return:
            l_markers: List of markers.

    """

    l_test_ct = [['AC_B', 'AC_Y'], ['HC_01', 'HC_02'], ['CdBC', 'ChBC', 'RBC'], ['RPE', 'CM'], ['RPE', 'FB']]

    for test_ct in l_test_ct:

        l_markers = l_markers + get_markers(adata, 'cell_type_group_mod', n_genes, path_results, test_ct)
        l_markers = sorted(set(l_markers))
    
    return l_markers


def get_expr_all(adata_p, genes_p):
    
    """
    Get expression matrix of all cell types. Each cell is normalized using L2 norm.
    
    :param adata_p: Adata.
    :param genes_p: Genes to keep.
    
    :retrun df_in_f: Expression matrix.
            m_norm_in_f: L2 normalized expression matrix.
            
    """
    

    genes_not_present = sorted(list(set(genes_p).difference(set(adata_p.var_names))))
    
    if len(genes_not_present) > 0:

        print_var(genes_not_present, 'Error: Marker genes not present in adata.')
        sys.exit()


    adata_in_f = adata_p[:, genes_p].copy()
    
    df_in_f = pd.DataFrame(data=adata_in_f.X.todense(), columns=genes_p)

    m_norm_in_f = preprocessing.normalize(df_in_f.to_numpy(), norm='l2', axis=1)
    
    m_norm_in_f = pd.DataFrame(data=m_norm_in_f, columns=df_in_f.columns)
    
    len_in_f = np.square(m_norm_in_f.loc[0,:].to_numpy()).sum()
   

    return df_in_f, m_norm_in_f

def get_sd(adata, tab_test, tab_ref, k):
    
    """
    Get standardized distances.
    
    Parameters:
        adata: Adata.
        tab_test: Data from the dataset to test.
        tab_ref: Data from the model.
        k: Number of random cells to use in the distance calculation.
    """
    
    l_ct = sorted(tab_ref.index.unique().to_list())
    
    l_sd = []
    l_name_sd = []

    for ct in l_ct:
        
        print('{} standardized distance calculation.'.format(ct))
        print()

        tab_ref_ct = tab_ref.loc[ct, :].copy(deep=True)

        mean_ct_ref, mad_ct_ref = get_min_dist(tab_ref_ct, tab_ref_ct, k, True)

        min_dist_ct_test = get_min_dist(tab_test, tab_ref_ct, k, False)

        sd_ct = (min_dist_ct_test - mean_ct_ref) / mad_ct_ref
        
        l_sd.append(sd_ct)
        l_name_sd.append('{}_sd'.format(ct))
    
    l_sd = np.array(l_sd)
    l_sd = np.transpose(l_sd)
    
    adata.obsm['X_sd'] = l_sd
    adata.uns['sd_names'] = l_name_sd
    print('added o adata.obsm and adata.uns')

def get_min_dist(test_m, ref_m, k, is_reference):

    """
    Get minimum distances among each cell in test_m and a set of cells in ref_m.
    
    Parameters:
        test_m: Test sample.
        ref_m: Control sample.
        k: Number of cells to select.
        is_reference: True is used to build the model: the cells to test and the cells to test come from the same population,
                      it excludes the same cell in distance calulation since the distance would be 0.
  
    Return:
        mean_l_min_dist, mad_l_min_dist if is_reference is True. l_min_dist if is_reference is False.
    """
    
    ref_m = ref_m.to_numpy()
 
    test_m = test_m.to_numpy()

    l_dist = []
    l_min_dist = []
    check = True
    
    for i in range(0, test_m.shape[0]):
        
        if k >= ref_m.shape[0] and is_reference:
            ref_m_k_cells = np.delete(ref_m, i, axis=0)
            if check:
                check = False

        else:
            if k >= ref_m.shape[0]:
                ref_m_k_cells = ref_m
            else:
                idx = rng.choice(ref_m.shape[0], replace=False, shuffle=True, size=k)
                if is_reference:
                    while i in idx:
                        idx = rng.choice(ref_m.shape[0], replace=False, shuffle=True, size=k)

                ref_m_k_cells = ref_m[idx, :]
            if check:
                check = False
        
        
        l_dist = []
        for j in range(0, ref_m_k_cells.shape[0]):
            
            l_dist.append(distance.euclidean(ref_m_k_cells[j, :], test_m[i, :]))
            
        l_dist = np.array(l_dist)
        
        min_dist_in_f = np.min(l_dist)
        
        l_min_dist.append(min_dist_in_f)
    
    l_min_dist = np.array(l_min_dist)
        
    if is_reference:
        
        mean_l_min_dist = np.mean(l_min_dist)
        mad_l_min_dist = stats.median_abs_deviation(l_min_dist, scale=1/1.4826)
        
        
        return mean_l_min_dist, mad_l_min_dist
    
    else:

        return l_min_dist
        
def get_sd_centroids(adata, tab_test,tab_ref):
       
    l_ct = sorted(tab_ref.index.unique().to_list())
    l_sd_ct = []
    c_name = []
    centroids = []
    tab_test = tab_test.to_numpy()
    
    for ct in l_ct:
           
        tab_ref_ct = tab_ref.loc[ct, :].copy(deep=True)
        m_ref = tab_ref_ct.mean(axis=0)
        l_dist = []
        for  j in range(0, tab_test.shape[0]):
            l_dist.append(distance.euclidean(tab_test[j, :], m_ref))
            
        l_dist = np.array(l_dist)
        
        l_dist_ref = []
        tab_ref_ct = tab_ref_ct.to_numpy()
        for  j in range(0, tab_ref_ct.shape[0]):
            l_dist_ref.append(distance.euclidean(tab_ref_ct[j, :], m_ref))
            
        l_dist_ref = np.array(l_dist_ref)
        u_l_dist_ref = np.mean(l_dist_ref)
        mad_l_dist_ref = stats.median_abs_deviation(l_dist_ref, scale=1/1.4826)
        sd_ct = (l_dist - u_l_dist_ref)/ mad_l_dist_ref
        l_sd_ct.append(sd_ct)
        centroids.append(m_ref)
        c_name.append('{}'.format(ct))
    l_centroids = np.array(centroids)
    l_centroids = np.transpose(l_centroids)
    l_dist_ct = np.array(l_sd_ct)
    l_dist_ct = np.transpose(l_sd_ct)
    adata.obsm['sd_centroids'] = l_dist_ct

def sd_calculation(adata, dataset_name, path_results, k, path_model=None, model_name=None, atlas_test=None):
    
    """
    Save standardized distances in adata.
    
    Parameters:
        adata: adata.
        dataset_name: Name to add to path files.
        path_results: Path results.
        k: Number of random cells to use in the distance calculation.
        path_model: Path where to find the model. If it is None create the model.
        model_name: Name of the model: peripheral_retina 
        
    Return:
        adata
    """
    
    adata = adata.copy()

    lab_ct_group = 'cell_type_group'
    
    if path_model is None:

        adata_model = adata

        mask_markers = adata_model.var.loc[:, 'markers']
        l_markers = adata_model.var.loc[mask_markers, :].index.to_list()
        l_markers_new = l_markers
        genes_not_present = sorted(list(set(l_markers).difference(set(l_markers_new))))

        expr, l2_norm_m = get_expr_all(adata, l_markers)
        
        path_adata_out = os.path.join(path_results, 'adata_{}.h5ad'.format(dataset_name))
        
        path_scaler_model = os.path.join(path_results, 'scaler_{}.pickle'.format(dataset_name))
        
        scaler = StandardScaler()
        scaler.fit(l2_norm_m)
        
        with open(path_scaler_model, 'wb') as f:

            pickle.dump(scaler, f, pickle.HIGHEST_PROTOCOL)

    else:

        adata_model = sc.read_h5ad(os.path.join(path_model, 'adata_{}.h5ad'.format(model_name)))
        
        mask_markers = adata_model.var.loc[:, 'markers']
        l_markers = adata_model.var.loc[mask_markers, :].index.to_list()
        
        adata.var.loc[:, 'markers'] = adata.var.index.isin(l_markers)
        mask_markers_new = adata.var.loc[:, 'markers'] 
        l_markers_new = adata.var.loc[mask_markers_new, :].index.to_list()
         
        

        expr, l2_norm_m = get_expr_all(adata, l_markers_new)

        path_adata_out = os.path.join(path_results, 'adata_{}.h5ad'.format(dataset_name))
        
        genes_not_present = sorted(list(set(l_markers).difference(set(l_markers_new))))
         
        if len(genes_not_present) == 0:
            print("marker genes match in model and test data")
            path_scaler_model = os.path.join(path_model, 'scaler_{}.pickle'.format(model_name))

            with open(path_scaler_model, 'rb') as f:

                scaler = pickle.load(f)
        else:
            print("marker genes do NOT match in model and test data")
            expr_model, l2_norm_m_model = get_expr_all(adata_model, l_markers_new)
            print("done get expression")
            scaler = StandardScaler()
            print("done standardscaler")
            scaler.fit(l2_norm_m_model)
            print("done fit")

    l2_norm_stand_scal_m = scaler.transform(l2_norm_m)

    df_l2_norm_stand_scal_m = pd.DataFrame(l2_norm_stand_scal_m, columns=l_markers_new)
    
    print("scaling done") 
    if path_model is None:
        
        path_pca_model = os.path.join(path_results, 'pca_{}.pickle'.format(dataset_name))

        pca = PCA(n_components=10, svd_solver='full', random_state=0)
        pca.fit(l2_norm_stand_scal_m)

        with open(path_pca_model, 'wb') as f:
            pickle.dump(pca, f, pickle.HIGHEST_PROTOCOL)
    else:
        if len(genes_not_present) == 0:

            path_pca = os.path.join(path_model, 'pca_{}.pickle'.format(model_name))

            with open(path_pca, 'rb') as f:

                pca = pickle.load(f)
        else:
            pca = PCA(n_components=10, svd_solver='full', random_state=0)
            l2_norm_stand_scal_m_model = scaler.transform(l2_norm_m_model)
            pca.fit(l2_norm_stand_scal_m_model)
            pca_scores_model = pca.transform(l2_norm_stand_scal_m_model)
            df_pca_scores_model = pd.DataFrame(data=pca_scores_model)
            df_pca_scores_model = df_pca_scores_model.set_index(adata_model.obs.loc[:, lab_ct_group])
    
    pca_scores = pca.transform(l2_norm_stand_scal_m)
    print("done pca scoring") 
    adata.obsm['X_pca_sd'] = pca_scores
    
    df_pca_scores = pd.DataFrame(data=pca_scores)
    if len(genes_not_present) == 0:
        df_pca_scores_model = pd.DataFrame(adata_model.obsm['X_pca_sd'])
    
        df_pca_scores_model = df_pca_scores_model.set_index(adata_model.obs.loc[:, lab_ct_group])
    
    get_sd(adata, df_pca_scores, df_pca_scores_model, k)
    print('get sd done')
    get_sd_centroids(adata, df_pca_scores,df_pca_scores_model)
    print('get_sd_centroids')
    if path_model is not None and atlas_test is None:
        adata.var.mt = adata.var.mt.astype('str') 
        adata.var.ribosomal = adata.var.ribosomal.astype('str') 
        adata.var.highly_variable = adata.var.highly_variable.astype('str')
    
        mt_var = np.array(adata.var.mt)
        mt_var[mt_var == 'nan'] = "False"
        adata.var.mt = mt_var
        adata.var.mt = adata.var.mt.astype('bool')

        rb_var = np.array(adata.var.ribosomal)
        rb_var[rb_var == 'nan'] = "False"
        adata.var.ribosomal = rb_var
        adata.var.ribosomal = adata.var.ribosomal.astype('bool')

        hv_var = np.array(adata.var.highly_variable)
        hv_var[hv_var == 'nan'] = "False"
        adata.var.highly_variable = hv_var
        adata.var.highly_variable = adata.var.highly_variable.astype('bool')

    adata.write(path_adata_out, compression='gzip')
    print('writing adata file'+path_adata_out+'.gz done') 
    return adata

