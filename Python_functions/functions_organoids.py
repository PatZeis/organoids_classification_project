#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: patrice.zeis
"""
 
import pandas as pd
import anndata as adata
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

def create_new_adata(adata_old,batch_nam, ref=False, empty=False):
    if "batch" not in adata_old.obs.columns:
        adata_old.obs['batch'] = batch_nam
    if ref:  
        adata_old.var = adata_old.var.rename(columns={'gene': 'gene_ids'})
    else:
        adata_old.obs['cell_type'] = adata_old.obs['batch']
        adata_old.obs['cell_type_group'] = adata_old.obs['batch']
        adata_old.obs['cell_type'] = adata_old.obs['cell_type'].astype('category')
        adata_old.obs['cell_type_group'] = adata_old.obs['cell_type_group'].astype('category')
        adata_old.obs['batch'] = adata_old.obs['batch'].astype('category')
     
    if empty:
        new_adata = adata.AnnData(X = adata_old.X)
        new_adata.var_names = adata_old.var_names
        new_adata.obs_names = adata_old.obs_names
        new_adata.obs["batch"] = adata_old.obs["batch"]
        new_adata.obs["cell_type"] = adata_old.obs["cell_type"]
        new_adata.obs["cell_type_group"] = adata_old.obs["cell_type_group"]
    
    else:
        new_adata = adata.AnnData(X = adata_old.X)
        new_adata.var_names = adata_old.var_names
        new_adata.obs_names = adata_old.obs_names
        new_adata.var["gene_ids"] = adata_old.var["gene_ids"]
        new_adata.obs["batch"] = adata_old.obs["batch"]
        new_adata.obs["cell_type"] = adata_old.obs["cell_type"]
        new_adata.obs["cell_type_group"] = adata_old.obs["cell_type_group"]
        if not ref:
            new_adata.obs["sd_class_mature"] = adata_old.obs["sd_class_mature"] 
            new_adata.obs["sd_class_immature"] = adata_old.obs["sd_class_immature"]
            new_adata.obs["sd_class"] = adata_old.obs["sd_class"]
            new_adata.obs["sd_region"] = adata_old.obs["sd_region"]
            new_adata.obs["pseudotime"] = adata_old.obs["pseudotime"]
    
    return new_adata

def bknn_integration(reference, explant,ref_nam, query_nam, path, ref_plot=False, markers=False,setmarkers=None, reftype="peripheral", leiden=1, masking=False, mask_type="tissue", mask_to_filt="chrRPE", drop_query=True, knn_key='batch'):
    sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(4, 4), facecolor='white')      
    if masking == True:
        mask_fib = reference.obs.loc[:,mask_type] == mask_to_filt
        reference = reference[mask_fib,:]
        explant_batch = explant.copy()
        explant_batch.obs['cell_type_group'] =  explant.obs["batch"]
        explant_batch.obs['batch'] = query_nam      
        reference_batch = reference.copy()
        reference_batch.obs['batch'] = ref_nam
        data_concat_for_sd = reference_batch.concatenate(explant_batch, batch_key="batch",batch_categories=[ref_nam, query_nam])
    if setmarkers is None:    
        if markers==True:
            if reftype=="peripheral":
                mask = np.array(explant.var.loc[:,"peripheral_markers"])
                explant_new = explant[:,mask]
                explant_new = create_new_adata(explant, query_nam, False, True)
            elif reftype=="foveal":
                mask = np.array(explant.var.loc[:,"foveal_markers"])
                explant_new = explant[:,mask]
                explant_new = create_new_adata(explant, query_nam, False, True)
            else:
                raise ValueError("set either peripheral or foveal")
    
        elif markers ==False:
            explant_new = create_new_adata(explant, query_nam, False, True)
        else:
            raise ValueError("set markers either True or False")
    else: 
        reference = reference[:,setmarkers]
        explant_new = create_new_adata(explant, query_nam, False, True)
 
    reference = create_new_adata(reference, ref_nam,True,False)  
        
    var_names = reference.var_names.intersection(explant_new.var_names)
    reference_new = reference[:, var_names]
    explant_new = explant_new[:, var_names]
    data_concat = reference_new.concatenate(explant_new, batch_key="batch",batch_categories=[ref_nam , query_nam,])
    data_concat.obs.cell_type_group = data_concat.obs.cell_type_group.astype('category')
    if setmarkers is None:
        if markers == True:
            sc.pp.highly_variable_genes(data_concat, n_top_genes=data_concat.shape[1]-2)
        else:    
            sc.pp.highly_variable_genes(data_concat, n_top_genes=3000)
    else:
        sc.pp.highly_variable_genes(data_concat, n_top_genes=data_concat.shape[1]-2)
        
    sc.tl.pca(data_concat, svd_solver='arpack')
    sc.external.pp.bbknn(data_concat, batch_key=knn_key)
    print("bknn done")
    sc.tl.umap(data_concat)
    sc.tl.leiden(data_concat, resolution=leiden)
    data_concat.obs.cell_type_group = data_concat.obs.cell_type_group.astype('str')   
    data_concat.obs.cell_type = data_concat.obs.cell_type.astype('str')
    mask_na = data_concat.obs.loc[:,"cell_type_group"] == "nan"
    data_concat.obs.loc[mask_na,"cell_type_group"] = query_nam
    data_concat.obs.loc[mask_na,"cell_type"] = query_nam
    data_concat.obs.cell_type_group = data_concat.obs.cell_type_group.astype('category')
    data_concat.obs.cell_type = data_concat.obs.cell_type.astype('category')  
    df = data_concat.obs[["cell_type_group", "leiden"]]
    contigency = pd.crosstab(index=df['cell_type_group'], columns=df['leiden'])
    if drop_query:
        contigency = contigency.drop(query_nam)
    celltypes_sum = contigency.apply(sum, axis=1)
    norm_contigency = contigency.div(celltypes_sum, axis="index")
    cluster_sum_norm = norm_contigency.apply(sum, axis=0)
    rel_norm_contigency = norm_contigency.div(cluster_sum_norm)
    cell_type_for_clus = pd.DataFrame({'celltype':rel_norm_contigency.idxmax(axis=0)})
    cell_type_for_clus.celltype = cell_type_for_clus.celltype.astype('str')    
    if drop_query:
        mask_na_leiden = cell_type_for_clus.loc[:,'celltype'] == "nan"
        cell_type_for_clus.loc[mask_na_leiden,'celltype'] = query_nam
     
    leiden_to_celltype = [cell_type_for_clus.loc[e][0] for e in list(df["leiden"]) ]
    data_concat.obs["leiden_to_celltype"] = leiden_to_celltype
    data_concat.obs.leiden_to_celltype = data_concat.obs.leiden_to_celltype.astype('category')
    sc.pl.umap(data_concat, color=['batch', 'cell_type_group', 'leiden_to_celltype','leiden'], return_fig=True)
    nam = path+"/"+query_nam+"_knn_integration_"+ref_nam+"_leiden"+str(leiden)+"_common_embedding_knnkey_"+knn_key
    if setmarkers is not None:
        nam = nam + "_" +reftype +"_markers_new_tf"
    elif markers == True:
        nam = nam + "_" +reftype +"_markers_new_tf"
    if masking:
        nam = nam + "_masking_" + mask_type + "_filtby_" + mask_to_filt
    if drop_query:
        nam = nam + "_querydropped"
    plt.savefig(nam + ".pdf")
    return data_concat, data_concat_for_sd


#### get standardized distance from classification 


def filter_cell_type_classified_centroids_sd(path_query,  to_filt, samples_list, samples_keys_list,filt_target=True,ref=None):
    samples = samples_list
    samples_key = samples_keys_list
    samples_key_new = [e.split("_")[0] for e in samples_key]  
    path_extend = "/classifier/results/"
    tmp = "/tmp/"
    adata_set = []
    for i, samp in enumerate(samples):
        adata_nam = "adata_classified_"+samp+"_sd_from_peripheral_retina.h5ad"
        sample_to_load = path_query+path_extend+"classified_"+samp+tmp+adata_nam
        to_load = sc.read_h5ad(sample_to_load)
        to_load.obs["batch"] = samples_key_new[i]
        adata_set.append(to_load)
        
    adata_concat = sc.AnnData.concatenate( *adata_set , batch_key = 'batch', batch_categories = samples_key_new ).copy() 
    mask_type="cell_type"
    mask_to_filt=to_filt
    reference = sc.read_h5ad(path_query+"/models/results/peripheral_retina/adata_peripheral_retina.h5ad")
    query_nam= "query_cells"
    ref_nam = "retina.peri"   
    mask_fib = reference.obs.loc[:,mask_type] == mask_to_filt
    reference = reference[mask_fib,:]
    explant_batch = adata_concat.copy()
    explant_batch.obs['cell_type_group'] = adata_concat.obs["batch"]
    explant_batch.obs['batch'] = query_nam      
    reference_batch = reference.copy()
    reference_batch.obs['batch'] = ref_nam
    data_concat_for_sd = reference_batch.concatenate(explant_batch, batch_key="batch",batch_categories=[ref_nam, query_nam])
    sd_query = data_concat_for_sd.obsm["X_sd"]
    ct_names_sd = [s.split("_")[0] for s in reference.uns["sd_names"]]
    sd_concat = pd.DataFrame(sd_query, columns=ct_names_sd)
    sd_concat.index = data_concat_for_sd.obs.index
    sd_df = sd_concat.copy(deep=True)     
    classified = sd_df.min(axis=1) < 3
    idx_min = sd_df.idxmin(axis=1).to_list()
    sd_df["mature_classified"] = list(classified.astype('str'))
    sd_df["mature_class"] = idx_min
    mask_mature = sd_df.loc[:,"mature_classified"] == "False"
    sd_df.loc[mask_mature,"mature_class"] = "immature"
    if filt_target == True:
        if ref is None:
            print("ref is none")
            mask_RPE = sd_df["mature_class"] == to_filt
            sd_df["batch"] = data_concat_for_sd.obs["cell_type_group"]
            sd_rpe = sd_df.loc[mask_RPE,:]
        else:
            sd_ref = pd.read_csv(ref, index_col=0) 
            mask_ref = sd_concat.index.isin(sd_ref.index)
            sd_df["batch"] = data_concat_for_sd.obs["cell_type_group"]
            sd_rpe = sd_df.loc[mask_ref,:]         
        
        sd_centroids = data_concat_for_sd.obsm["sd_centroids"]
        sd_centroids = pd.DataFrame(sd_centroids, columns=ct_names_sd)
        sd_centroids.index = data_concat_for_sd.obs.index
        filt_rpe = sd_centroids.index.isin(sd_rpe.index)
        sd_centroids_rpe_filt = sd_centroids.loc[filt_rpe,:]
        organoids = sc.read_h5ad(path_query+"/classifier/results/classified_celltype_organoids/tmp/adata_classified_celltypes_organoids_sd_from_peripheral_retina.h5ad")
        mask_org_fib = organoids.obs.loc[:,mask_type] == mask_to_filt
        organoids = organoids[mask_org_fib,:]
        organoids_sd_knn = pd.DataFrame(organoids.obsm["X_sd"], columns=ct_names_sd)
        organoids_sd_knn.index = organoids.obs.index
        classified = organoids_sd_knn.min(axis=1) < 3
        idx_min = organoids_sd_knn.idxmin(axis=1).to_list()
        organoids_sd_knn["mature_classified"] = list(classified.astype('str'))
        organoids_sd_knn["mature_class"] = idx_min
        mask_mature = organoids_sd_knn.loc[:,"mature_classified"] == "False"
        organoids_sd_knn.loc[mask_mature,"mature_class"] = "immature"
        organoids_sd_knn["batch"] = "celltype_organoids"
        sd_rpe_full = pd.concat([sd_rpe, organoids_sd_knn], axis=0)
        organoids_sd_cent = pd.DataFrame(organoids.obsm["sd_centroids"], columns=ct_names_sd)
        organoids_sd_cent.index = organoids.obs.index
        sd_centroids_full = pd.concat([sd_centroids_rpe_filt, organoids_sd_cent], axis=0)
        sd_rpe_full.to_csv(path_query+"/standardised_distance_data_knn_distance.csv")
        sd_centroids_full.to_csv(path_query+"/standardised_distance_data_centroid_distance.csv")
    else:
        sd_df["batch"] = data_concat_for_sd.obs["cell_type_group"]
        sd_df.to_csv(path_query+"/classification_standardised_distance_data_centroid_distance.csv")
        


