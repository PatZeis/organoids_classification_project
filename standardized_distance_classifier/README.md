# Single-cell classifier

## Authors
Riccardo Panero
Patrice Zeis
## Features of the classifier:

* Uses peripheral retina as reference cell atlas
* Uses a quantitative metric - standardized distances to compare new cells to adult reference cell types.
* Identifies mature cell types

## Description of the method:
The goal of this method is to classify retinal cells from single-cell RNA-seq culture experiments such as organoids.\
We utilise models the peripheral retina from Cowan and Renner et al.[^1] , as reference.\
In order to classifiy cells of a single-cell experiment we put its cells in the context of these models (PCA embedding) and we calculate the standardized distance (see `standardized distance calculation`) of each cell with respect to a subset of cells in the peripheral retina.\
In addition,in order to assess the similarity of query cells to reference cell type manifold, for each reference cell type the centroid is calculated and for each query cell the standardize distance to these centroids is calculated.\ 

## Standardized distance calculation:
Distances are calculated by measuring the minimum Euclidean distance of each cell with respect to k random adult cells or to the reference cell type centroids.\
After this, they are standardized considering the mean and the MAD of the adult cell type (k random cells or reference cells centroid distances).
These standardized distances are saved in 'adata.obsm.X_sd' or 'adata.obsm.sd_centroids'and their names are saved in `adata.uns.sd_names`.
If the standarized distances of 'adata.obsm.X.sd' is <3 then the cell is classified to its closest reference cell type. If the standardised distance of a cell is >= 3 then the cells is classified as immature..\

     
* ### Create new models:
     This option allows the generation of a adult reference model that can be used to perform the standardized distance calculation\
     Users can change the default number of markers to generate the peripheral model.\
     To create the model run the script `create_model.py` in `models/scripts`.
     \
     \
     Command lines used to generate the pre-generated model:
     ```
     python create_model.py -input_adata ../../adata_final_clean_periphery_tf.h5ad -k 1000 -name_output_adata=peripheral_retina  

     ```
     The `peripheral_retina_model` can be found in `models\results\peripheral_retina`.\
     The command line uses the default number of markers, to change this number check the help of the function: \
     ```python create_model.py --help```
     
     Content of the `models/results/peripheral_retina` folder:.
     
     * **adata_peripheral_retina.h5ad:** AnnData file
     * **pca_peripheral_retina.pickle:** PCA used for the embedding.
     * **scaler_peripheras_retina.pickle:** Standard scaler used to normalize the data before the embedding.
     
* ### Classify new cells using pre-generated models:
     To use the classifier, model generated before is used which is present in the `models\results` folder.\
     This models is named `peripheral_retina_model` 
     To run the classification use the `calculate_sd.py` script in `classifier_mature\scripts`
     \
     \
     Example:
     ```
     python calculate_sd.py -input_adata /path_to_dataset/dataset_to_classify.h5ad
                            -name_peripheral_model peripheral_retina_model
                            -k 1000
                            -name_output_adata classified_dataset
     ```

### References: ###
[^1]:https://doi.org/10.1016/j.cell.2020.08.013
