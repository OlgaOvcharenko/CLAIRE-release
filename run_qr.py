import argparse
import scanpy as sc
import anndata as ad
from os import listdir
from os.path import isfile, join
import pandas as pd
from os import walk
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, f1_score
import numpy as np
import faiss

only_files = ['results1/', 'results2/']
datasets = ['Pancreas_Mutaro', 'Pancreas_Segerstolpe', 'Pancreas_Wang', 'Pancreas_Xin']
batches = ['Mutaro_b2', 'Segerstolpe_b3', 'Wang_b4', 'Xin_b5']


cell_type_label = "celltype"
batch = "batch"

for i, dataset in enumerate(datasets):
    accs, f1_macros = [], []
    accs_pca, f1_macros_pca = [], []
    for folder in only_files:
        adata = sc.read_h5ad(f'/Users/olga_ovcharenko/Downloads/Olga_Data/{dataset.split("_")[0]}.h5ad')
        adata_ref = adata[adata.obs["batchlb"] != batches[i], :].copy()
        adata_query = adata[adata.obs["batchlb"] == batches[i],:].copy()

        embedding_ref = sc.read_h5ad("/Users/olga_ovcharenko/Downloads/outputs/" + f"{dataset}/{folder}/ad_100_ref.h5ad").X
        embedding_query = sc.read_h5ad("/Users/olga_ovcharenko/Downloads/outputs/" + f"{dataset}/{folder}/ad_100_query.h5ad").X
        
        exclude_bool = list(set(
            list(set(adata_query.obs[cell_type_label].tolist()) - set(adata_ref.obs[cell_type_label].tolist())) + 
            list(set(adata_ref.obs[cell_type_label].tolist()) - set(adata_query.obs[cell_type_label].tolist()))))
        exclude_ref = adata_ref.obs[cell_type_label].isin(exclude_bool)!= True
        adata_ref = adata_ref[exclude_ref]
        embedding_ref = embedding_ref[exclude_ref]

        exclude_query = adata_query.obs[cell_type_label].isin(exclude_bool)!= True
        adata_query = adata_query[exclude_query]
        embedding_query = embedding_query[exclude_query]

        knn = KNeighborsClassifier(n_neighbors=5)
        knn.fit(embedding_ref, adata_ref.obs[cell_type_label].tolist())
        cat_preds = knn.predict(embedding_query)
        f1_macro = f1_score(adata_query.obs[cell_type_label].tolist(), cat_preds, average='macro')
        acc = accuracy_score(adata_query.obs[cell_type_label].tolist(), cat_preds)
    
        # print(f1_score(adata_query.obs[cell_type_label].tolist(), cat_preds, labels=cell_types_list, average=None))
        
        accs.append(acc)
        f1_macros.append(f1_macro)

    print(f'{dataset} avg Accuracy {np.mean(accs).round(3)}, F1 macro {np.mean(f1_macro).round(3)}')
    print(f'{dataset} std Accuracy {np.std(accs).round(3)}, F1 macro {np.std(f1_macro).round(3)}')


only_files = ['results1/', 'results2/']
datasets = ['ImmuneAtlas']

cell_type_label = "cell_type"
batch = "batchlb"

for i, dataset in enumerate(datasets):
    accs, f1_macros, f1_weighteds = [], [], []
    
    adata = sc.read_h5ad(f'/Users/olga_ovcharenko/Downloads/Olga_Data/{dataset}.h5ad')
    adata_ref = adata[adata.obs["batchlb"] != "10x 5' v2", :].copy()
    adata_query = adata[adata.obs["batchlb"] == "10x 5' v2",:].copy()
    for folder in only_files:
        
        embedding_ref = sc.read_h5ad(f"/Users/olga_ovcharenko/Downloads/outputs/{dataset}_our_ref/{folder}ad_100_ref.h5ad").X
        embedding_query = sc.read_h5ad(f"/Users/olga_ovcharenko/Downloads/outputs/{dataset}_our_ref/{folder}/ad_100_query.h5ad").X

        exclude_bool = list(set(
            list(set(adata_query.obs[cell_type_label].tolist()) - set(adata_ref.obs[cell_type_label].tolist())) + 
            list(set(adata_ref.obs[cell_type_label].tolist()) - set(adata_query.obs[cell_type_label].tolist()))))

        exclude_ref = adata_ref.obs[cell_type_label].isin(exclude_bool)!= True
        adata_ref = adata_ref[exclude_ref]
        embedding_ref = embedding_ref[exclude_ref]

        exclude_query = adata_query.obs[cell_type_label].isin(exclude_bool)!= True
        adata_query = adata_query[exclude_query]
        embedding_query = embedding_query[exclude_query]

        knn = KNeighborsClassifier(n_neighbors=5)
        knn.fit(embedding_ref, adata_ref.obs[cell_type_label].tolist())
        cat_preds = knn.predict(embedding_query)

        f1_macro = f1_score(adata_query.obs[cell_type_label].tolist(), cat_preds, average='macro')
        acc = accuracy_score(adata_query.obs[cell_type_label].tolist(), cat_preds)

        accs.append(acc)
        f1_macros.append(f1_macro)

    print(f'{dataset} avg Accuracy {np.mean(accs)}, F1 macro {np.mean(f1_macro)}')
    print(f'{dataset} std Accuracy {np.std(accs)}, F1 macro {np.std(f1_macro)}')
