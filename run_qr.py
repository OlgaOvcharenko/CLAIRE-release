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

only_files = ['results1/', 'results2/']
datasets = ['Pancreas_Xin'] # 'Pancreas_Mutaro', 'Pancreas_Segerstolpe', 'Pancreas_Wang', 'Pancreas_Xin'


cell_type_label = "celltype"
batch = "batch"

for dataset in datasets:
    accs, f1_macros = [], []
    accs_pca, f1_macros_pca = [], []
    for folder in only_files:
        adata_ref = sc.read_h5ad(f'/cluster/home/oovcharenko/Olga_Data/{dataset}_qr.h5ad')
        adata_query = sc.read_h5ad(f'/cluster/home/oovcharenko/Olga_Data/{dataset}_qr_query.h5ad')

        embedding_ref = sc.read_h5ad("outputs/" + f"{dataset}/{folder}/ad_100_ref.h5ad").X
        embedding_query = sc.read_h5ad("outputs/" + f"{dataset}/{folder}/ad_100_query.h5ad").X
        
        # exclude_bool = list(set(
        #     list(set(adata_query.obs[cell_type_label].tolist()) - set(adata_ref.obs[cell_type_label].tolist())) + 
        #     list(set(adata_ref.obs[cell_type_label].tolist()) - set(adata_query.obs[cell_type_label].tolist()))))
        # exclude_ref = adata_ref.obs[cell_type_label].isin(exclude_bool)!= True
        # adata_ref = adata_ref[exclude_ref]
        # embedding_ref = embedding_ref[exclude_ref]

        # exclude_query = adata_query.obs[cell_type_label].isin(exclude_bool)!= True
        # adata_query = adata_query[exclude_query]
        # embedding_query = embedding_query[exclude_query]

        knn = KNeighborsClassifier(n_neighbors=5)
        knn.fit(embedding_ref, adata_ref.obs[cell_type_label].tolist())
        cat_preds = knn.predict(embedding_query)
        cell_types_list = pd.unique(adata_ref.obs[cell_type_label]).tolist()
        acc = accuracy_score(adata_query.obs[cell_type_label].tolist(), cat_preds)
        f1_macro = f1_score(adata_query.obs[cell_type_label].tolist(), cat_preds, labels=cell_types_list, average='macro')

        # print(f1_score(adata_query.obs[cell_type_label].tolist(), cat_preds, labels=cell_types_list, average=None))
        
        accs.append(acc)
        f1_macros.append(f1_macro)

    print(f'{dataset} avg Accuracy {np.mean(accs)}, F1 macro {np.mean(f1_macro)}')
    print(f'{dataset} std Accuracy {np.std(accs)}, F1 macro {np.std(f1_macro)}')

exit()

only_files = ['results1/', 'results2/']
datasets = ['ImmuneAtlas']

cell_type_label = "cell_type"
batch = "batchlb"

for dataset in datasets:
    accs, f1_macros, f1_weighteds = [], [], []
    for folder in only_files:
        adata_ref = sc.read_h5ad(f'/cluster/home/oovcharenko/Olga_Data/{dataset}_qr.h5ad')
        adata_query = sc.read_h5ad(f'/cluster/home/oovcharenko/Olga_Data/{dataset}_qr_query.h5ad')
            
        embedding_ref = sc.read_h5ad("outputs/" + f"{dataset}/{folder}/ad_100_ref.h5ad").X
        embedding_query = sc.read_h5ad("outputs/" + f"{dataset}/{folder}/ad_100_query.h5ad")

        # exclude_bool = list(set(
        #     list(set(adata_query.obs[cell_type_label].tolist()) - set(adata_ref.obs[cell_type_label].tolist())) + 
        #     list(set(adata_ref.obs[cell_type_label].tolist()) - set(adata_query.obs[cell_type_label].tolist()))))
        # exclude_ref = adata_ref.obs[cell_type_label].isin(exclude_bool)!= True
        # adata_ref = adata_ref[exclude_ref]
        # embedding_ref = embedding_ref[exclude_ref]

        # exclude_query = adata_query.obs[cell_type_label].isin(exclude_bool)!= True
        # adata_query = adata_query[exclude_query]
        # embedding_query = embedding_query[exclude_query]

        knn = KNeighborsClassifier(n_neighbors=5)
        knn.fit(embedding_ref, adata_ref.obs[cell_type_label].tolist())
        cat_preds = knn.predict(embedding_query)

        cell_types_list = pd.unique(adata_ref.obs[cell_type_label]).tolist()
        acc = accuracy_score(adata_query.obs[cell_type_label].tolist(), cat_preds)
        f1_macro = f1_score(adata_query.obs[cell_type_label].tolist(), cat_preds, labels=cell_types_list, average='macro')
        
        accs.append(acc)
        f1_macros.append(f1_macro)

    print(f'{dataset} avg Accuracy {np.mean(accs)}, F1 macro {np.mean(f1_macro)}')
    print(f'{dataset} std Accuracy {np.std(accs)}, F1 macro {np.std(f1_macro)}')
