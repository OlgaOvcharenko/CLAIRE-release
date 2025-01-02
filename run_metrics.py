import argparse
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import scanpy as sc
import anndata as ad
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from os import listdir
import os
from os.path import isfile, join
import pandas as pd
from os import walk


_BIO_METRICS = BioConservation(isolated_labels=True, 
                               nmi_ari_cluster_labels_leiden=True, 
                               nmi_ari_cluster_labels_kmeans=False, 
                               silhouette_label=True, 
                               clisi_knn=True
                               )
_BATCH_METRICS = BatchCorrection(graph_connectivity=True, 
                                 kbet_per_label=True, 
                                 ilisi_knn=True, 
                                 pcr_comparison=True, 
                                 silhouette_batch=True
                                 )


def evaluate_model(adata, batch_key="batch", cell_type_label="cell_type_l1"):
    names_obs = ['new']
    print(adata)
    bm = Benchmarker(
                adata,
                batch_key=batch_key,
                label_key=cell_type_label,
                embedding_obsm_keys=names_obs,
                bio_conservation_metrics=_BIO_METRICS,
                batch_correction_metrics=_BATCH_METRICS,
                n_jobs=4,
            )
    bm.benchmark()
    a = bm.get_results(False, True)
    results = a.round(decimals=4)
    return results

def scale_result(df):
    cols = ['Isolated labels', 'Leiden NMI', 'Leiden ARI', 'Silhouette label', 'cLISI', 'Silhouette batch', 'iLISI', 'KBET', 'Graph connectivity', 'PCR comparison']
    df[cols] = MinMaxScaler().fit_transform(df[cols])
    scaled= pd.DataFrame(df, columns=df.columns, index=df.index)
    
    biometrics = [i for i in ['Isolated labels', 'Leiden NMI', 'Leiden ARI', 'Silhouette label', 'cLISI']]
    batchmetrics = [i for i in ['Silhouette batch', 'iLISI', 'KBET', 'Graph connectivity', 'PCR comparison']]
    scaled[f"Batch correction final"] = scaled[batchmetrics].mean(1)
    scaled[f"Bio conservation final"] = scaled[biometrics].mean(1)
    scaled[f"Total final"] = 0.6 * scaled[f"Bio conservation final"] + 0.4 * scaled[f"Batch correction final"]
    df[f"Bio conservation final"] = scaled[f"Bio conservation final"].copy()
    df[f"Batch correction final"] = scaled[f"Batch correction final"].copy()
    df[f"Total final"] = scaled[f"Total final"].copy()
    return df

repeat = 0


# only_files = ['result/CLEAR/', 'result/CLEAR_0/']
only_files = ['outputs/']

for mypath in only_files:
    #filenames = next(walk(mypath), (None, None, []))[2]  # [] if no file 
    filenames = [os.path.join(dp, f) for dp, dn, filenames in os.walk(mypath) for f in filenames if os.path.splitext(f)[1] == '.h5ad']
    
    # filenames = ["outputs/ImmHuman_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[2]_filter=gmm_yita=0.5_eps=120_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results1/ad_20.h5ad", "outputs/ImmHuman_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[2]_filter=gmm_yita=0.5_eps=120_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results1/ad_10.h5ad", "outputs/ImmHuman_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[2]_filter=gmm_yita=0.5_eps=120_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results1/ad_120.h5ad", "outputs/ImmHuman_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[2]_filter=gmm_yita=0.5_eps=120_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results1/ad_40.h5ad", "outputs/ImmHuman_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[2]_filter=gmm_yita=0.5_eps=120_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results1/ad_80.h5ad", "outputs/ImmHuman_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[2]_filter=gmm_yita=0.5_eps=120_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results2/ad_20.h5ad", "outputs/ImmHuman_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[2]_filter=gmm_yita=0.5_eps=120_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results2/ad_10.h5ad", "outputs/ImmHuman_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[2]_filter=gmm_yita=0.5_eps=120_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results2/ad_120.h5ad", "outputs/ImmHuman_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[2]_filter=gmm_yita=0.5_eps=120_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results2/ad_40.h5ad", "outputs/ImmHuman_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[2]_filter=gmm_yita=0.5_eps=120_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results2/ad_80.h5ad", "outputs/PBMC_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[10]_filter=gmm_yita=0.5_eps=100_lr=0.0001_batch-size=256_adjustLr=True_schedule=[10]/results1/ad_100.h5ad", "outputs/PBMC_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[10]_filter=gmm_yita=0.5_eps=100_lr=0.0001_batch-size=256_adjustLr=True_schedule=[10]/results2/ad_100.h5ad", "outputs/Lung_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[5]_filter=gmm_yita=0.5_eps=100_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results1/ad_100.h5ad", "outputs/Lung_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[5]_filter=gmm_yita=0.5_eps=100_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results2/ad_100.h5ad", "outputs/MCA_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[4]_filter=gmm_yita=0.5_eps=100_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results1/ad_100.h5ad", "outputs/MCA_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[4]_filter=gmm_yita=0.5_eps=100_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results2/ad_100.h5ad", "outputs/Pancreas_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[2]_filter=gmm_yita=0.5_eps=100_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results1/ad_100.h5ad", "outputs/Pancreas_our/knn=10_alpha=0.5_augSet=['int']_anc-schedl=[2]_filter=gmm_yita=0.5_eps=100_lr=0.0001_batch-size=256_adjustLr=False_schedule=[10]/results2/ad_100.h5ad"]

    for file in filenames:
        data = file.split("/")[1].split("_")[0]

        if data == "ImmuneAtlas":
            cell_type_label = "cell_type"
            batch = "batchlb"
            adata_RNA = sc.read_h5ad('/cluster/home/oovcharenko/Olga_Data/ImmuneAtlas.h5ad')

        elif data == "Lung":
            cell_type_label = "cell_type"
            batch = "batch"
            adata_RNA = sc.read_h5ad('/cluster/home/oovcharenko/Olga_Data/Lung.h5ad')

        elif data == "MCA":
            cell_type_label = "CellType"
            batch = "batch"
            adata_RNA = sc.read_h5ad('/cluster/home/oovcharenko/Olga_Data/MCA.h5ad')

        elif data == "Pancreas":
            cell_type_label = "celltype"
            batch = "batch"
            adata_RNA = sc.read_h5ad('/cluster/home/oovcharenko/Olga_Data/Pancreas.h5ad')

        elif data == "PBMC":
            cell_type_label = "CellType"
            batch = "batch"
            adata_RNA = sc.read_h5ad('/cluster/home/oovcharenko/Olga_Data/PBMC.h5ad')

        elif data == "ImmHuman":
            cell_type_label = "CellType"
            batch = "batch"
            adata_RNA = sc.read_h5ad('/cluster/home/oovcharenko/Olga_Data/ImmHuman.h5ad')
            
        else:
            continue

        print(data)
        
        adata_RNA.obsm["new"] = sc.read_h5ad(file).X

        df = evaluate_model(adata=adata_RNA, batch_key=batch, cell_type_label=cell_type_label)
        
        print(df)
        if file.split("/")[3] == 'results1':
            df.to_csv(f'results/bc1/{data}_unscaled.csv')
        else:
            df.to_csv(f'results/bc2/{data}_unscaled.csv')

