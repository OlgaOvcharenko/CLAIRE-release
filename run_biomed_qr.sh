#!/bin/sh
#SBATCH -o logs/clair-%j.out
#SBATCH --nodes=1
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16G

# TODO Put filter to False
python main_optim.py --dname 'ImmuneAtlas_our_ref' --dname_query 'ImmuneAtlas_our_query' --n_repeat 2 --moco_k 2048 --moco_m 0.999 --moco_t 0.07 \
                    --block_level 1 --lat_dim 128 --symmetric True \
                    --select_hvg 4000 \
                    --knn 10 --alpha 0.5 --augment_set 'int' \
                    --lr 1e-4 --optim Adam --weight_decay 1e-5 --epochs 100 --batch_size 256\
                    --anchor_schedule 5 --fltr 'gmm' --yita 0.5 \
                    --save_freq 10 --start_epoch 0 --workers 6 --init 'uniform' \
                    --visualize_ckpts 100


# python main_optim.py --dname 'Pancreas_Mutaro' --dname_query 'Pancreas_Mutaro_query' --n_repeat 1 --moco_k 2048 --moco_m 0.999 --moco_t 0.07 \
#                      --block_level 1 --lat_dim 128 --symmetric True \
#                      --select_hvg 4000 \
#                      --knn 10 --alpha 0.5 --augment_set 'int' \
#                      --lr 1e-4 --optim Adam --weight_decay 1e-5 --epochs 100 --batch_size 256\
#                      --anchor_schedule 2 --fltr 'gmm' --yita 0.5 \
#                      --save_freq 10 --start_epoch 0 --workers 6 --init 'uniform' \
#                      --visualize_ckpts 100  

# python main_optim.py --dname 'Pancreas_Segerstolpe' --dname_query 'Pancreas_Segerstolpe_query' --n_repeat 1 --moco_k 2048 --moco_m 0.999 --moco_t 0.07 \
#                      --block_level 1 --lat_dim 128 --symmetric True \
#                      --select_hvg 4000 \
#                      --knn 10 --alpha 0.5 --augment_set 'int' \
#                      --lr 1e-4 --optim Adam --weight_decay 1e-5 --epochs 100 --batch_size 256\
#                      --anchor_schedule 2 --fltr 'gmm' --yita 0.5 \
#                      --save_freq 10 --start_epoch 0 --workers 6 --init 'uniform' \
#                      --visualize_ckpts 100  

# python main_optim.py --dname 'Pancreas_Wang' --dname_query 'Pancreas_Wang_query' --n_repeat 1 --moco_k 2048 --moco_m 0.999 --moco_t 0.07 \
#                      --block_level 1 --lat_dim 128 --symmetric True \
#                      --select_hvg 4000 \
#                      --knn 10 --alpha 0.5 --augment_set 'int' \
#                      --lr 1e-4 --optim Adam --weight_decay 1e-5 --epochs 100 --batch_size 256\
#                      --anchor_schedule 2 --fltr 'gmm' --yita 0.5 \
#                      --save_freq 10 --start_epoch 0 --workers 6 --init 'uniform' \
#                      --visualize_ckpts 100  

# python main_optim.py --dname 'Pancreas_Xin' --dname_query 'Pancreas_Xin_query' --n_repeat 1 --moco_k 2048 --moco_m 0.999 --moco_t 0.07 \
#                      --block_level 1 --lat_dim 128 --symmetric True \
#                      --select_hvg 4000 \
#                      --knn 10 --alpha 0.5 --augment_set 'int' \
#                      --lr 1e-4 --optim Adam --weight_decay 1e-5 --epochs 100 --batch_size 256\
#                      --anchor_schedule 2 --fltr 'gmm' --yita 0.5 \
#                      --save_freq 10 --start_epoch 0 --workers 6 --init 'uniform' \
#                      --visualize_ckpts 100  
