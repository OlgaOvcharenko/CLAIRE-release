#!/bin/sh
#SBATCH -o logs/clair-%j.out
#SBATCH --nodes=1
#SBATCH -p gpu
#SBATCH --gres=gpu:rtx3090:1
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-type=END,FAIL

# mca
#python main_optim.py --dname 'MCA_our' --n_repeat 2 --moco_k 2048 --moco_m 0.999 --moco_t 0.07 \
#                    --block_level 1 --lat_dim 128 --symmetric True \
#                    --select_hvg 4000 \
#                    --knn 10 --alpha 0.5 --augment_set 'int' \
#                    --anchor_schedule 4 --fltr 'gmm' --yita 0.5 \
#                    --lr 1e-4 --optim Adam --weight_decay 1e-5 --epochs 100 --batch_size 256\
#                    --save_freq 10 --start_epoch 0 --workers 6 --init 'uniform' \
#                    --visualize_ckpts 100

# # pbmc
#python main_optim.py --dname 'PBMC_our' --n_repeat 2 --moco_k 2048 --moco_m 0.999 --moco_t 0.07 \
#                     --block_level 1 --lat_dim 128 --symmetric True \
#                     --select_hvg 4000 \
#                     --knn 10 --alpha 0.5 --augment_set 'int' \
#   		     --lr 1e-4 --adjustLr --schedule 10 --optim Adam --weight_decay 1e-5 --epochs 100 --batch_size 256\
#                     --anchor_schedule 10 --fltr 'gmm' --yita 0.5 \
#                     --save_freq 10 --start_epoch 0 --workers 6 --init 'uniform' \
#                     --visualize_ckpts 100

# # pancreas
# python main_optim.py --dname 'Pancreas_our' --n_repeat 1 --moco_k 2048 --moco_m 0.999 --moco_t 0.07 \
#                      --block_level 1 --lat_dim 128 --symmetric True \
#                      --select_hvg 4000 \
#                      --knn 10 --alpha 0.5 --augment_set 'int' \
#                      --lr 1e-4 --optim Adam --weight_decay 1e-5 --epochs 100 --batch_size 256\
#                      --anchor_schedule 2 --fltr 'gmm' --yita 0.5 \
#                      --save_freq 10 --start_epoch 0 --workers 6 --init 'uniform' \
#                      --visualize_ckpts 99   

# # immune
# python main.py --dname 'ImmHuman_our' --n_repeat 2 --moco_k 2048 --moco_m 0.999 --moco_t 0.07 \
#                     --block_level 1 --lat_dim 128 --symmetric True \
#                     --select_hvg 4000 \
#                     --knn 10 --alpha 0.5 --augment_set 'int' \
#                     --lr 1e-4 --optim Adam --weight_decay 1e-5 --epochs 120 --batch_size 256\
#                     --anchor_schedule 2 --fltr 'gmm' --yita 0.5 \
#                     --save_freq 10 --start_epoch 0 --workers 6 --init 'uniform' \
#                     --visualize_ckpts 10 20 40 80 120 

# lung
#python main.py --dname 'Lung_our' --n_repeat 2 --moco_k 2048 --moco_m 0.999 --moco_t 0.07 \
#                     --block_level 1 --lat_dim 128 --symmetric True \
#                     --select_hvg 4000 \
#                     --knn 10 --alpha 0.5 --augment_set 'int' \
#                    --lr 1e-4 --optim Adam --weight_decay 1e-5 --epochs 100 --batch_size 256\
#                     --anchor_schedule 5 --fltr 'gmm' --yita 0.5 \
#                     --save_freq 10 --start_epoch 0 --workers 6 --init 'uniform' \
#                     --visualize_ckpts 100

# immune atlas
python main.py --dname 'ImmuneAtlas_our' --n_repeat 2 --moco_k 2048 --moco_m 0.999 --moco_t 0.07 \
                    --block_level 1 --lat_dim 128 --symmetric True \
                    --select_hvg 4000 \
                    --knn 10 --alpha 0.5 --augment_set 'int' \
                    --lr 1e-4 --optim Adam --weight_decay 1e-5 --epochs 100 --batch_size 256\
                    --anchor_schedule 5 --fltr 'gmm' --yita 0.5 \
                    --save_freq 10 --start_epoch 0 --workers 6 --init 'uniform' \
                    --visualize_ckpts 100

