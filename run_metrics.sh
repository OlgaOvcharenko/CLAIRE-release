#!/bin/sh
#SBATCH -o logs/clear-%j.out
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=15:00:00
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=END,FAIL

python run_metrics.py
