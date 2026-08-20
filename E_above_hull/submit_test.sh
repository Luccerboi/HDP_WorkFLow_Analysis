#!/bin/bash
#SBATCH --job-name=hdp_ehull_MACE
#SBATCH --time=7-00:00:00
#SBATCH --output=results.out.%j
#SBATCH --error=results.err.%j
#SBATCH --nodes=1
#SBATCH --gres=gpu:2g.20gb:1
#SBATCH --cpus-per-task=4
##SBATCH --gpus=1
 
eval "$(/home/lwalterb/.local/bin/micromamba shell hook --shell bash)"
mlip="SevenNet-omat24"
export JOBFLOW_CONFIG_FILE="/home/lwalterb/hdp_ehull/jobflow_settings/7net_omat_jobflow.yaml"
micromamba activate ehull_sevennet


export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python run_ehull_batch.py $mlip 10 5 
