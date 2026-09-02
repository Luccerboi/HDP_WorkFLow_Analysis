#!/bin/bash
#SBATCH --job-name=hdp_ehull_MACE
#SBATCH --time=7-00:00:00
#SBATCH --output=results.out.%j
#SBATCH --error=results.err.%j
#SBATCH --nodes=1
##SBATCH --gres=gpu:2g.20gb:1
#SBATCH --cpus-per-task=4
#SBATCH --gpus=1
 
eval "$(/home/lwalterb/.local/bin/micromamba shell hook --shell bash)"
mlip="GRACE-3L-OMAT-large-ft-AM"
export JOBFLOW_CONFIG_FILE="/home/lwalterb/hdp_ehull/jobflow_settings/grace_jobflow.yaml"
micromamba activate ehull_grace

export TF_USE_LEGACY_KERAS=1

export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

N_CONCURRENT=1  # start here, watch nvidia-smi + squeue timing
python drive_grace_calcs.py  0 1  $mlip &
wait
