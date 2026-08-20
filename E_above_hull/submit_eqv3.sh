#!/bin/bash
#SBATCH --job-name=hdp_ehull_EQv3
#SBATCH --time=7-00:00:00
#SBATCH --output=results.out.%j
#SBATCH --error=results.err.%j
#SBATCH --nodes=1
##SBATCH --gres=gpu:4g.40gb:1
#SBATCH --cpus-per-task=4
#SBATCH --gpus=1

 
eval "$(/home/lwalterb/.local/bin/micromamba shell hook --shell bash)"
mlip="Equiformer_v3"
export JOBFLOW_CONFIG_FILE="/home/lwalterb/hdp_ehull/jobflow_settings/eqv3_jobflow.yaml"
micromamba activate ehull_eqv3


export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python drive_eqv3_calcs.py 0 1 $mlip

#N_CONCURRENT=6  # start here, watch nvidia-smi + squeue timing
#for i in $(seq 0 $((N_CONCURRENT-1))); do
#	python drive_ehull_calcs.py  $i  $N_CONCURRENT $mlip &
#done
#wait
