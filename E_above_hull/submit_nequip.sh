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
mlip="Nequip-OAM-XL"
export JOBFLOW_CONFIG_FILE="/home/lwalterb/hdp_ehull/jobflow_settings/nequip_jobflow.yaml"
micromamba activate ehull_nequip


export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

N_CONCURRENT=4  # start here, watch nvidia-smi + squeue timing
for i in $(seq 0 $((N_CONCURRENT-1))); do
    python drive_nequip_calcs.py  $i  $N_CONCURRENT $mlip &
    sleep 2 
   # python drive_nequip_calcs.py  0 3  $mlip &
done
wait
