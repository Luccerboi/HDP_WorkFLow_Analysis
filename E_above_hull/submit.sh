#!/bin/bash
#SBATCH --job-name=hdp_ehull_MACE
#SBATCH --time=2-10:30:00
#SBATCH --output=results.out.%j
#SBATCH --error=results.err.%j
#SBATCH --nodes=1
##SBATCH --gres=gpu:4g.40gb:1
#SBATCH --cpus-per-task=16
#SBATCH --gpus=1
 
eval "$(/home/lwalterb/.local/bin/micromamba shell hook --shell bash)"
mlip="Nequip"
micromamba activate ehull_nequip

export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

N_CONCURRENT=1  # start here, watch nvidia-smi + squeue timing
for i in $(seq 0 $((N_CONCURRENT-1))); do
    python drive_ehull_calcs.py  $i  $N_CONCURRENT $mlip &
done
wait
