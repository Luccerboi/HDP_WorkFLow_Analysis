#!/bin/bash
#SBATCH --job-name=test_Macemp0b3
#SBATCH --time=3-10:30:00
#SBATCH --output=test48.out.%j
#SBATCH --error=test48.err.%j
#SBATCH --nodes=1
##SBATCH --gres=gpu:1g.10gb:1
#SBATCH --gpus=1
 
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
eval "$(/home/lwalterb/.local/bin/micromamba shell hook --shell bash)"
micromamba activate hdp_ehull

python /home/lwalterb/hdp_ehull/run_ehull_calcs.py MACE-MPA-0 200 48
