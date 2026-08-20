#!/bin/bash
#SBATCH --job-name=nequip_compide
#SBATCH --time=0-10:30:00
#SBATCH --output=results.out.%j
#SBATCH --error=results.err.%j
#SBATCH --nodes=1
##SBATCH --gres=gpu:2g.20gb:1
#SBATCH --gpus=1
 
eval "$(/home/lwalterb/.local/bin/micromamba shell hook --shell bash)"
micromamba activate hdp_ehulls

nequip-compile nequip.net:mir-group/NequIP-OAM-XL:0.1 mir-group__NequIP-OAM-XL__0.1.nequip.pt2 --mode aotinductor --device cuda --target ase
