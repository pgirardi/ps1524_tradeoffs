#!/bin/bash 
#SBATCH --job-name=dRep-syr
#SBATCH --cpus-per-task=32 
#SBATCH --mem=100G 
#SBATCH -A mpag-np 
#SBATCH -p mpag-np 
#SBATCH --time=168:00:00 
#SBATCH --output=dRep-syr_%j.out 
#SBATCH --error=dRep-syr_%j.err 

set -euo pipefail

# Always start from the submission directory 
cd "$SLURM_SUBMIT_DIR" 
ENV_DIR="$SLURM_SUBMIT_DIR/drep_env" 

# Use the env's binaries without activating 
export PATH="$ENV_DIR/bin:$PATH" 

# Run dRep (skip CheckM) 
dRep dereplicate out_syr\ 
    -g syringae/*.fna \ 
    --ignoreGenomeQuality \ 
    -pa 0.95 -sa 0.999 -nc 0.50 \ 
    -p "$SLURM_CPUS_PER_TASK"