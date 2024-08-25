#!/bin/bash

#SBATCH --account=mel_333
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --array=1-2000

module purge
eval "$(conda shell.bash hook)"
conda activate /project/mel_333/crjEnv

cd /project/mel_333/crjCode

python slurm-test.py
