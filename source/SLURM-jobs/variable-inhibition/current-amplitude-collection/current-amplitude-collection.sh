#!/bin/bash

#SBATCH --account=mel_333
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=16G
#SBATCH --time=5:00:00
#SBATCH --array=1-240

module purge
eval "$(conda shell.bash hook)"
conda activate /project/mel_333/crjEnv

cd /project/mel_333/crjCode

python current-amplitude-collection.py
