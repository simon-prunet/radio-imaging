#!/bin/bash
#SBATCH --job-name=hltaularge_msc
#SBATCH --partition=x40
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=hltaularge_msc.out
#SBATCH --error=hltaularge_msc.err 
#SBATCH --exclusive

module load intel-oneapi/2023.1.0
module load intel-oneapi/2023.1.0-gcc13
module load anaconda3/2022.10
module load julia/1.8.3

conda activate test2

python ip_msclean.py hltau_large_test.config
