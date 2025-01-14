#!/bin/bash
#SBATCH --job-name=sgra_msc
#SBATCH --partition=cpu_p1
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=sgra_msc.out
#SBATCH --error=sgra_msc.err 
#SBATCH --exclusive
#SBATCH -A klu@cpu

module load intel-oneapi-all/2023.1
module load gcc/12.2.0
module load anaconda-py3/2024.06
module load julia/1.10.4
conda activate ri_parallelization

python ip_msclean.py sgra.config
