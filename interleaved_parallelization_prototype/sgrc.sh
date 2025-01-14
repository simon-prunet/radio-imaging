#!/bin/bash
#SBATCH --job-name=sgrc
#SBATCH --partition=cpu_p1
#SBATCH --time=2:00:00
#SBATCH --ntasks=3
#SBATCH --ntasks-per-node=1
#SBATCH --output=sgrc.out
#SBATCH --error=sgrc.err 
#SBATCH --exclusive
#SBATCH -A klu@cpu

module load intel-oneapi-all/2023.1
module load gcc/12.2.0
module load anaconda-py3/2024.06
module load julia/1.10.4
conda activate ri_parallelization

mpirun -np 3 python ip_prototype.py sgrc.config
