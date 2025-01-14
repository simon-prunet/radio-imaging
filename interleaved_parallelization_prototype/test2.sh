#!/bin/bash
#SBATCH --job-name=interleaved_test
#SBATCH --partition=x40
#SBATCH --time=10:00:00
#SBATCH --ntasks=3
#SBATCH --ntasks-per-node=1
#SBATCH --output=test2.out
#SBATCH --error=test2.err 
#SBATCH --exclusive

module load intel-oneapi/2023.1.0
module load gcc/12.2.0
module load openmpi/4.1.5
module load anaconda3/2022.10
module load julia

source ~/.bashrc

conda activate test2

mpirun -np 3 python ip_prototype.py
