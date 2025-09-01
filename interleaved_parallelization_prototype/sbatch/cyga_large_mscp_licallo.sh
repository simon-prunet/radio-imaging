#!/bin/bash
#SBATCH --job-name=cyga_large_mscp
#SBATCH --partition=x40
#SBATCH --time=24:00:00
#SBATCH --ntasks=3
#SBATCH --ntasks-per-node=1
#SBATCH --output=cyga_large_mscp.out
#SBATCH --error=cyga_large_mscp.err 
#SBATCH --exclusive

source ~/.bashrc

module load intel-oneapi/2023.1.0
module load openmpi/4.1.5/gcc-12.2.0
module load anaconda3/2022.10
module load julia/1.8.3

conda activate test2

mpiexec -np 3 python ip_msclean_p.py configs/cyg_a_large_test.config