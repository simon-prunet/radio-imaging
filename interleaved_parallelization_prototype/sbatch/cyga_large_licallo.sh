#!/bin/bash
#SBATCH --job-name=cygalarge_msc
#SBATCH --partition=x40
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=cygalarge_msc.out
#SBATCH --error=cygalarge_msc.err 
#SBATCH --exclusive

source ~/.bashrc

module load intel-oneapi/2023.1.0
module load intel-oneapi/2023.1.0-gcc13
module load anaconda3/2022.10
module load julia/1.8.3

conda activate test2

python ip_msclean.py configs/cyg_a_large_test.config
