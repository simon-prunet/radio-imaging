#!/bin/bash
#SBATCH --job-name=l1large_test
#SBATCH --partition=x40
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=l1large_test.out
#SBATCH --error=l1large_test.err 
#SBATCH --exclusive

source ~/.bashrc

module load intel-oneapi/2023.1.0
module load intel-oneapi/2023.1.0-gcc13
module load anaconda3/2022.10
module load julia/1.8.3

conda activate test2

/usr/bin/time -v python test_deconv.py
