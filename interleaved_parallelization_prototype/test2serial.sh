#!/bin/bash
#SBATCH --job-name=interleaved_serial_test
#SBATCH --partition=x40
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=test2serial.out
#SBATCH --error=test2serial.err 
#SBATCH --exclusive

module load intel-oneapi/2023.1.0
module load gcc/12.2.0
module load openmpi/4.1.5
module load anaconda3/2022.10
module load julia

source ~/.bashrc

conda activate test2

python ip_serial.py
