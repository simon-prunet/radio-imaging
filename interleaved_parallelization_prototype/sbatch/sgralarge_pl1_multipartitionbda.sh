#!/bin/bash
#SBATCH --job-name=sgralarge_pl1bda
#SBATCH --partition=cpu_p1
#SBATCH --time=20:00:00
#SBATCH --ntasks=6
#SBATCH --ntasks-per-node=1
#SBATCH --output=sgralarge_pl1bda.out
#SBATCH --error=sgralarge_pl1bda.err 
#SBATCH --exclusive
#SBATCH -A klu@cpu

source ~/.bashrc

module load intel-oneapi-all/2023.1
module load gcc/12.2.0
module load anaconda-py3/2024.06
module load julia/1.10.4
conda activate /lustre/fswork/projects/rech/klu/ulc65eb/conda_envs/ri_parallelization2

export JULIA_DEPOT_PATH=/lustre/fswork/projects/rech/klu/ulc65eb/libs/julia

mpirun -np 6 python ip_l1_p_multipartition.py configs/sgra_multipartition_bda.config
