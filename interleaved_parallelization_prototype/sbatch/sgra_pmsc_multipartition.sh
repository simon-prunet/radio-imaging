#!/bin/bash
#SBATCH --job-name=sgrapmsc
#SBATCH --partition=cpu_p1
#SBATCH --time=01:00:00
#SBATCH --ntasks=6
#SBATCH --ntasks-per-node=1
#SBATCH --output=sgrapmsc.out
#SBATCH --error=sgrapmsc.err 
#SBATCH --exclusive
#SBATCH -A klu@cpu

source ~/.bashrc

module load intel-oneapi-all/2023.1
module load gcc/12.2.0
module load anaconda-py3/2024.06
module load julia/1.10.4
conda activate /lustre/fswork/projects/rech/klu/ulc65eb/conda_envs/ri_parallelization2

export JULIA_DEPOT_PATH=/lustre/fswork/projects/rech/klu/ulc65eb/libs/julia

mpirun -np 6 python ip_msclean_p_multipartition.py configs/sgra_multipartition_jz.config
