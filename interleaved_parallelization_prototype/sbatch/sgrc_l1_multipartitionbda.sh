#!/bin/bash
#SBATCH --job-name=sgrcl1bda
#SBATCH --partition=cpu_p1
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=sgrcl1bda.out
#SBATCH --error=sgrcl1bda.err 
#SBATCH --exclusive
#SBATCH -A klu@cpu

source ~/.bashrc

module load intel-oneapi-all/2023.1
module load gcc/12.2.0
module load anaconda-py3/2024.06
module load julia/1.10.4
conda activate /lustre/fswork/projects/rech/klu/ulc65eb/conda_envs/ri_parallelization2

export JULIA_DEPOT_PATH=/lustre/fswork/projects/rech/klu/ulc65eb/libs/julia

python ip_l1_mp.py configs/sgrc_multipartition_bda.config
