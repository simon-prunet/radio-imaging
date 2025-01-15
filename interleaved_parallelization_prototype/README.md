# Parallelization Prototype

The code in this directory pertains to a radio-interferometric prototype that parallelizes by baseline length, and is described in the paper:
Wang, S., Mignot, S., Prunet, S., Di Mascolo, L., Spinelli, M. and Ferrari, A., 2025. A Decentralized Framework for Radio-Interferometric Image Reconstruction. (tbd)

We make use of code from several external repositories. Specifically, for the L1 method, we use a Julia implementation available at https://github.com/andferrari/DeconvMultiStep.jl. For other aspects, such as de/gridding, multi-scale CLEAN etc. we use RASCIL 1.1.0, available at https://gitlab.com/ska-telescope/external/rascil-main with the branch tag 1.1.0. We also use mpi4py (https://mpi4py.readthedocs.io/en/stable/) to distribute the system across multiple processes/nodes.

## Code structure

There are four different executable python scripts, each running a different method:
1. ip_prototype.py, which implements a parallelized version of the L1 method. This is distributed using mpi4py, thus, to run, you would execute: 
```
mpirun -n 3 python ip_prototype.py <config_file>
```
where the config_file is a json file that describes the ingested dataset and parameters. There are several examples of these files in this directory.

2. ip_serial.py, which implements a serial version of the L1 method. This is run using:
```
python ip_serial.py <config_file>
```

3. ip_msclean_p.py, which implements a parallelized version of multi-scale CLEAN. This is run using:
```
mpirun -n 3 python ip_msclean_p.py <config_file>
```

4. ip_msclean.py, which implements the serial multi-scale CLEAN and is run using:
```
python ip_msclean.py <config_file>
```

## Dataset parameters
The parameters for all of our datasets can be found in their respective (non "test") config files. The datasets are stored in the data6/ directory of this repository. Running the get_data.sh bash script will automatically download and extract them here.

## Slurm sbatch files
We executed our code on the Jean Zay cluster, which uses the Slurm to distribute jobs. Examples of sbatch files used for our results can be found in this directory.

## Results
Our results are presented in two jupyter notebooks, under the directory notebooks_and_pyscripts. The first is in ri_split_results.ipynb, which shows results for our paper. The second is in parallel_single_dataset.ipynb, which performs the same experiments but in a reduced environment, freezing all dataset parameters except for the dataset size across three datasets. To run the notebooks, you will need to download the actual results, which can be done by running the get_results.sh script under the results directory.