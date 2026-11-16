# Overview
The code in this directory pertains to a radio-interferometric prototype that parallelizes by spatial frequency, and is described in papers:
Wang, S., Mignot, S., Prunet, S., Di Mascolo, L., Spinelli, M. and Ferrari, A., 2025. A Decentralized Framework for Radio-interferometric Image Reconstruction. The Astronomical Journal 169.6 (2025): 289 and Wang, S., Prunet, S., Mignot, S. and Ferrari, A., 2026. SMURFIT: Spatial-frequency Multi-class Radio Fourier Imaging Technique. Astronomy and Computing (2026): 101212.

# A Decentralized Framework for Radio-interferometric Image Reconstruction
The code for this specific paper uses several external repositories. Specifically, for the L1 method, we use a Julia implementation available at https://github.com/andferrari/DeconvMultiStep.jl. For other aspects, such as de/gridding, multi-scale CLEAN etc. we use RASCIL 1.1.0. We also use mpi4pyto distribute the system across multiple processes/nodes. All these dependencies should be handled with the setup_env.sh file in the root directory

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

# SMURFIT: Spatial-frequency Multi-class Radio Fourier Imaging Technique
The code for this paper is relies on very similar external repositiries, with the exception of the julia code, which has now been moved into the python library.

## Code structure
There are four different executable python scripts, each running a different method:
1. ip_l1_p_multipartition.py, which implements a parallelized version of the L1 method. This is distributed using mpi4py, thus, to run, you would execute, for k-1 partitions: 
```
mpirun -n <k> python ip_l1_p_multipartition.py <config_file>
```
where the config_file is a json file that describes the ingested dataset and parameters. There are several examples of these files in this directory. This differs from the previous as it can run on an arbitrary number of nodes

2. ip_l1.py, which implements a serial version of the L1 method. This differs from ip_serial.py as a different regularization strategy is used, and should work better with less parameter tuning. This is run using:
```
python ip_l1.py <config_file>
```

3. ip_msclean_p_multipartition.py, which implements a parallelized version of multi-scale CLEAN. This is run using:
```
mpirun -n <k> python ip_msclean_p_multipartition.py <config_file>
```
Like the above, this also runs on an aribitrary number of partitions

4. ip_msclean.py, which implements the serial multi-scale CLEAN and is run using:
```
python ip_msclean_mp.py <config_file>
```
Largely the same as ip_msclean.py, just has some slightly different orchestration for the new paper.

## Datasets
Our experiment datasets can be downloaded using the get_data.sh script in the data8 directory. Parameters for their reconstruction can be found in the .config files which are downloaded together with the datasets. The slurm batch files for obtaining our results are also available with the datasets in the form of l1.sh, msc.sh, pl1.sh, pmsc.sh

## Results
Our results can be found in the following notebooks under the directory notebooks_and_pyscripts: 
- multipartition_bda2_results.ipynb and multipartition_bda_realdata_results.ipynb have the results for the simulated and real datasets
- bda.ipynb performs preliminary tests on our bda implementation
- filters_prototype.ipynb shows some 1d preliminary experiments with filters
- loadbal.ipynb shows shows results for our load-balancing algorithm
- sim_largems.ipynb simulates a large measurement set, performing bda on it before writing.