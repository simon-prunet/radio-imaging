# radio-imaging

This repository contains the implementation and results for the paper entitled "Multi-Step Deconvolution of Radio-Interferometric Images". The code is implemented in Julia, and integrated in rascil. These are not stored directly on this repository but rather in a variety of others, specifically:
1. Julia code: https://github.com/andferrari/DeconvMultiStep.jl
2. RASCIL integration: https://gitlab.com/prunet1/ska-sdp-func-python and https://gitlab.com/igorawratu/rascil-main

## Instructions for running system
The entire system is contained within a docker image, for which the Dockerfile is within this repository. To run: 
1. Clone repository
2. Install docker
3. Build docker image using the command > docker build -t <image_name> .
4. Run the image using > docker run -p 8888:8888 -v <host_directory>:<container_directory> -it <image_name>
5. The docker container will start. You can then run RASCIL from its home directory, as all the required julia files will be copied there.
6. Example commands for execution, datasets can be found in the various data directories. Measurement sets will have to be copied to the home directory, or the full path need to be used:
    - Single-step LASSO: > python /rascil-main/rascil/apps/rascil_imager.py --ingest_msname SGRA.ms --ingest_vis_nchan 1 --imaging_npixel 512 --imaging_cellsize 0.00001849451 --imaging_weighting uniform --clean_nmajor 5 --clean_algorithm mstep --clean_fractional_threshold 0.3 --clean_threshold 1e-3 --clean_restored_output integrated --mstep_output_intermediate True --mstep_mode full --mstep_lambda 0.01 --mstep_lambda_mul 1 --mstep_wavelet daubechies
    - Low-resolution step: > python /rascil-main/rascil/apps/rascil_imager.py --ingest_msname SGRA_small_baselines_60.ms --ingest_vis_nchan 1 --imaging_npixel 512 --imaging_cellsize 0.00001849451 --imaging_weighting uniform --clean_nmajor 5 --clean_algorithm mstep --clean_fractional_threshold 0.3 --clean_threshold 1e-3 --clean_restored_output integrated --mstep_output_intermediate True --mstep_mode low --mstep_lambda 0.05 --mstep_lambda_mul 2 --mstep_wavelet daubechies --mstep_cut_center 55 --mstep_cut_hw 5
    - Full-resolution step: > python /rascil-main/rascil/apps/rascil_imager.py --ingest_msname SGRA_long_baselines_50.ms --ingest_vis_nchan 1 --imaging_npixel 512 --imaging_cellsize 0.00001849451 --imaging_weighting uniform --clean_nmajor 5 --clean_algorithm mstep --clean_fractional_threshold 0.3 --clean_threshold 1e-3 --clean_restored_output integrated --mstep_output_intermediate True --mstep_mode multi-step --mstep_lambda 0.05 --mstep_lambda_mul 2 --mstep_lowfreq_image small_baselines_dataset_deconvolved.fits --mstep_wavelet daubechies --mstep_cut_center 55 --mstep_cut_hw 5

## Splitting measurement sets
Measurement sets can be split using casa, using the uvrange parameter of the split function. A ready-made docker for this can be found at rxastro/casa6. Check README.casatools for the exact commands to perform the split.

## Results
The results for the paper are presented in a variety of jupyter notebooks. These are:
- dirty.ipynb for presenting the original datasets and their split
- lambda.ipynb for preliminary experiments on $\lambda$
- noise.ipynb for preliminary experiments and estimation of $\sigma^2$ and $\eta^2$
- split.ipynb for results regarding partition configurations, and against the single-step LASSO reconstruction

The raw data for our results can be accessed in the results/ directory