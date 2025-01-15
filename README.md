# radio-imaging

This repository contains the implementation and results for the following papers:
1. Wang, S., Prunet, S., Mignot, S. and Ferrari, A., 2024. Multi-step reconstruction of radio-interferometric images. Astronomy & Astrophysics, 692, p.A61
2. Wang, S., Mignot, S., Prunet, S., Di Mascolo, L., Spinelli, M. and Ferrari, A., 2025. A Decentralized Framework for Radio-Interferometric Image Reconstruction. (tbd)

The instructions below are for the paper "Multi-step reconstruction of radio-interferometric images". Instructions for running the code for "A Decentralized Framework for Radio-Interferometric Image Reconstruction" can be found under the subdirectory "interleaved_parallelization_prototype". 

The code is implemented in Julia, and integrated in rascil. These are not stored directly on this repository but rather in a variety of others, specifically:
1. Julia code: https://github.com/andferrari/DeconvMultiStep.jl
2. RASCIL integration: https://gitlab.com/prunet1/ska-sdp-func-python and https://gitlab.com/igorawratu/rascil-main

## Instructions for running system
The entire system is contained within a docker image, for which the Dockerfile is within this repository. To run: 
1. Clone repository
2. Install docker
3. Build docker image using the command > docker build -t <image_name> .
4. Run the image using > docker run -p 8888:8888 -v <host_directory>:<container_directory> -it <image_name>
5. The docker container will start. You can then run RASCIL from its home directory, as all the required julia files will be copied there.
6. Example commands for execution can be found in the README.rascil file. Simulated datasets can be found in the various data directories. 

## Splitting measurement sets
Measurement sets can be split using casa. Specifically, one can use casatools.ms, for which documentation can be found: https://casadocs.readthedocs.io/en/stable/api/tt/casatools.ms.html#casatools.ms
Update: The dirty.ipynb notebook has been updated to also be able to split the datasets

## Results
The results for the paper are presented in a variety of jupyter notebooks. These are:
- dirty.ipynb for presenting the original datasets and their split
- lambda.ipynb for preliminary experiments on $\lambda$
- noise.ipynb for preliminary experiments and estimation of $\sigma^2$ and $\eta^2$
- split.ipynb for results regarding partition configurations, and against the single-step LASSO reconstruction

The raw data for our results will be made available in the near future