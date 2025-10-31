#!/bin/bash


rm -r environment/radioimaging
mkdir environment

#use the below if using conda instead

# conda create -p environment/radioimaging python=3.10.12
# conda activate environment/radioimaging

virtualenv environment/radioimaging
source environment/radioimaging/bin/activate

rm -r external_dependencies
mkdir external_dependencies

cd external_dependencies
git clone https://gitlab.com/ska-telescope/external/rascil-main
cd rascil-main
git checkout tags/1.1.0

cp ../../hotfixes/rascil-main/Makefile .
cp ../../hotfixes/rascil-main/requirements.in .
cp ../../hotfixes/rascil-main/simulation_helpers.py rascil/processing_components/simulation/
cp ../../hotfixes/rascil-main/msv2.py rascil/processing_components/visibility/

make requirements
make install_requirements
cd ../../

pip install -r requirements.txt
pip install -e .

pip install notebook

#uncomment if you wish to set a custom package install path. Often necessary for clusters since the home directory is somewhat limited in space

#export JULIA_DEPOT_PATH=path/to/julia/packages

julia -e 'import Pkg; Pkg.update()' && \
julia -e 'import Pkg; Pkg.add("FITSIO")' && \
julia -e 'import Pkg; Pkg.add(url="https://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@github.com/andferrari/IUWT.jl")' && \
julia -e 'import Pkg; Pkg.add(url="https://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@github.com/andferrari/DeconvMultiStep.jl")' 