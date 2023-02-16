# FROM "artefact.skao.int/rascil-full:1.0.0"
FROM "python:3.10"

# install julia

WORKDIR /
USER root
  
ENV JULIA_VERSION=1.8.5

RUN mkdir /opt/julia-${JULIA_VERSION} && \
    cd /tmp && \
    wget -q https://julialang-s3.julialang.org/bin/linux/x64/`echo ${JULIA_VERSION} | cut -d. -f 1,2`/julia-${JULIA_VERSION}-linux-x86_64.tar.gz && \
    tar xzf julia-${JULIA_VERSION}-linux-x86_64.tar.gz -C /opt/julia-${JULIA_VERSION} --strip-components=1 && \
    rm /tmp/julia-${JULIA_VERSION}-linux-x86_64.tar.gz

RUN ln -fs /opt/julia-*/bin/julia /usr/local/bin/julia

# Add julia packages
 
RUN julia -e 'import Pkg; Pkg.update()' && \
    julia -e 'import Pkg; Pkg.add("PyCall")' && \
    julia -e 'ENV["PYTHON"] = "/usr/local/bin/python"; import Pkg; Pkg.build("PyCall")'  && \
    julia -e 'import Pkg; Pkg.add(url="https://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@github.com/andferrari/DeconvMultiStep.jl")'

# PyJulia

RUN pip install --upgrade pip && python3 -m pip install julia

# Add tailored version of ska-sdp-func-python from gitlab fork. Build and pip install

# Install poetry

RUN curl -sSL https://install.python-poetry.org | python3 -
ENV PATH=$PATH:/root/.local/bin

# Install compilers, cmake, etc.
RUN apt-get update && \
    apt-get install -y cmake && \
    apt-get install -y g++ 

# Add tailored version of ska-sdp-datamodels from gitlab fork, for Meerkat-only config. Build and pip install

RUN cd /tmp && \
    git clone https://oauth2:glpat-sRXuumo5FrvWzznszuLw@gitlab.com/prunet1/ska-sdp-datamodels.git && \
    cd ska-sdp-datamodels && \ 
    poetry build && \
    pip install dist/ska_sdp_datamodels-0.1.3-py3-none-any.whl 

# Add tailored version of ska-sdp-func-python from gitlab fork. Build and pip install

RUN cd /tmp && \
    git clone https://oauth2:glpat-sRXuumo5FrvWzznszuLw@gitlab.com/prunet1/ska-sdp-func-python.git && \
    cd ska-sdp-func-python && \
    poetry build && \
    pip install dist/ska_sdp_func_python-0.1.4-py3-none-any.whl

# Rascil, built from sources
RUN cd / && \
    git clone https://gitlab.com/ska-telescope/external/rascil-main.git && \
    cd /rascil-main && \
    pip install pip-tools && pip-compile --resolver=backtracking requirements.in && \
    pip install -r requirements.txt

ENV RASCIL=/rascil-main
ENV PYTHONPATH=$RASCIL:$PYTHONPATH

#RUN cd $RASCIL && apt-get install git-lfs && \
#    git-lfs install && \
#    git-lfs pull


