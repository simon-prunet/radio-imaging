FROM "artefact.skao.int/rascil-full:1.0.0"

USER rascil

# install julia

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

RUN python3 -m pip install julia 
