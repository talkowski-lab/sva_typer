FROM rust:1.86

# WORKDIR /usr/src
#
# RUN git clone https://github.com/talkowski-lab/sva_typer.git 
#
# RUN cd sva_typer && cargo install --path .

# install conda
RUN wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
    bash Miniforge3.sh -b -p "/opt/conda"

RUN . "/opt/conda/etc/profile.d/conda.sh"
# For mamba support also run the following command
RUN . "/opt/conda/etc/profile.d/mamba.sh"

ENV PATH="/opt/conda/bin:${PATH}"
ARG PATH="/opt/conda/bin:${PATH}"

RUN conda update -y -n base conda && \
    conda install -y -c conda-forge conda-pack libmamba && \
    conda config --set solver libmamba && \
    conda clean --all --yes

WORKDIR /usr/src/sva_typer

COPY . .

RUN cargo install --path .

# copy other resources
COPY ./environment.yml /
# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH="/opt/conda/envs/sva_typer/bin:${PATH}"
ENV SVA_HMM_PATH="/usr/src/sva_typer/ref"

# copy python scripts
COPY scripts/* /scripts/

# activate conda environment
RUN echo "source activate seq_extract" > ~/.bashrc

# CMD ["sva_typer"]
