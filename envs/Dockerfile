FROM continuumio/miniconda:4.7.12

WORKDIR /home/snakemake/

COPY global_env.yaml ./

RUN conda env create -f global_env.yaml

RUN echo "conda activate BSanalysis" > ~/.bashrc
ENV PATH /opt/conda/envs/BSanalysis/bin:$PATH

ENTRYPOINT ["snakemake"]
