FROM continuumio/miniconda3:4.9.2

# Install hhsuite, r-base=4.2.2, and aspring==1.0.0
RUN conda install -y -c conda-forge -c bioconda hhsuite && \
    conda install -y -c conda-forge r-base=4.2.2 && \
    pip install aspring==1.0.0

# Install renv package for R and restore renv environment for aspring R script
RUN Rscript -e 'options(repos = c(CRAN = "https://cran.rstudio.com/")); install.packages("renv")' && \
    Rscript -e 'renv::restore(project="/opt/conda/lib/python3.8/site-packages/aspring/R_script")'

ENV HHSUITE_SCRIPTS /opt/conda/scripts

VOLUME /data

WORKDIR /data

CMD ["aspring"]
