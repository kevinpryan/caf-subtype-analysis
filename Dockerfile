FROM nfcore/base:1.14
LABEL authors="Kevin Ryan" \
      description="Docker container for CAF subtype analysis"
WORKDIR ./
COPY environment.yml ./
RUN conda env create -f environment.yml && conda clean -a
ENV PATH /opt/conda/envs/subtype/bin:$PATH
RUN R -e "devtools::install_github('Danko-Lab/BayesPrism/BayesPrism')"
