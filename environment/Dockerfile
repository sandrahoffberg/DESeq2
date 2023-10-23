# hash:sha256:b288481ab70ed89e052ab03729a55fa16245f078a4007622ac92735210eae4fe
FROM registry.apps.codeocean.com/codeocean/py-r:python3.10.12-R4.2.3-IRkernel-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys \
        0xAD2323F17326AE31401037733E05EBFF05441C52 \
    && apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && add-apt-repository -y 'deb http://deb.codeocean.com/shiny-server-trusty/ ubuntu main' \
    && apt-get purge -y --autoremove software-properties-common

RUN mamba install -y --channel=conda-forge --channel=r \
        r-base==4.3.1 \
        r-essentials==4.3 \
        r-htmltools==0.5.6 \
        r-shiny==1.7.5 \
    && mamba clean -ya

RUN Rscript -e 'remotes::install_version("BiocManager", "1.30.22")' \
    && Rscript -e 'remotes::install_version("svglite", "2.1.1")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "BiocGenerics", \
        "BiocVersion", \
        "DESeq2", \
        "rhdf5", \
        "tximport" \
    ))' # Original versions: 0.46.0 3.17.1 1.40.2 2.44.0 1.28.0
