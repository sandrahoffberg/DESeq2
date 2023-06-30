# hash:sha256:5442e25f5ae104bd619d58ab254e2b07ca9c4ea58872decea2f0fc4909ebd8e7
FROM registry.apps.codeocean.com/codeocean/r-studio:2022.07.0-548-r4.2.1-ubuntu18.04

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

RUN Rscript -e 'remotes::install_version("BiocManager", "1.30.21")' \
    && Rscript -e 'remotes::install_version("svglite", "2.1.1")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "BiocGenerics", \
        "BiocVersion", \
        "DESeq2", \
        "tximport" \
    ))' # Original versions: 0.44.0 3.16.0 1.38.3 1.26.1
