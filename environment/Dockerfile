# hash:sha256:c217645501972c55188215d7a5bcba526badec7ba8907891a27e863b28b34ad6
FROM registry.apps.codeocean.com/codeocean/r-studio:2022.07.0-548-r4.2.1-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "BiocGenerics", \
        "DESeq2" \
    ))' # Original versions: 0.44.0 1.38.3
