# hash:sha256:f5046b9250bc54513b87effdb4a8bc57d8e9f09e5e57dd2b1de1df3e78a7a73d
FROM registry.apps.codeocean.com/codeocean/r-studio:2022.07.0-548-r4.2.1-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        pandoc=1.19.2.4~dfsg-1build4 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("BiocManager", "1.30.21")' \
    && Rscript -e 'remotes::install_version("markdown", "1.7")' \
    && Rscript -e 'remotes::install_version("svglite", "2.1.1")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "BiocGenerics", \
        "BiocVersion", \
        "DESeq2" \
    ))' # Original versions: 0.44.0 3.16.0 1.38.3
