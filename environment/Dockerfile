# hash:sha256:08da8081b6fb70b8069fd0330f2abda44c11767c62b29f189169f9d461b36644
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
    && Rscript -e 'remotes::install_version("markdown", "1.7")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "BiocGenerics", \
        "BiocVersion", \
        "DESeq2" \
    ))' # Original versions: 0.44.0 3.16.0 1.38.3
