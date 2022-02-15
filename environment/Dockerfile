# hash:sha256:6e270042f603a99f370c1b52d44e52086d93688972fb0d331a1a458781aa2430
FROM registry.codeocean.com/codeocean/r-base:4.0.3-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys \
        0xAD2323F17326AE31401037733E05EBFF05441C52 \
    && apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && add-apt-repository -y 'deb http://deb.codeocean.com/rstudio-server-bionic/ ubuntu main' \
    && apt-get purge -y --autoremove software-properties-common \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        pandoc=1.19.2.4~dfsg-1build4 \
        rstudio-server=1.2.5033 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("MCMCpack", "1.6-0")' \
    && Rscript -e 'remotes::install_version("abind", "1.4-5")' \
    && Rscript -e 'remotes::install_version("boot", "1.3-28")' \
    && Rscript -e 'remotes::install_version("data.table", "1.14.2")' \
    && Rscript -e 'remotes::install_version("doParallel", "1.0.17")' \
    && Rscript -e 'remotes::install_version("ggExtra", "0.9")' \
    && Rscript -e 'remotes::install_version("ggplot2", "3.3.5")' \
    && Rscript -e 'remotes::install_version("iterators", "1.0.14")' \
    && Rscript -e 'remotes::install_version("mvtnorm", "1.1-3")' \
    && Rscript -e 'remotes::install_version("profmem", "0.6.0")' \
    && Rscript -e 'remotes::install_version("pwr", "1.3-0")' \
    && Rscript -e 'remotes::install_version("reshape2", "1.4.4")' \
    && Rscript -e 'remotes::install_version("rmutil", "1.1.8")' \
    && Rscript -e 'remotes::install_version("tictoc", "1.0.1")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "PharmacoGx" \
    ))' # Original versions: 2.2.4
