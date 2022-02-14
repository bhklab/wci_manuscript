# hash:sha256:92d615b573564d897c25f45c5450c5a6e4ffaf7c111e926fd223ed999ae768e1
FROM registry.codeocean.com/codeocean/r-base:4.0.3-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'remotes::install_version("MCMCpack")' \
    && Rscript -e 'remotes::install_version("abind")' \
    && Rscript -e 'remotes::install_version("boot")' \
    && Rscript -e 'remotes::install_version("data.table")' \
    && Rscript -e 'remotes::install_version("doParallel")' \
    && Rscript -e 'remotes::install_version("ggExtra")' \
    && Rscript -e 'remotes::install_version("ggplot2")' \
    && Rscript -e 'remotes::install_version("iterators")' \
    && Rscript -e 'remotes::install_version("mvtnorm")' \
    && Rscript -e 'remotes::install_version("profmem")' \
    && Rscript -e 'remotes::install_version("pwr")' \
    && Rscript -e 'remotes::install_version("reshape2")' \
    && Rscript -e 'remotes::install_version("rmutil")' \
    && Rscript -e 'remotes::install_version("tictoc")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "PharmacoGx" \
    ))' # Original versions: latest
