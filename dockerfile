FROM rocker/geospatial
LABEL maintainer="Jesse Piburn <piburnjo@ornl.gov>"

# Install 
RUN apt-get -y --no-install-recommends install \
    ed \
    clang  \
    ccache \
    && install2.r --error --deps TRUE \
        StanHeaders \
        rstan \
        tidybayes \
        data.table \
        foreach \
        doParallel \
        furrr \
    && R -e "remotes::install_github('nset-ornl/covidmodeldata');"
    

COPY . /home/rstudio/covid-model/
