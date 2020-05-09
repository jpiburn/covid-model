FROM rocker/geospatial
LABEL maintainer="Jesse Piburn <piburnjo@ornl.gov>"

# ------------------------------
# Install rstan and friends
# ------------------------------
# Docker Hub (and Docker in general) chokes on memory issues when compiling
# with gcc, so copy custom CXX settings to /root/.R/Makevars and use ccache and
# clang++ instead

# Make ~/.R
RUN mkdir -p $HOME/.R

# $HOME doesn't exist in the COPY shell, so be explicit
COPY R/Makevars /root/.R/Makevars

# Install 
RUN apt-get -y --no-install-recommends install \
    ed \
    clang  \
    ccache \
    && install2.r --error --deps TRUE \
        StanHeaders \
        rstan \
    && R -e "remotes::install_github('nset-ornl/covidmodeldata');"