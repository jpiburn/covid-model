FROM rocker/geospatial:3.6.3
LABEL AUTHORS="Piburn & Grant" \
  MAINTAINER="Joshua Grant <grantjn@ornl.gov>"

# Install required packages
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    apt-utils ed libnlopt-dev liblapack-dev libblas-dev libboost-dev \
    libopenblas-dev libeigen3-dev \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/

COPY covidmodeldata.tar.gz /tmp/covidmodeldata.tar.gz

# Install R packages
RUN install2.r --error --deps true StanHeaders rstan loo tidybayes data.table \
      foreach doParallel bayesplot furrr rstanarm rstantools shinystan \
      ggmcmc snakecase janitor \
    && R CMD INSTALL /tmp/covidmodeldata.tar.gz \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds /tmp/covidmodeldata.tar.gz

COPY ./scripts /tmp/scripts

ENTRYPOINT /tmp/scripts/check_nyt_update.py
