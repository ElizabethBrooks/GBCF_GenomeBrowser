FROM rocker/r-ver:4.3.0

RUN apt-get update && apt-get install -y \
  libharfbuzz-dev \
  libfribidi-dev \
  libxml2-dev \
  libtiff-dev \
  wget

RUN apt-get install -y --no-install-recommends software-properties-common dirmngr

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN apt-get install -y --no-install-recommends r-base

