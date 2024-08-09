#!/usr/bin/env bash

# Install the required R packages
Rscript -e "remotes::install_version(package = 'PRECAST', version = '1.6.3', repos = 'https://cran.uni-muenster.de/')"
