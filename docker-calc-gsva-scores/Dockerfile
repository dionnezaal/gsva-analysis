FROM r-base

## install necessary Rpackages for GSVA analysis
# Install packages to be able to run GSVA analysis in parallel
RUN apt-get update && apt-get install -y \
	libcurl4-openssl-dev \
	libxml2-dev
RUN R -e "install.packages(c('snow', 'DBI', 'RSQLite', 'RCurl', 'xtable', 'openssl', 'Rmpi', 'covr', 'rlecuyer'), repos='http://cran-mirror.cs.uu.nl/', dependencies=TRUE)"
RUN R -e "install.packages(c('qusage', 'GSEABase', 'GSVA'), repos='http://bioconductor.org/packages/3.5/bioc')"

COPY . /usr/local/src/scripts/
WORKDIR /usr/local/src/scripts/

RUN ["chmod", "+x", "calculate_GSVA_scores.R"]
