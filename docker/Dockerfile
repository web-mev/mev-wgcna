FROM r-base:3.6.1

RUN apt-get update && \
  apt-get install -y libxml2-dev \
  libcurl4-openssl-dev \
  libssl-dev

# Install all software under /opt/software
RUN mkdir -p /opt/software
ADD wgcna.R /opt/software/

# Install directly from github maintained by the developer (instead of Bioc repo)
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("WGCNA")'
RUN R -e 'install.packages("optparse")'

ENTRYPOINT ["/bin/bash"]
