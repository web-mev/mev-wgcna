# mev-wgcna

This repository contains a WebMeV-compatible tool for performing weighted correlation network analysis using Bioconductor's WGCNA package (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)

The outputs of this analysis include
- a file providing the gene modules resulting from the clustering
- a QC file which summarizes which "network threshold" was chosen. 

Note that WGCNA relies on a user-defined network connectivity parameter which dictates the clustering on gene co-expression and the resulting network topology. By default (if unspecified), we use a heuristic algorithm to estimate this parameter based on the location of the "elbow" of a plot (see below). Roughtly, the y-axis value of the plot is related to a fit to a "scale-free" network topology and the authors recommend that the $\beta$ value is chosen to be such that we select the elbow. Since this heuristic may not work well in all cases, the second output file provides data to create the plot below which can help users identify if the heuristic algorithm chose an appropriate threshold for their data. In most cases the heuristic selection works pretty well. However, if the selected point is not close to the elbow, the plot can be used to manually pick a threshold and re-run the analysis.

![Alt text](./elbow.svg)
<img src="./elbow.svg">

---

### To run external of WebMeV:

Either:
- build the Docker image using the contents of the `docker/` folder (e.g. `docker build -t myuser/wgcna:v1 .`) 
- pull the docker image from the GitHub container repository (see https://github.com/web-mev/mev-wgcna/pkgs/container/mev-wgcna)

To run, change into the directory where you have your expression matrix. Then:
```
Rscript /usr/local/bin/wgcna.R \
    -f /work/<expression matrix> \
    -b (optional) <network connectivity parameter>
```
Note that we mount the current directory at `/work`, hence the expression matrix is relative to `/work`.

The call to the script assumes the following:
- The input file of expression counts is tab-delimited format and contains *normalized* expression data.