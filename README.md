
# pMPBA

This repository contains pipeline for data processing and Jupyter notebooks scripts used for data analysis.

## Running pipeline (snakemake)
download all fastq files from GEO accession number: GSE272404 and put them in a folder termed 'data', and then run pipeline:
`./workflow/Snakefile.sh` 

This will result in csv files for each experiment in folder 'res_files'. These processed files are already provided.

## data analysis
scripts for data analysis can be found in folder 'analysis', containing jupyter notebook scripts and additional data used in analysis. 

The script `./analysis/createTables.ipynb` will generate all supplemental tables which are used in all other scripts.


