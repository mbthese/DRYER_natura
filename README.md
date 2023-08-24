# About

This is a snakemake workflow based on the obitools suite of programs, that analyzes DNA metabarcoding data. Sequence analysis is performed with the obitools (Boyer et al. 2016) and sumaclust (Mercier et al. 2013) through a Snakemake pipeline (Mölder et al. 2021).

It is also where all analyses scripts are stored for the article.

# Directories and files structure

The repository contains five folders:
- `config/`: contains the configuration file of the Snakemake workflow (`config.yaml`).
- `log/`: where log files of each rule are written.
- `resources/`: where raw data is stored.
- `results/`: where all output files are written.
- `script/`: where all scripts are stored for the article analyses.
- `workflow/`: contains the Snakemake workflow (`Snakefile`), the configuration file of the submission parameters on the cluster (`cluster.yaml`), the script to submit the workflow on the cluster (`sub_smk.sh`). 


# Acknowledgements

This github was initally forked from the following repository for high throughput sequencing analyses workflow:
Anne-Sophie Benoiston. (2022). AnneSoBen/obitools_workflow: v1.0.2. GitHub. https://doi.org/10.5281/zenodo.6676577.

