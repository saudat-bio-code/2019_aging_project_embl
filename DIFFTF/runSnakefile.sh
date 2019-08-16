#!/bin/bash

# Cluster-specific settings, customized for the EMBL SLURM cluster.
# Important: Only needed if conda is NOT used; otherwise, conda provides all necessary tools so no moduels are needed except for Snakemake
# Unload all modules and load the necessary ones for the pipeline

# module purge
# module load GCCcore ncurses BEDTools SAMtools R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.0 Autoconf FastQC Trimmomatic snakemake MACS2 deepTools

########################
# PATHS AND PARAMETERS #
########################

# All parameters and paths are defined here:
. "./params.sh"




########################
# RUN AUTOMATED SCRIPT #
########################
. "/g/scb/zaugg/zaugg_shared/scripts/Christian/src/Snakemake/runSnakemakeWrapper.sh"
