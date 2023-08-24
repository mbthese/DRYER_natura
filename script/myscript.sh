#!/bin/bash
#SBATCH -J lauchRscript
#SBATCH -o output.out

#SBATCH --mail-type=BEGIN,END,FAIL

module purge

module load system/R-3.5.1

Rscript NDV_script.R