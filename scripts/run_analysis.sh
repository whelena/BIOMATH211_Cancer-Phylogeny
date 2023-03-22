#!/bin/sh
#SBATCH --time=23:59:00
#SBATCH --job-name SRC-algo-MSK
#SBATCH --mail-type=NONE
#SBATCH --partition=F16
###################################################################################################
source ~/.bashrc && conda activate r_env

SIMS=$(basename -a /hot/user/hwinata/BIOMATH211_Cancer-Phylogeny/data/input/simulation/*src_input* | sed 's\-src_input.tsv\\')
LOG=/hot/user/hwinata/BIOMATH211_Cancer-Phylogeny/data
touch $LOG/simulation.log

for simid in $SIMS; do
    echo "Running $simid" >> $LOG/simulation.log
    Rscript /hot/user/hwinata/BIOMATH211_Cancer-Phylogeny/scripts/simulate.single.cell.R \
    --name $simid >> $LOG/simulation.log
    done