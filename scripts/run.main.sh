#!/bin/sh
###################################################################################################
#SBATCH --job-name SRC-algo-MSK
#SBATCH --mail-type=NONE
#SBATCH --partition=F16
#SBATCH --job-name sc-sim
#SBATCH --array=1-3
#SBATCH --mem=5G #memory allocation per job
#SBATCH -c 1 #cpu allocation per job
#SBATCH --partition=F32
#SBATCH --nodes=2
#SBATCH --nodelist=F32-[9-10]
###################################################################################################
source ~/.bashrc && conda activate r_env

NCELL=(50 100 200)
NCLONE=(3 5 10)
NSNV=(50 100 200)
noise=(0.2 0.1 0.05)

for i in "${NCELL[${SLURM_ARRAY_TASK_ID}]}"; do
# for i in "${NCELL[@]}"; do
  for j in "${NCLONE[@]}"; do
    for a in "${NSNV[@]}"; do
      for b in "${noise[@]}"; do
        echo "Running sim_K${j}_M${a}_C${i}_N${b}"
        /usr/bin/time -f 'real(s)\tusr(s)\tsys(s)\tmem(kB)\tcpu(%)\n%e\t%U\t%S\t%M\t%P\n' \
        Rscript /hot/user/hwinata/BIOMATH211_Cancer-Phylogeny/scripts/main.R \
        --nclone $j \
        --ncell $i \
        --nsnv $a \
        --noise $b \
        &> "/hot/user/hwinata/BIOMATH211_Cancer-Phylogeny/data/sim-output/log/sim_K${j}_M${a}_C${i}_N${b}.timemem"
        # touch "/hot/user/hwinata/BIOMATH211_Cancer-Phylogeny/data/sim_data/log/sim_K${i}_M${a}_C${j}_N${b}.timemem"
      done
    done
  done
done