#!/bin/sh
###################################################################################################
#SBATCH --job-name SRC-algo-MSK
#SBATCH --mail-type=NONE
#SBATCH --partition=F16
#SBATCH --job-name sc-sim
#SBATCH --array=1-4
#SBATCH --mem=5G #memory allocation per job
#SBATCH -c 1 #cpu allocation per job
#SBATCH --partition=F32
#SBATCH --nodes=2
#SBATCH --nodelist=F32-[9-10]
###################################################################################################
source ~/.bashrc && conda activate r_env

NCELL=(10 50 100 500)
NCLONE=(3 5 10)
NSNV=(50 100 500)
noise=(0.2 0.1 0.05)

# function run_main {
#   /usr/bin/time -f 'real(s)\tusr(s)\tsys(s)\tmem(kB)\tcpu(%)\n%e\t%U\t%S\t%M\t%P\n' \
#     Rscript /hot/user/hwinata/BIOMATH211_Cancer-Phylogeny/scripts/main.R \
#     --nclone $i \
#     --ncell $j \
#     --nsnv $a \
#     --noise $b \
#     &> "/hot/user/hwinata/BIOMATH211_Cancer-Phylogeny/data/log/sim_K${i}_M${a}_C${j}_N${b}.timemem"
#     }

# for i in "${NCELL[${SLURM_ARRAY_TASK_ID}]}"; do
for i in "${NCELL[@]}"; do
  for j in "${NCLONE[@]}"; do
    for a in "${NSNV[@]}"; do
      for b in "${noise[@]}"; do
        echo "Running sim_K${j}_M${a}_C${i}_N${b}"
        attempt=1
        while [ $attempt -le 3 ]
          do
            /usr/bin/time -f 'real(s)\tusr(s)\tsys(s)\tmem(kB)\tcpu(%)\n%e\t%U\t%S\t%M\t%P\n' \
            Rscript /hot/user/hwinata/BIOMATH211_Cancer-Phylogeny/scripts/main.R \
            --nclone $j \
            --ncell $i \
            --nsnv $a \
            --noise $b \
            &> "/hot/user/hwinata/BIOMATH211_Cancer-Phylogeny/data/log/sim_K${j}_M${a}_C${i}_N${b}.timemem"

            if [ $? -eq 0 ]
            then
              echo "Script completed successfully"
            else
              echo "Script failed, attempt $attempt of 3"
              attempt=$((attempt+1))
            fi
          done
        # touch "/hot/user/hwinata/BIOMATH211_Cancer-Phylogeny/data/log/sim_K${i}_M${a}_C${j}_N${b}.timemem"
      done
    done
  done
done