#! /bin/bash
#SBATCH --job-name=powerNormal  # Job name
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=petr.smirnov@mail.utoronto.ca     # Where to send mail
#SBATCH --nodes=1
#SBATCH --array=0#-1000#0-816
#SBATCH --ntasks-per-node=40
#SBATCH --time=12:00:00
#SBATCH --output=/cluster/projects/bhklab/tmp/psmirnov_scratch/runRCIPowerNormal_%A_%a.log   # Standard output and error log
#SBATCH --mem=60GB
#SBATCH -p superhimem


source ~/.bashrc


module load R

cd ~/Github/wci_manuscript

Rscript $HOME/Github/wci_manuscript/SimulatePowerNormal.R 

