#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --partition=thin
#SBATCH --cpus-per-task=4
#SBATCH --cpu-freq=highm1
#SBATCH --partition=thin
#SBATCH --time=15:00:00

source ${HOME}/load_modules.sh
source ${HOME}/python-venv/bin/activate

echo "$(date) | Starting runs"
python -u check/check_all.py $@ --ncores=128
echo "$(date) | Finished"
