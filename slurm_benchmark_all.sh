#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --partition=thin
#SBATCH --cpus-per-task=1
#SBATCH --cpu-freq=highm1
#SBATCH --partition=thin
#SBATCH --time=12:00:00

source ${HOME}/load_modules.sh
source ${HOME}/python-venv/bin/activate

echo "Starting runs"
python check/check_all.py
echo "Finished"
