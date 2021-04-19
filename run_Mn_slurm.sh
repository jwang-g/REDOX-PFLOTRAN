#!/bin/bash

#SBATCH -A ccsi
#SBATCH -p batch 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 6:00:00
#SBATCH -o Mn_output/Mn_out_%j.txt
#SBATCH -e Mn_output/Mn_error_%j.txt

source ~/.bashrc

module load anaconda3
conda activate myanaconda3

python Manganese_profile.py