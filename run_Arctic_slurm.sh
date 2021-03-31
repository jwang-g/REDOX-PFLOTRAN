#!/bin/bash

#SBATCH -A ccsi
#SBATCH -p batch 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 2:00:00
#SBATCH -o Arctic_Fe_output/Arctic_out_%j.txt
#SBATCH -e Arctic_Fe_output/Arctic_error_%j.txt

source ~/.bashrc

module load anaconda3
conda activate myanaconda3

python Arctic_Fe_CH4.py