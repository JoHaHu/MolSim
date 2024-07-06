#!/bin/bash
#SBATCH -J molsim_d
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./MolSim/build
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --time=01:00:00
module load slurm_setup
module load xerces-c

./MolSim -i ../input/eingabe-collision-ws2.xml