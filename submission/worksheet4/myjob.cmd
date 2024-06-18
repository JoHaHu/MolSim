#!/bin/bash
#SBATCH -J molsim_testrun_group_d
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --time=01:00:00

./myprog.exe