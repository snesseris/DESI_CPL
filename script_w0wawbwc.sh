#!/bin/sh
#SBATCH -J SPARTA
#SBATCH --partition=batch4
#SBATCH --get-user-env
#SBATCH -o ./job.%j.out
#SBATCH -e ./job.%j.err
#SBATCH -n 1
#SBATCH -t 2-01:59:00

wolframscript ./script_MCMC_w0wawbwc.wls