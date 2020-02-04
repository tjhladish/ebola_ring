#!/bin/bash
#SBATCH --job-name=ebola_ring
#SBATCH --output=./auto_output/ebola_ring_%A_%a.out
#SBATCH --error=./auto_output/ebola_ring_%A_%a.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --qos=epi
#SBATCH --chdir=/home/tjhladish/work/ebola
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=24:00:00
#SBATCH --array=0-999
#SBATCH --partition=hpg2-compute

module load gcc/7.3.0 gsl

for i in `seq 1 20`;
do
    ./ebola_net abc_ebola_rings-test.json --simulate
    #./ebola_net abc_ebola_rings.json --simulate -n 20
done
