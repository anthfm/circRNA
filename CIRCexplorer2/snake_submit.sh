#!/bin/bash

#SBATCH --job-name=snake_circ
#SBATCH --mem=3000M
#SBATCH -c 1
#SBATCH -t 3:00:00

source /cluster/home/amammoli/.bashrc
module load python3




output_dir="/cluster/projects/bhklab/Data/ncRNA_detect/circRNA/CIRCexplorer2/test"  snakemake -s /cluster/projects/bhklab/Data/ncRNA_detect/circRNA/CIRCexplorer2/test/Snakefile --latency-wait 100 -j 200 --cluster 'sbatch -t {params.runtime} --mem={resources.mem_mb} -c {threads}'
