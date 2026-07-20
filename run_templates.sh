#!/bin/bash
#SBATCH -A b1042             ## account (b1012 for buyin, genomics = b1042)
#SBATCH -p genomics          ## "-p" instead of "-q"
#SBATCH -J run_seurat_pipeline        ## job name
#SBATCH --mail-type=FAIL,TIME_LIMIT_90
#SBATCH --mail-user=your_email@northwestern.edu
#SBATCH -o "%x.o%j"
#SBATCH -N 1                 ## number of nodes
#SBATCH -n 1                 ## number of cores
#SBATCH -t 48:00:00          ## walltime
#SBATCH --mem=80G
export MC_CORES=${SLURM_NTASKS}
module purge all
module load R/4.4.0
module load geos/3.8.1
module load hdf5/1.8.19-serial
module load pandoc/2.2.1

# Pass one or more .Rmd files. They are rendered in order against one config
# (config.yaml by default; override with --export=ALL,CONFIG_FILE=other.yaml).
#
# Examples:
#   sbatch run_templates.sh scRNA_part1_QC.Rmd
#   sbatch run_templates.sh scRNA_part2a_integration.Rmd scRNA_part2b_clustering.Rmd
#   sbatch run_templates.sh scRNA_part2b_clustering.Rmd        # re-cluster only
#   sbatch run_templates.sh scRNA_part3_summary.Rmd
#   sbatch --export=ALL,CONFIG_FILE=config_test.yaml run_templates.sh scRNA_part1_QC.Rmd # for alt config file


echo -n "Running pipeline file(s): $* at "
date

Rscript execute_pipeline-Ex.R "$@"

echo -n "Done at "
date
