#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=phylogenetics
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --tasks-per-node=8
#SBATCH --mem=60G
#SBATCH --time=500:00:00
#SBATCH --mail-user=jc33471@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/merge_vcf.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate wes_env

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor/phylogenetics"
mkdir -p ${run_dir}
cd ${run_dir}

mkdir -p vcf/ merge_vcf/

config="/home/jc33471/canine_tumor_wes/scripts/phylogenetic/config.json"

snakemake \
    -p \
    --cores ${SLURM_NTASKS} \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --use-conda \
    --configfile ${config} \
    --snakefile "${project_dir}/scripts/phylogenetic/Snakefile"
