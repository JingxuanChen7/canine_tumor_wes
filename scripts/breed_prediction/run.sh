#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=breed
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=60G
#SBATCH --time=100:00:00
#SBATCH --mail-user=jc33471@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/scratch/jc33471/canine_tumor_test/breed.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate breed_prediction

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor_test"
mkdir -p ${run_dir}
cd ${run_dir}

Rscript --vanilla ${project_dir}/scripts/breed_prediction/germline_VAF_reset_low_coverage.R

Rscript --vanilla ${project_dir}/scripts/breed_prediction/breed_specific_variants.R
