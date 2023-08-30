#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=wes_mapping
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=60G
#SBATCH --time=100:00:00
#SBATCH --mail-user=jc33471@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/scratch/jc33471/canine_tumor/CMT-100.out

# install conda env
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/wes_env.yml --name wes_env
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/annovar_env.yml --name annovar_env
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/mutect2_env.yml --name mutect2_env
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/strelka_env.yml --name strelka_env
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/java17.yml --name java17
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/data_collection.yml --name data_collection

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate wes_env

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor"
mkdir -p ${run_dir}
cd ${run_dir}
mkdir -p logs/ config/ data/ out/

config="/home/jc33471/canine_tumor_wes/scripts/per_case/config.json"
# config=$run_dir/config/test.json

python ${project_dir}/scripts/per_case/make_snakemake_config.py \
    --project_dir ${project_dir} \
    --out ${run_dir}/config/test.json \
    --outdir ${run_dir} \
    --Bioproject "PRJNA489159" \
    --Normal_Run "SRR7780976" \
    --Tumor_Run "SRR7780979" \
    --CaseName "CMT-100" \
    --threads 8 \
    --memory "60G"

snakemake \
    --dry-run \
    --cores ${SLURM_NTASKS} \
    --use-conda \
    --configfile ${config} \
    --snakefile "${project_dir}/scripts/per_case/Snakefile"

# generate dag image
# snakemake \
#     --dag \
#     --cores ${SLURM_NTASKS} \
#     --use-conda \
#     --configfile ${config} \
#     --snakefile "${project_dir}/scripts/per_case/Snakefile" | dot -Tpdf > dag.pdf

