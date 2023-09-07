#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=master_job
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=60G
#SBATCH --time=500:00:00
#SBATCH --mail-user=jc33471@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/scratch/jc33471/canine_tumor/master_job.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate wes_env

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor"
mkdir -p ${run_dir}
cd ${run_dir}
mkdir -p logs/ config/ data/ out/ results

config="/home/jc33471/canine_tumor_wes/scripts/sum_cases/config.json"
# config=$run_dir/config/test.json

python ${project_dir}/scripts/sum_cases/make_snakemake_config.py \
    --project_dir ${project_dir} \
    --out ${run_dir}/config/test.json \
    --outdir ${run_dir} \
    --metadata ${project_dir}/metadata/data_collection_new.csv\
    --threads 8 \
    --memory "60G"

snakemake \
    --dry-run \
    --jobs 200 \
    --use-conda \
    --latency-wait 60 \
    --keep-going \
    --restart-times 0 \
    --configfile ${config} \
    --snakefile "${project_dir}/scripts/sum_cases/Snakefile" \
    --cluster '
        sbatch \
            --partition=batch \
            --job-name={wildcards.Bioproject}_{wildcards.CaseName}_smk \
            --nodes=1 \
            --tasks-per-node={threads} \
            --mem={resources.mem} \
            --time=72:00:00 \
            --mail-user=jc33471@uga.edu \
            --mail-type=FAIL \
            --output=logs/{wildcards.Bioproject}_{wildcards.CaseName}_log.o \
            --error=logs/{wildcards.Bioproject}_{wildcards.CaseName}_log.e'


