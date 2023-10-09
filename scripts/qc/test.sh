#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=test_qc
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --tasks-per-node=8
#SBATCH --mem=60G
#SBATCH --time=500:00:00
#SBATCH --mail-user=jc33471@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/scratch/jc33471/canine_tumor/test_qc.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate wes_env

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor"
mkdir -p ${run_dir}
cd ${run_dir}
mkdir -p logs/ config/ data/ out/ results/

config="/home/jc33471/canine_tumor_wes/scripts/qc/config.json"
# config=$run_dir/config/test_qc.json

# python ${project_dir}/scripts/qc/make_snakemake_config.py \
#     --project_dir ${project_dir} \
#     --out ${config} \
#     --outdir ${run_dir} \
#     --Bioproject "PRJNA489159" \
#     --Normal_Run "SRR7780976" \
#     --Tumor_Run "SRR7780979" \
#     --CaseName "CMT-100" \
#     --metadata ${project_dir}/metadata/data_collection_old.csv \
#     --readlength ${project_dir}/metadata/data_new_readlength.csv \
#     --threads 8 \
#     --memory "60G"

snakemake \
    --cores ${SLURM_NTASKS} \
    --use-conda \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --configfile ${config} \
    --snakefile ${project_dir}"/scripts/qc/Snakefile"

cd $run_dir
java -Xmx"60G" "/home/jc33471/canine_tumor_wes/scripts/qc/Line_by_line_Total_dict_Get_exon_reads.java" \
    "/scratch/jc33471/canine_tumor/results/QC/PRJNA489159/CMT-100/SRR7780976.sam" \
    "SRR7780976" \
    "/scratch/jc33471/canine_tumor/results/QC/PRJNA489159/CMT-100/SRR7780976_MT_PRJNA489159-CDS_Mapping_summary.txt;" \
    101 \
    "MT" \
    "Normal" 