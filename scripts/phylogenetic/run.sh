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
#SBATCH --output=/scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/sub.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate wes_env

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor/phylogenetics"
breed_dir="/home/${USER}/breed_prediction"
mkdir -p ${run_dir}
cd ${run_dir}

mkdir -p vcf/ merge_vcf/

# config="/home/jc33471/canine_tumor_wes/scripts/phylogenetic/config.json"
config=${run_dir}/config.json

python ${project_dir}/scripts/phylogenetic/make_snakemake_config.py \
    --project_dir ${project_dir} \
    --breed_dir ${breed_dir} \
    --out ${config} \
    --outdir ${run_dir}/merge_vcf \
    --metadata "/scratch/jc33471/canine_tumor/breed_prediction/this_wunpaired_meta.csv" \
    --vcffilelist "/scratch/jc33471/canine_tumor/breed_prediction/vcf_file_list.txt" \
    --breedSpecific "/scratch/jc33471/canine_tumor/breed_prediction/output_exclude_WGS/all_breed_specific_variants.txt" \
    --somaticMutation ${project_dir}"/metadata/Pass_QC_Final_Total_withGene_Burair_Filtering4_VAF_Mutect_orientBiasModified_04_02.txt" \
    --threads 8 \
    --memory "60G"

snakemake \
    -pn \
    --latency-wait 60 \
    --cores ${SLURM_NTASKS} \
    --rerun-incomplete \
    --use-conda \
    --rerun-triggers mtime \
    --configfile ${config} \
    --snakefile "${project_dir}/scripts/phylogenetic/Snakefile"
