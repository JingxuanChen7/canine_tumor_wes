#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=wes_mapping
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=60G
#SBATCH --time=100:00:00
#SBATCH --mail-user=jc33471@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/scratch/jc33471/canine_tumor/test.out

# install conda env
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/wes_env.yml --name wes_env
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/annovar_env.yml --name annovar_env
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/mutect2_env.yml --name mutect2_env
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/strelka_env.yml --name strelka_env

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate wes_env

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor"
mkdir -p ${run_dir}
cd ${run_dir}
mkdir -p logs/ config/ data/ out/

config="/home/jc33471/canine_tumor_wes/scripts/per_case/config.json"
# config=$run_dir/config/${prefix}.json

snakemake \
    --dry-run \
    --cores ${SLURM_NTASKS} \
    --use-conda \
    --configfile ${config} \
    --snakefile "${project_dir}/scripts/per_case/Snakefile"

/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf/annotate_variation.pl
/work/szlab/kh31516_Lab_Share_script/Add_GeneName_N_Signature.py

python /work/szlab/kh31516_Lab_Share_script/Filter_MutectStat_5steps.py             /scratch/jc33471/canine_tumor/results/Mutect/PRJNA552905/004/004_PASS.stat             /scratch/jc33471/canine_tumor/results/Mutect/PRJNA552905/004/004_rg_added_sorted_dedupped_removed.realigned.MuTect.vcf-PASS             /scratch/jc33471/canine_tumor/results/Mutect/PRJNA552905/004/004_vaf_before.txt             /scratch/jc33471/canine_tumor/results/Mutect/PRJNA552905/004/004_vaf_after.txt             /scratch/jc33471/canine_tumor/results/Mutect/PRJNA552905/004/004_whyout.txt             004