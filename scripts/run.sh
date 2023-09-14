#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=60G
#SBATCH --time=100:00:00
#SBATCH --mail-user=jc33471@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/scratch/jc33471/canine_tumor/test.out

# install conda env
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/wes_env.yml --name wes_env
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/annovar_env.yml --name annovar_env
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/mutect2_env.yml --name mutect2_env
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/strelka_env.yml --name strelka_env
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/java17.yml --name java17
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/data_collection.yml --name data_collection
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/breed_prediction.yml --name breed_prediction

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate wes_env

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor"
mkdir -p ${run_dir}
cd ${run_dir}
mkdir -p logs/ config/ data/ out/

# config="/home/jc33471/canine_tumor_wes/scripts/per_case/config.json"
config=$run_dir/config/test.json

python ${project_dir}/scripts/per_case/make_snakemake_config.py \
    --project_dir ${project_dir} \
    --out ${config} \
    --outdir ${run_dir} \
    --Bioproject "PRJNA677995" \
    --Normal_Run "SRR13050752" \
    --Tumor_Run "SRR13050748-SRR13050750-SRR13050751" \
    --CaseName "Pt03" \
    --threads 8 \
    --memory "60G"

snakemake \
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

## validate 10 samples
for case in "CMT-100" "CMT-102" "CMT-103" "CMT-105" "CMT-106" "CMT-107" "CMT-109" "CMT-111" "CMT-112" "CMT-114"; do
    wc -l /project/szlab/Kun_Lin/Pan_Cancer/Mammary_Cancer/Germline/${case}/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName
    wc -l /scratch/jc33471/canine_tumor/results/Germline/PRJNA489159/${case}/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName

done 

for case in "CMT-100" "CMT-102" "CMT-103" "CMT-105" "CMT-106" "CMT-107" "CMT-109" "CMT-111" "CMT-112" "CMT-114"; do
    wc -l /project/szlab/Kun_Lin/Pan_Cancer/Mammary_Cancer/DepthOfCoverage/${case}/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName
    wc -l /scratch/jc33471/canine_tumor/results/DepthOfCoverage/PRJNA489159/${case}/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName

done
diff -i -E -Z -w -y /project/szlab/Kun_Lin/Pan_Cancer/Mammary_Cancer/Germline/CMT-100/SRR7780976_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName \
    /scratch/jc33471/canine_tumor/results/Germline/PRJNA489159/CMT-100/SRR7780976_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName

####### phylogenetics ########
conda activate phylogenetics
project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor_test/breed_prediction"
cd $run_dir

python ${project_dir}/scripts/phylogenetic/vaf2fasta.py \
    -i PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy.txt \
    --output-folder /scratch/jc33471/canine_tumor_test/breed_prediction \
    --resolve-IUPAC

python ${project_dir}/scripts/phylogenetic/vaf2fasta.py \
    -i PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy.txt \
    --output-folder /scratch/jc33471/canine_tumor_test/breed_prediction \
    --output-prefix "breed_specific" \
    --select_sites output_exclude_WGS/57_WGS_all_breed_specific_variants.txt \
    --resolve-IUPAC

awk -F, 'BEGIN{split("Shih Tzu,Schnauzer,Golden Retriever,Rottweiler,Greyhound,Maltese,Yorkshire Terrier,Boxer,Poodle,Cocker Spaniel", breed, ",")}{
    for (i in breed) if ( $5=="Normal" && $7==breed[i]) print $2;
        }'\
    ${project_dir}/metadata/data_collection_old.csv > ${run_dir}/breed_sample.list

# samples with breed info, breed specific sites
seqkit grep -n -f ${run_dir}/breed_sample.list ${run_dir}/breed_specific.min4.fasta | seqkit rmdup -n > ${run_dir}/breed_specific_breed_sample.min4.fasta

# samples with breed info, all sites
seqkit grep -n -f ${run_dir}/breed_sample.list ${run_dir}/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy.txt.min4.fasta | seqkit rmdup -n > ${run_dir}/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy_breed_sample.min4.fasta

