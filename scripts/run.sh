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

# list of samples with selected pure breed info
awk -F, 'BEGIN{split("Shih Tzu,Schnauzer,Golden Retriever,Rottweiler,Greyhound,Maltese,Yorkshire Terrier,Boxer,Poodle,Cocker Spaniel", breed, ",")}{
    for (i in breed) if ( $5=="Normal" && $7==breed[i] && $10~/Pass QC/) print $2;
        }'\
    ${project_dir}/metadata/data_collection_old.csv > ${run_dir}/breed_sample.list

# list of samples with selected pure breed info + missing
awk -F, 'BEGIN{split("Shih Tzu,Schnauzer,Golden Retriever,Rottweiler,Greyhound,Maltese,Yorkshire Terrier,Boxer,Poodle,Cocker Spaniel,No breed provided", breed, ",")}{
    for (i in breed) if ( $5=="Normal" && $7==breed[i] && $10~/Pass QC/) print $2;
        }'\
    ${project_dir}/metadata/data_collection_old.csv > ${run_dir}/breed_plus_unknown_sample.list

# convert MAF matrix to sequence, including all samples & all sites
python ${project_dir}/scripts/phylogenetic/vaf2fasta.py \
    -i ${run_dir}/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage.txt.gz \
    --output-folder ${run_dir} \
    --output-prefix "all_sample_site" \
    --resolve-IUPAC

# samples with breed info, all sites
seqkit grep -n -f ${run_dir}/breed_sample.list ${run_dir}/all_sample_site.min4.fasta | seqkit rmdup -n > ${run_dir}/all_breedSample_site.min4.fasta

# samples with breed+missing info, all sites
seqkit grep -n -f ${run_dir}/breed_plus_unknown_sample.list ${run_dir}/all_sample_site.min4.fasta | seqkit rmdup -n > ${run_dir}/all_breedPlusMissingSample_site.min4.fasta


# convert MAF matrix to sequence, including all samples & breed specific sites
python ${project_dir}/scripts/phylogenetic/vaf2fasta.py \
    -i ${run_dir}/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage.txt.gz \
    --output-folder ${run_dir} \
    --output-prefix "all_sample_breedSpecific" \
    --select_sites ${run_dir}/output_exclude_WGS/57_WGS_all_breed_specific_variants.txt \
    --resolve-IUPAC

# samples with breed info, breed specific sites
seqkit grep -n -f ${run_dir}/breed_sample.list ${run_dir}/all_sample_breedSpecific.min4.fasta | seqkit rmdup -n > ${run_dir}/all_breedSample_breedSpecific.min4.fasta

# samples with breed+missing info, breed specific sites
seqkit grep -n -f ${run_dir}/breed_plus_unknown_sample.list ${run_dir}/all_sample_breedSpecific.min4.fasta | seqkit rmdup -n > ${run_dir}/all_breedPlusMissingSample_breedSpecific.min4.fasta

### generate NJ tree and figures
# samples with breed info, all sites
Rscript --vanilla ${project_dir}/scripts/phylogenetic/nj_tree.R \
    $SLURM_NTASKS \
    ${run_dir}/all_breedSample_site.min4.fasta \
    "all_breedSample_site" \
    ${run_dir} \
    ${project_dir}/metadata/data_collection_old.csv

# samples with breed+missing info, all sites
Rscript --vanilla ${project_dir}/scripts/phylogenetic/nj_tree.R \
    $SLURM_NTASKS \
    ${run_dir}/all_breedPlusMissingSample_site.min4.fasta \
    "all_breedPlusMissingSample_site" \
    ${run_dir} \
    ${project_dir}/metadata/data_collection_old.csv

# samples with breed info, breed specific sites
Rscript --vanilla ${project_dir}/scripts/phylogenetic/nj_tree.R \
    $SLURM_NTASKS \
    ${run_dir}/all_breedSample_breedSpecific.min4.fasta \
    "all_breedSample_breedSpecific" \
    ${run_dir} \
    ${project_dir}/metadata/data_collection_old.csv

# samples with breed+missing info, breed specific sites
Rscript --vanilla ${project_dir}/scripts/phylogenetic/nj_tree.R \
    $SLURM_NTASKS \
    ${run_dir}/all_breedPlusMissingSample_breedSpecific.min4.fasta \
    "all_breedPlusMissingSample_breedSpecific" \
    ${run_dir} \
    ${project_dir}/metadata/data_collection_old.csv

awk '{if($1=="Gene" || $1=="DLA-12") print}' ${run_dir}/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy.txt > ${run_dir}/DLA-12_maf.txt
awk '{if($1=="Gene" || $1=="DLA88") print}' ${run_dir}/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy.txt > ${run_dir}/DLA88_maf.txt
python ${project_dir}/scripts/phylogenetic/vaf2fasta.py \
    -i ${run_dir}/DLA-12_maf.txt \
    --output-folder ${run_dir} \
    --output-prefix "DLA-12" \
    --resolve-IUPAC
seqkit grep -n -f ${run_dir}/breed_sample.list ${run_dir}/DLA-12.min4.fasta | seqkit rmdup -n > ${run_dir}/DLA-12_breedSample.min4.fasta

python ${project_dir}/scripts/phylogenetic/vaf2fasta.py \
    -i ${run_dir}/DLA88_maf.txt \
    --output-folder ${run_dir} \
    --output-prefix "DLA88" \
    --resolve-IUPAC
seqkit grep -n -f ${run_dir}/breed_sample.list ${run_dir}/DLA88.min4.fasta | seqkit rmdup -n > ${run_dir}/DLA88_breedSample.min4.fasta

# list of site for DLA-12
awk  '{if($1=="Gene" || $1=="DLA-12") print}' ${run_dir}/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy.txt | \
    awk 'BEGIN{OFS="\t"}{gene=$1;locus=$2":"$3;mutation=$4">"$5; print gene,locus,mutation}' \
    > ${run_dir}/DLA-12_sites.txt
awk  '{if($1=="Gene" || $1=="DLA88") print}' ${run_dir}/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy.txt | \
    awk 'BEGIN{OFS="\t"}{gene=$1;locus=$2":"$3;mutation=$4">"$5; print gene,locus,mutation}' \
    > ${run_dir}/DLA88_sites.txt

### including low depth variants
python ${project_dir}/scripts/phylogenetic/vaf2fasta.py \
    -i ${run_dir}/PanCancer_disc_val_merged_germline_VAF_01_01_2021.txt.gz \
    --output-folder ${run_dir} \
    --output-prefix "DLA-12" \
    --select_sites ${run_dir}/DLA-12_sites.txt \
    --resolve-IUPAC
seqkit grep -n -f ${run_dir}/breed_sample.list ${run_dir}/DLA-12.min4.fasta | seqkit rmdup -n > ${run_dir}/DLA-12_breedSample.min4.fasta

python ${project_dir}/scripts/phylogenetic/vaf2fasta.py \
    -i ${run_dir}/PanCancer_disc_val_merged_germline_VAF_01_01_2021.txt.gz \
    --output-folder ${run_dir} \
    --output-prefix "DLA88" \
    --select_sites ${run_dir}/DLA88_sites.txt \
    --resolve-IUPAC
seqkit grep -n -f ${run_dir}/breed_sample.list ${run_dir}/DLA88.min4.fasta | seqkit rmdup -n > ${run_dir}/DLA88_breedSample.min4.fasta

# DRB1 no variants???
awk '{if($1=="Gene" || $1=="DRB1") print}' ${run_dir}/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy.txt > ${run_dir}/DRB1_maf.txt
cut -f1 PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy.txt | grep "DRB1" | sort | uniq -c


# DQA1
cut -f1 PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy.txt | grep -E "LA-DRB1|DQA1|DQB1" | sort | uniq -c
awk '{if($1=="Gene" || $1=="DLA-DQA1") print}' ${run_dir}/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy.txt > ${run_dir}/DLA-DQA1_maf.txt
python ${project_dir}/scripts/phylogenetic/vaf2fasta.py \
    -i ${run_dir}/DLA-DQA1_maf.txt \
    --output-folder ${run_dir} \
    --output-prefix "DLA-DQA1" \
    --resolve-IUPAC
seqkit grep -n -f ${run_dir}/breed_sample.list ${run_dir}/DLA-DQA1.min4.fasta | seqkit rmdup -n > ${run_dir}/DLA-DQA1_breedSample.min4.fasta

Rscript --vanilla ${project_dir}/scripts/phylogenetic/nj_tree.R \
    $SLURM_NTASKS \
    ${run_dir}/DLA-DQA1_breedSample.min4.fasta \
    "DLA-DQA1_breedSample" \
    ${run_dir} \
    ${project_dir}/metadata/data_collection_old.csv


# DQB1
awk '{if($1=="Gene" || $1=="HLA-DQB1") print}' ${run_dir}/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy.txt > ${run_dir}/HLA-DQB1_maf.txt
python ${project_dir}/scripts/phylogenetic/vaf2fasta.py \
    -i ${run_dir}/HLA-DQB1_maf.txt \
    --output-folder ${run_dir} \
    --output-prefix "HLA-DQB1" \
    --resolve-IUPAC
seqkit grep -n -f ${run_dir}/breed_sample.list ${run_dir}/HLA-DQB1.min4.fasta | seqkit rmdup -n > ${run_dir}/HLA-DQB1_breedSample.min4.fasta

Rscript --vanilla ${project_dir}/scripts/phylogenetic/nj_tree.R \
    $SLURM_NTASKS \
    ${run_dir}/HLA-DQB1_breedSample.min4.fasta \
    "HLA-DQB1_breedSample" \
    ${run_dir} \
    ${project_dir}/metadata/data_collection_old.csv

