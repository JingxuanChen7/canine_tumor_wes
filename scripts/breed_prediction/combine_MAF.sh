#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=breed
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --tasks-per-node=8
#SBATCH --mem=60G
#SBATCH --time=100:00:00
#SBATCH --mail-user=jc33471@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/scratch/jc33471/canine_tumor/breed_prediction/breed.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate breed_prediction

project_dir="/home/${USER}/canine_tumor_wes"
breed_dir="/home/${USER}/breed_prediction"
run_dir="/scratch/${USER}/canine_tumor/breed_prediction"
mkdir -p ${run_dir}
cd ${run_dir}

#### prepare metadata for this run ####
metadata=${run_dir}/this_paired_meta.csv
metadata_wunpaired=${run_dir}/this_wunpaired_meta.csv

# only keep paired samples, note that doulble tumor/normal are excluded
grep -f <(cut -d, -f1 ${project_dir}/metadata/data_collection_new_renamed.csv | sort | uniq -d) ${project_dir}/metadata/data_collection_new_renamed.csv |\
  grep -v -f <(cut -d, -f1,5 ${project_dir}/metadata/data_collection_new_renamed.csv | sort | uniq -d | cut -f1 -d,) > ${run_dir}/paired_new_meta.csv
# tumor or normal only samples
grep -f <(cut -d, -f1 ${project_dir}/metadata/data_collection_new_renamed.csv | sort | uniq -u) ${project_dir}/metadata/data_collection_new_renamed.csv > ${run_dir}/unpaired_new_meta.csv

cat ${project_dir}/metadata/data_collection_old.csv ${run_dir}/paired_new_meta.csv  > ${metadata}
cat ${project_dir}/metadata/data_collection_old.csv ${run_dir}/paired_new_meta.csv ${run_dir}/unpaired_new_meta.csv > ${metadata_wunpaired}

##### create input file list #####
## backup file list (runs before 2021)
backup="/home/jc33471/canine_tumor_wes/metadata/backup_vcf_file.list"
## annovar
grep "_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName" ${backup} |\
  awk '{print "/scratch/jc33471/canine_tumor/phylogenetics/vcf"$0}' \
  > ${run_dir}/backup_annovar.txt
# PRJNA525883 re-run
ls -1 /scratch/jc33471/canine_tumor_rerun/results/Germline/PRJNA525883/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName > ${run_dir}/rerun_annovar.txt
# new datasets
ls -1 /scratch/jc33471/canine_tumor_0908/results/Germline/*/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName > ${run_dir}/newdata_annovar.txt
cat ${run_dir}/backup_annovar.txt ${run_dir}/rerun_annovar.txt ${run_dir}/newdata_annovar.txt > ${run_dir}/all_annovar.txt
## depth of coverage
grep "_DepthofCoverage_CDS.bed" ${backup} |\
  awk '{print "/scratch/jc33471/canine_tumor/phylogenetics/vcf"$0}' \
  > ${run_dir}/backup_depthofcoverage.txt
# PRJNA525883 re-run
ls -1 /scratch/jc33471/canine_tumor_rerun/results/DepthOfCoverage/PRJNA525883/*/*_DepthofCoverage_CDS.bed > ${run_dir}/rerun_depthofcoverage.txt
# new datasets
ls -1 /scratch/jc33471/canine_tumor_0908/results/DepthOfCoverage/*/*/*_DepthofCoverage_CDS.bed > ${run_dir}/newdata_depthofcoverage.txt
cat ${run_dir}/backup_depthofcoverage.txt ${run_dir}/rerun_depthofcoverage.txt ${run_dir}/newdata_depthofcoverage.txt > ${run_dir}/all_depthofcoverage.txt
## vcf
grep -E "_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf$" ${backup} |\
  awk '{print "/scratch/jc33471/canine_tumor/phylogenetics/vcf"$0}' \
  > ${run_dir}/backup_vcf.txt
# PRJNA525883 re-run
ls -1 /scratch/jc33471/canine_tumor_rerun/results/Germline/PRJNA525883/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf > ${run_dir}/rerun_vcf.txt
# new datasets
ls -1 /scratch/jc33471/canine_tumor_0908/results/Germline/*/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf > ${run_dir}/newdata_vcf.txt

cat ${run_dir}/backup_vcf.txt ${run_dir}/rerun_vcf.txt ${run_dir}/newdata_vcf.txt > ${run_dir}/all_vcf.txt

## reformat annovar file list, only paired
python ${project_dir}/scripts/breed_prediction/create_file_list.py \
  --metatable ${metadata} \
  --filelist ${run_dir}/all_annovar.txt \
  --vcflist ${run_dir}/all_vcf.txt \
  --depthlist ${run_dir}/all_depthofcoverage.txt \
  --output ${run_dir}/annovar_file_list_paired.txt \
  --vcfout ${run_dir}/vcf_file_list_paired.txt

## reformat annovar file list, with unpaired
python ${project_dir}/scripts/breed_prediction/create_file_list.py \
  --metatable ${metadata_wunpaired} \
  --filelist ${run_dir}/all_annovar.txt \
  --vcflist ${run_dir}/all_vcf.txt \
  --depthlist ${run_dir}/all_depthofcoverage.txt \
  --output ${run_dir}/annovar_file_list.txt \
  --vcfout ${run_dir}/vcf_file_list.txt

######## expected input arguments ########
# A file with one row (three columns) per normal/tumor sample including file path for germline variants called by HaplotypeCaller and annotated by ANNOVAR and Pancancer.jar (see sample files)
annovar_file_list=${run_dir}/annovar_file_list_paired.txt
# A file with one row (three columns) per normal/tumor sample including file paths for VCF format files for germline variants called by HaplotypeCaller and base coverage information output by DepthOfCoverage (see sample files)
vcf_file_list=${run_dir}/vcf_file_list.txt

# output file generated by util.CreateVariantList with a combined list of all germline variants (one row, seven columns per variant) observed in both normal and tumor of at least one case (see sample file)
output_variant_file=${run_dir}/variant_list.txt
# output file generated by util.CreateVariantList with a combined list of all germline variants (one row per variant) and the VAF in each normal and tumor sample (one column per sample, see sample files)
output_VAF_matrix=${run_dir}/germline_VAF_matrix.txt
# output file generated by util.CreateVariantList and updated by util.AddDepthOfCoverage with a combined list of all germline variants (one row per variant) and the depth values (total depth from DepthOfCoverage, total depth DP and allele depth AD from HaplotypeCaller) in each normal and tumor sample (one column per sample, see sample files)
output_depths_matrix=${run_dir}/germline_depth_matrix.txt
# output file generated by util.CreateVariantList and updated by germline_VAF_reset_low_coverage.R 
out_reset_low_coverage=${run_dir}/germline_VAF_matrix.reset_low_coverage.txt

# Create a new variant file list
echo "Creating new variant list"
java -Xmx60g -cp ${breed_dir}/Pancancer.jar util.CreateVariantList ${annovar_file_list} ${output_variant_file}.list include_synonymous

# Sort and add a header line
sort ${output_variant_file}.list > ${output_variant_file}.list.sorted
printf "Gene\tChromosome\tPosition\tRef\tAlt\tProtein\tMutation_type\n" > ${output_variant_file}.header
cat ${output_variant_file}.header ${output_variant_file}.list.sorted > ${output_variant_file}.sorted
rm ${output_variant_file}.header ${output_variant_file}.list ${output_variant_file}.list.sorted

# create variant VAF and depths matrix files
echo ""
echo "Creating variant VAF matrix"
java -Xmx60g -cp ${breed_dir}/Pancancer.jar util.GetVariantInfoFromVCF ${output_variant_file}.sorted 7 ${vcf_file_list} ${output_VAF_matrix} VAF

echo ""
echo "Creating variant depths matrix"
java -Xmx60g -cp ${breed_dir}/Pancancer.jar util.GetVariantInfoFromVCF ${output_variant_file}.sorted 7 ${vcf_file_list} ${output_depths_matrix} depths

echo ""
echo "Adding depth from DepthOfCoverage bed files"
java -Xmx60g -cp ${breed_dir}/Pancancer.jar util.AddDepthOfCoverage ${output_depths_matrix} 7 ${vcf_file_list} overwrite

gzip ${output_VAF_matrix}
gzip ${output_depths_matrix}

# convert VAF at low coverage regions to NA
Rscript --vanilla ${project_dir}/scripts/breed_prediction/germline_VAF_reset_low_coverage.R \
  ${output_VAF_matrix}".gz"\
  ${output_depths_matrix}".gz" \
  ${out_reset_low_coverage}".gz"


# identify breed specific variants
# prepare metadata
Rscript --vanilla ${project_dir}/scripts/breed_prediction/reformat_meta_text.R \
  ${metadata_wunpaired} \
  ${run_dir}/breed_prediction_metadata.txt \
  ${run_dir}/breed_prediction_cases.txt

# awk 'BEGIN{FS=OFS="\t"}{if ($6 ~ "Pass QC" && $5 == "Normal") print $4}' ${run_dir}/breed_prediction_cases.txt | sort | uniq -c | sort -n

mkdir -p ${run_dir}"/output_exclude_WGS"

# Rscript --vanilla ${project_dir}/scripts/breed_prediction/breed_specific_variants.R \
#   ${project_dir}/scripts/breed_prediction/build_sample_meta_data.R \
#   ${out_reset_low_coverage}".gz" \
#   ${run_dir}"/output_exclude_WGS/breed_unique_variants_10.txt" \
#   ${run_dir}"/output_exclude_WGS/breed_enriched_variants_10.txt" \
#   ${run_dir}"/output_exclude_WGS/all_breed_specific_variants_10.txt" \
#   ${run_dir}"/breed_prediction_metadata.txt" \
#   "Yorkshire Terrier" "Shih Tzu" "Poodle" "Maltese" "Rottweiler" "Greyhound" "Golden Retriever" "Cocker Spaniel" "Boxer" "Schnauzer"

Rscript --vanilla ${project_dir}/scripts/breed_prediction/breed_specific_variants.R \
  ${project_dir}/scripts/breed_prediction/build_sample_meta_data.R \
  ${out_reset_low_coverage}".gz" \
  ${run_dir}"/output_exclude_WGS/breed_unique_variants.txt" \
  ${run_dir}"/output_exclude_WGS/breed_enriched_variants.txt" \
  ${run_dir}"/output_exclude_WGS/all_breed_specific_variants.txt" \
  ${run_dir}"/breed_prediction_metadata.txt" \
  "Yorkshire Terrier" "Shih Tzu" "Poodle" "Maltese" "Rottweiler" "Greyhound" "Golden Retriever" "Cocker Spaniel" "Boxer" "Schnauzer" "Labrador Retriever" "Boston Terrier"
  