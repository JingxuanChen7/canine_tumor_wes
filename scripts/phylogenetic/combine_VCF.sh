#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=merge_vcf
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --tasks-per-node=8
#SBATCH --mem=60G
#SBATCH --time=500:00:00
#SBATCH --mail-user=jc33471@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/scratch/jc33471/canine_tumor/phylogenetics/merge_vcf.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate phylogenetics

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor/phylogenetics"
mkdir -p ${run_dir}
cd ${run_dir}

mkdir -p vcf/ merge_vcf/


#### on xfer ####
# list of files on project
out_list="/scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list"
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/Mammary_Cancer/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf > ${out_list}
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/ValidationData/MC/store/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> ${out_list}
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/Melanoma/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/osteosarcoma/results/normal/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/osteosarcoma/results/tumor/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/lymphoma/results/normal/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/lymphoma/results/tumor/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/other/results/normal/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/other/results/tumor/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> ${out_list}
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/Glioma/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> ${out_list}
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/ValidationData/HSA/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> ${out_list}
ls -1 /project/szlab/tw71066/8-25-20-scratch-backup/All_NEW_WES/PRJNA552034-HSA/Mutation/new_mutect_results/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> ${out_list}
# Annovar output with gene names
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/Mammary_Cancer/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName >> ${out_list}
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/ValidationData/MC/store/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName >> ${out_list}
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/Melanoma/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/osteosarcoma/results/normal/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/osteosarcoma/results/tumor/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/lymphoma/results/normal/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/lymphoma/results/tumor/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/other/results/normal/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/other/results/tumor/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName >> ${out_list}
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/Glioma/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName >> ${out_list}
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/ValidationData/HSA/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName >> ${out_list}
ls -1 /project/szlab/tw71066/8-25-20-scratch-backup/All_NEW_WES/PRJNA552034-HSA/Mutation/new_mutect_results/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName >> ${out_list}
# depthOfCoverage out bed files
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/Mammary_Cancer/DepthOfCoverage/*_DepthofCoverage_CDS.bed >> ${out_list}
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/ValidationData/MC/store/DepthOfCoverage/*/*_DepthofCoverage_CDS.bed >> ${out_list}
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/Melanoma/DepthOfCoverage/*_DepthofCoverage_CDS.bed >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/osteosarcoma/results/normal/*/*_DepthofCoverage_CDS.bed >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/osteosarcoma/results/tumor/*/*_DepthofCoverage_CDS.bed >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/lymphoma/results/normal/*/*_DepthofCoverage_CDS.bed >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/lymphoma/results/tumor/*/*_DepthofCoverage_CDS.bed >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/other/results/normal/*/*_DepthofCoverage_CDS.bed >> ${out_list}
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/other/results/tumor/*/*_DepthofCoverage_CDS.bed >> ${out_list}
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/Glioma/DepthOfCoverage/*/*_DepthofCoverage_CDS.bed >> ${out_list}
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/ValidationData/HSA/DepthOfCoverage/*/*_DepthofCoverage_CDS.bed >> ${out_list}
ls -1 /project/szlab/tw71066/8-25-20-scratch-backup/All_NEW_WES/PRJNA552034-HSA/Mutation/new_mutect_results/*/*_DepthofCoverage_CDS.bed >> ${out_list}

# copy backup results to scratch
rsync -av --files-from=/scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list / /scratch/jc33471/canine_tumor/phylogenetics/vcf/


# file path/name on scratch
awk '{print "/scratch/jc33471/canine_tumor/phylogenetics/vcf"$0}' ${run_dir}/merge_vcf/backup_vcf_file.list > ${run_dir}/merge_vcf/old_merged_germline_variants.list
awk '{print "/scratch/jc33471/canine_tumor/phylogenetics/vcf"$0".reheader.gz"}' ${run_dir}/merge_vcf/backup_vcf_file.list > ${run_dir}/merge_vcf/old_merged_germline_variants_gz.list

# PRJNA525883 was not backedup, re-run
ls -1 /scratch/jc33471/canine_tumor_rerun/results/Germline/PRJNA525883/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf > ${run_dir}/merge_vcf/rerun_merged_germline_variants.list
awk '{print $0".reheader.gz"}' ${run_dir}/merge_vcf/rerun_merged_germline_variants.list > ${run_dir}/merge_vcf/rerun_merged_germline_variants_gz.list

#### on work ####
cd ${run_dir}/merge_vcf
while read file; do
    rm -f ${file}".gz" ${file}".gz.csi"
    sample_name=`basename ${file} "_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf"`
    echo "${sample_name}" > ${file}"rename_samples_vcf.txt"
    bcftools reheader -s ${file}"rename_samples_vcf.txt" --threads ${SLURM_NTASKS} ${file} | awk '{ if($1 ~ /##/ || NF==10)print}' > ${file}".reheader"
    bgzip --force --threads ${SLURM_NTASKS} ${file}".reheader"
    bcftools index --threads ${SLURM_NTASKS} ${file}".reheader.gz"
    echo "Compressed and indexed ${file}"
done < ${run_dir}/merge_vcf/old_merged_germline_variants.list

while read file; do
    rm -f ${file}".gz" ${file}".gz.csi"
    sample_name=`basename ${file} "_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf"`
    echo "${sample_name}" > ${file}"rename_samples_vcf.txt"
    bcftools reheader -s ${file}"rename_samples_vcf.txt" --threads ${SLURM_NTASKS} ${file} | awk '{ if($1 ~ /##/ || NF==10)print}' > ${file}".reheader"
    bgzip --force --threads ${SLURM_NTASKS} ${file}".reheader"
    bcftools index --threads ${SLURM_NTASKS} ${file}".reheader.gz"
    echo "Compressed and indexed ${file}"
done < ${run_dir}/merge_vcf/rerun_merged_germline_variants.list

# merge file lists from different run batches
# if other dataset should be added, combine here
cat ${run_dir}/merge_vcf/old_merged_germline_variants_gz.list \
    ${run_dir}/merge_vcf/rerun_merged_germline_variants_gz.list \
    > ${run_dir}/merge_vcf/all_germline_vairants_to_merge.list

# merge all runs
bcftools merge --threads ${SLURM_NTASKS} -O z --force-samples  --missing-to-ref \
    -o ${run_dir}/merge_vcf/all_merged_germline_variants.vcf.gz \
    --file-list ${run_dir}/merge_vcf/all_germline_vairants_to_merge.list

# only keep SNPs, only keep PASSed variants
bcftools view --threads ${SLURM_NTASKS} \
    -f PASS --types snps -O z \
    -o ${run_dir}/merge_vcf/all_merged_germline_variants_SNP.vcf.gz \
    ${run_dir}/merge_vcf/all_merged_germline_variants.vcf.gz


# select sites and/or samples from the VCF
# list of samples with selected pure breed info
awk -F, 'BEGIN{split("Shih Tzu,Schnauzer,Golden Retriever,Rottweiler,Greyhound,Maltese,Yorkshire Terrier,Boxer,Poodle,Cocker Spaniel", breed, ",")}{
    for (i in breed) if ( $5~"Normal" && $7~breed[i] && $10~/Pass QC/) print $2;
        }'\
    ${project_dir}/metadata/data_collection_old.csv | sed 's/\"//g' > ${run_dir}/merge_vcf/breed_sample.list

# list of samples with selected pure breed info + missing
awk -F, 'BEGIN{split("Shih Tzu,Schnauzer,Golden Retriever,Rottweiler,Greyhound,Maltese,Yorkshire Terrier,Boxer,Poodle,Cocker Spaniel,No breed provided", breed, ",")}{
    for (i in breed) if ( $5~"Normal" && $7~breed[i] && $10~/Pass QC/) print $2;
        }'\
    ${project_dir}/metadata/data_collection_old.csv | sed 's/\"//g' > ${run_dir}/merge_vcf/breed_plus_unknown_sample.list

# samples with selected pure breed info
bcftools view --threads ${SLURM_NTASKS} --force-samples \
    -S ${run_dir}/merge_vcf/breed_sample.list -O z \
    -o ${run_dir}/merge_vcf/all_merged_germline_variants_SNP_breedSample.vcf.gz \
    ${run_dir}/merge_vcf/all_merged_germline_variants_SNP.vcf.gz


# samples with selected pure breed info + missing
bcftools view --threads ${SLURM_NTASKS} --force-samplesx \
    -S ${run_dir}/merge_vcf/breed_plus_unknown_sample.list -O z \
    -o ${run_dir}/merge_vcf/all_merged_germline_variants_SNP_breedPlusMissingSample.vcf.gz \
    ${run_dir}/merge_vcf/all_merged_germline_variants_SNP.vcf.gz

# select breed specific sites
breedSpecific_sites="/scratch/jc33471/canine_tumor_test/breed_prediction/output_exclude_WGS/57_WGS_all_breed_specific_variants.txt"
awk 'BEGIN{OFS = "\t"}{if ($2 !~ /Locus/) {split($2,POS,":"); print POS[1],POS[2];} }' \
    ${breedSpecific_sites} > ${run_dir}/merge_vcf/breed_specific_sites.tab


# samples with breed info, breed specific sites
bcftools view --threads ${SLURM_NTASKS}  \
    -R ${run_dir}/merge_vcf/breed_specific_sites.tab -O z \
    -o ${run_dir}/merge_vcf/all_merged_germline_variants_SNP_breedSample_breedSpecific.vcf.gz \
    ${run_dir}/merge_vcf/all_merged_germline_variants_SNP_breedSample.vcf.gz

# samples with breed+missing, breed specific sites
bcftools view --threads ${SLURM_NTASKS}  \
    -R ${run_dir}/merge_vcf/breed_specific_sites.tab -O z \
    -o ${run_dir}/merge_vcf/all_merged_germline_variants_SNP_breedPlusMissingSample_breedSpecific.vcf.gz \
    ${run_dir}/merge_vcf/all_merged_germline_variants_SNP_breedPlusMissingSample.vcf.gz

# breed specific tree
# 25890926 SNPs, 397 dogs / 504 dogs
for prefix in "breedSample" "breedPlusMissingSample" "breedSample_breedSpecific" "breedPlusMissingSample_breedSpecific"; do
    if [[ ! -e ${run_dir}/merge_vcf/all_merged_germline_variants_SNP_${prefix}.mdist.gz ]]; then
        plink --threads ${SLURM_NTASKS} --memory "60000" --seed 123 --dog \
            --distance gz square '1-ibs' \
            --vcf ${run_dir}/merge_vcf/all_merged_germline_variants_SNP_${prefix}.vcf.gz \
            --out ${run_dir}/merge_vcf/all_merged_germline_variants_SNP_${prefix}
    fi
done


python ${project_dir}/scripts/phylogenetic/mdist2phydist.py \
    --mdist ${run_dir}/merge_vcf/breedSample_breedSpecific_merged_germline_variants_SNP.mdist.gz \
    --id ${run_dir}/merge_vcf/breedSample_breedSpecific_merged_germline_variants_SNP.mdist.id \
    --phydist ${run_dir}/merge_vcf/breedSample_breedSpecific_merged_germline_variants_SNP.phylip 
# neighbor < ${run_dir}/merge_vcf/all_merged_germline_variants_SNP_matrix.mdist.gz

cd ${run_dir}/merge_vcf
echo "breedSample_breedSpecific_merged_germline_variants_SNP.phylip" > input
echo "Y" >> input
rm out*
neighbor < input > output; mv outfile breedSample_breedSpecific.log; mv outtree breedSample_breedSpecific.nwk

# run vcf2phylip and phylip
python ${project_dir}/scripts/phylogenetic/vcf2philip.py \
    --input all_merged_germline_variants_SNP_breedSample_breedSpecific.vcf.gz \
    --fasta \
    --output-prefix "breedSample_breedSpecific" 


### MHCII genes region
# download annotation from refseq

wget -O ${run_dir}/merge_vcf/refGene.txt.gz "https://hgdownload.cse.ucsc.edu/goldenpath/canFam3/database/refGene.txt.gz"

wget -O ${run_dir}/merge_vcf/GCF_000002285.3_CanFam3.1_genomic.gtf.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9615/105/GCF_000002285.3_CanFam3.1/GCF_000002285.3_CanFam3.1_genomic.gtf.gz"
gunzip ${run_dir}/merge_vcf/GCF_000002285.3_CanFam3.1_genomic.gtf.gz

# DQA1 (ENSCAFG00000000812), DQB1 (ENSCAFG00000000814) and DRB1.
emsembl="/work/szlab/Lab_shared_PanCancer/source/Canis_familiaris.CanFam3.1.99.chr.gtf"
refseq="${run_dir}/merge_vcf/GCF_000002285.3_CanFam3.1_genomic.gtf"

grep -E 'gene_id "ENSCAFG00000000812"|gene_id "ENSCAFG00000000814"' ${emsembl} > ${run_dir}/merge_vcf/ensembl_MHCII.gtf
grep -E 'gene_id "HLA-DRB1"' ${refseq} > ${run_dir}/merge_vcf/refseq_MHCII.gtf
cat ${run_dir}/merge_vcf/ensembl_MHCII.gtf ${run_dir}/merge_vcf/refseq_MHCII.gtf |\
    awk '($3=="exon"){print}' |\
    sed 's/NC_006594.3/chr12/g' \
    > ${run_dir}/merge_vcf_old_test/threegenes_MHCII.gtf
awk 'BEGIN{FS=OFS="\t"}{$4=$4-1; print $1,$4,$5}' ${run_dir}/merge_vcf_old_test/threegenes_MHCII.gtf > ${run_dir}/merge_vcf_old_test/threegenes_MHCII.bed