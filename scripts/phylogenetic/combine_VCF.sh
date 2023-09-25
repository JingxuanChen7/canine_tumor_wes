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
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/Mammary_Cancer/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf > /scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/ValidationData/MC/store/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> /scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/Melanoma/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> /scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/osteosarcoma/results/normal/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> /scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/osteosarcoma/results/tumor/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> /scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/lymphoma/results/normal/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> /scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/lymphoma/results/tumor/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> /scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/other/results/normal/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> /scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list
ls -1 /project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/other/results/tumor/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> /scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/Glioma/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> /scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list
ls -1 /project/szlab/Kun_Lin/Pan_Cancer/ValidationData/HSA/Germline/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> /scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list
ls -1 /project/szlab/tw71066/8-25-20-scratch-backup/All_NEW_WES/PRJNA552034-HSA/Mutation/new_mutect_results/*/*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf >> /scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list
# copy backup results to scratch
rsync -av --files-from=/scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/backup_vcf_file.list / /scratch/jc33471/canine_tumor/phylogenetics/vcf/
# file path/name on scratch
awk '{print "/scratch/jc33471/canine_tumor/phylogenetics/vcf"$0}' ${run_dir}/merge_vcf/backup_vcf_file.list > ${run_dir}/merge_vcf/old_merged_germline_variants.list
awk '{print "/scratch/jc33471/canine_tumor/phylogenetics/vcf"$0".reheader.gz"}' ${run_dir}/merge_vcf/backup_vcf_file.list > ${run_dir}/merge_vcf/old_merged_germline_variants_gz.list

#### on work ####
cd ${run_dir}/merge_vcf
while read file; do
    rm -f ${file}".gz" ${file}".gz.csi"
    sample_name=`basename ${file} "_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf"`
    echo "${sample_name}" > ${file}"rename_samples_vcf.txt"
    bcftools reheader -s ${file}"rename_samples_vcf.txt" --threads ${SLURM_NTASKS} -o ${file}".reheader" ${file}
    bgzip --force --threads ${SLURM_NTASKS} ${file}".reheader"
    bcftools index --threads ${SLURM_NTASKS} ${file}".reheader.gz"
    echo "Compressed and indexed ${file}"
done < ${run_dir}/merge_vcf/old_merged_germline_variants.list

bcftools merge --threads ${SLURM_NTASKS} -O z \
    -o ${run_dir}/merge_vcf/old_merged_germline_variants.vcf.gz \
    --file-list ${run_dir}/merge_vcf/old_merged_germline_variants_gz.list

