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
#SBATCH --output=/scratch/jc33471/canine_tumor/wgs_breed_prediction/vcf/temp.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate phylogenetics

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor/wgs_breed_prediction"
breed_dir="/home/${USER}/breed_prediction"
mkdir -p ${run_dir}
cd ${run_dir}

mkdir -p vcf/ merge_vcf/

cd $run_dir/vcf

# dog10K
# wget -O $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs.vcf.gz "https://kiddlabshare.med.umich.edu/dog10K/SNP_and_indel_calls_2021-10-17/AutoAndXPAR.SNPs.vqsr99.vcf.gz"

# test dataset
zcat $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs.vcf.gz | head -n 3000 | bgzip -c > $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_test.vcf.gz
bcftools query -l $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs.vcf.gz > $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs.sample.txt
bcftools filter --threads $SLURM_NTASKS \
    -e ' FS > 30 || QD < 2 ' --output-type z \
    -o $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered.vcf.gz \
    $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs.vcf.gz

# set coverage < 10 to missing genotype, and calculate the fraction of missing
bcftools +setGT --threads $SLURM_NTASKS $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered.vcf.gz -- -t q -i 'FMT/DP<10' -n ./. |\
    bcftools +fill-tags - -O z --threads $SLURM_NTASKS -o $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added.vcf.gz -- -t INFO/F_MISSING

# filter VCF by fraction of missing
bcftools filter --threads $SLURM_NTASKS \
    -e ' F_MISSING > 0.2 ' --output-type z \
    -o $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing.vcf.gz \
    $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added.vcf.gz

# bcftools +fill-tags --threads $SLURM_NTASKS $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing.vcf.gz -- -t FORMAT/VAF |\
#     bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%VAF]\n' -o $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing.vaf.txt

# clean up chromosomes, making it consistent with canfam3 reference
grep ">" /work/szlab/Lab_shared_PanCancer/source/canFam3.fa | sed 's/>//g' | awk 'BEGIN{ORS=","}{print}' | sed 's/,$//g' > $run_dir/vcf/canfam3_chr.txt

bcftools index $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing.vcf.gz
bcftools view --threads $SLURM_NTASKS \
    -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29,chr30,chr31,chr32,chr33,chr34,chr35,chr36,chr37,chr38,chrX \
    -O z \
    -o $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_chr.vcf.gz \
    $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing.vcf.gz

# canfam4->canfam3 liftover
conda activate wes_env
# wget -O ${run_dir}/vcf/canFam4ToCanFam3.over.chain.gz "http://hgdownload.soe.ucsc.edu/goldenPath/canFam4/liftOver/canFam4ToCanFam3.over.chain.gz"
picard -Xmx60G LiftoverVcf \
    I=${run_dir}/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_chr.vcf.gz \
    O=${run_dir}/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_chr_liftover.vcf.gz \
    CHAIN=${run_dir}/vcf/canFam4ToCanFam3.over.chain.gz \
    REJECT=${run_dir}/vcf/rejected_variants.vcf \
    R="/work/szlab/Lab_shared_PanCancer/source/canFam3.fa"

# bcftools +liftover \
#   POPRES_Genotypes_QC2_v2_VCF.vcf \
#   -- \
#   -s hg18.fa \
#   -f hg38.fa \
#   -c hg18ToHg38.over.chain.gz
