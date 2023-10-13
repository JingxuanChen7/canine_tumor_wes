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
#SBATCH --output=/scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/merge_vcf.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate wes_env

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor/phylogenetics"
mkdir -p ${run_dir}
cd ${run_dir}

mkdir -p vcf/ merge_vcf/

config="/home/jc33471/canine_tumor_wes/scripts/phylogenetic/config.json"

snakemake \
    -p \
    --cores ${SLURM_NTASKS} \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --use-conda \
    --configfile ${config} \
    --snakefile "${project_dir}/scripts/phylogenetic/Snakefile"

breed_dir="/home/${USER}/breed_prediction"

# sites
zcat breedPlusMissingSample_breedSpecific_merged_germline_variants_SNP.vcf.gz |\
    awk 'BEGIN{OFS="\t"}/^[^#]/{print "NA",$1,$2,$4,$5}' |\
    sed '1 i\Gene\tChromosome\tPosition\tRef\tAlt' > test.txt

# get varaint sites
zcat breedPlusMissingSample_breedSpecific_merged_germline_variants_SNP.vcf.gz |\
    awk 'BEGIN{OFS=":"}/^[^#]/{print $1,$2}' > test.txt
# select variant sites for depth info
grep -f test.txt /scratch/jc33471/canine_tumor/phylogenetics/vcf/project/szlab/Kun_Lin/Pan_Cancer/ValidationData/MC/store/DepthOfCoverage/004/SRR9911377_DepthofCoverage_CDS.bed > test_DepthofCoverage_CDS.bed

# samples
grep -f breedSample.list ${project_dir}/metadata/vcf_file_list.txt > tmp_vcf_file_list.txt
# matrix for depth
java -Xmx60g -cp ${breed_dir}/Pancancer.jar util.GetVariantInfoFromVCF test.txt 5 tmp_vcf_file_list.txt germline_depth_matrix_for_phylo.txt depths
# mask vcf

