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
#SBATCH --output=/scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/test.out

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
    --breedSpecific "/scratch/jc33471/canine_tumor/breed_prediction/output_exclude_WGS/all_breed_specific_variants_13.txt" \
    --somaticMutation ${project_dir}"/metadata/Pass_QC_Final_Total_withGene_Burair_Filtering4_VAF_Mutect_orientBiasModified_04_02.txt" \
    --threads 8 \
    --memory "60G"

snakemake \
    -np \
    --latency-wait 60 \
    --cores ${SLURM_NTASKS} \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --use-conda \
    --configfile ${config} \
    --snakefile "${project_dir}/scripts/phylogenetic/Snakefile"


# # sites in order
# zcat breedPlusMissingSample_breedSpecific_merged_germline_variants_SNP.vcf.gz |\
#     awk 'BEGIN{OFS="\t"}/^[^#]/{print "NA",$1,$2,$4,$5}' |\
#     sed '1 i\Gene\tChromosome\tPosition\tRef\tAlt' > test.txt

# samples in order
grep -f breedSample.list ${project_dir}/metadata/vcf_file_list.txt |  awk 'NR==FNR { line[$1] = $0; next } $1 in line { print line[$1] }' - breedSample.list | cut -f3 > tmp_vcf_file_list.txt

# get varaint sites in order
zcat breedPlusMissingSample_breedSpecific_merged_germline_variants_SNP.vcf.gz |\
    awk 'BEGIN{OFS=":"}/^[^#]/{print $1,$2}' > test.txt

# select variant sites for depth info in order
while read file; do
    grep -f test.txt ${file} | awk 'NR==FNR { line[$1] = $0; next } $1 in line { print line[$1] }' - test.txt | cut -f2 > test_DepthofCoverage_CDS.txt
    if [[ ! -e germline_depth_matrix_for_phylo.txt ]]; then
        mv test_DepthofCoverage_CDS.txt germline_depth_matrix_for_phylo.txt
    else
        paste -d, germline_depth_matrix_for_phylo.txt test_DepthofCoverage_CDS.txt > tmp.txt
        mv tmp.txt germline_depth_matrix_for_phylo.txt
    fi
done < tmp_vcf_file_list.txt

# matrix for depth
file1="/scratch/jc33471/canine_tumor/phylogenetics/vcf/project/szlab/Kun_Lin/Pan_Cancer/ValidationData/MC/store/DepthOfCoverage/008/SRR9911360_DepthofCoverage_CDS.bed"
grep -f test.txt ${file1} | awk 'NR==FNR { line[$1] = $0; next } $1 in line { print line[$1] }' - test.txt | cut -f2 > test_DepthofCoverage_CDS1.txt

file2="/scratch/jc33471/canine_tumor/phylogenetics/vcf/project/szlab/Kun_Lin/Pan_Cancer/ValidationData/MC/store/DepthOfCoverage/009/SRR9911352_DepthofCoverage_CDS.bed"
grep -f test.txt ${file2} | awk 'NR==FNR { line[$1] = $0; next } $1 in line { print line[$1] }' - test.txt | cut -f2  > test_DepthofCoverage_CDS2.txt
paste -d, test_DepthofCoverage_CDS1.txt test_DepthofCoverage_CDS2.txt
> germline_depth_matrix_for_phylo.txt
if [[ ! -e germline_depth_matrix_for_phylo.txt ]]; then
    cat test_DepthofCoverage_CDS1.txt > germline_depth_matrix_for_phylo.txt
fi  
paste -d, germline_depth_matrix_for_phylo.txt test_DepthofCoverage_CDS2.txt > germline_depth_matrix_for_phylo.txt

# java -Xmx60g -cp ${breed_dir}/Pancancer.jar util.GetVariantInfoFromVCF test.txt 5 tmp_vcf_file_list.txt germline_depth_matrix_for_phylo.txt depths
# mask vcf
python /home/jc33471/canine_tumor_wes/scripts/phylogenetic/mask_vcf.py \
    germline_depth_matrix_for_phylo.txt \
    breedSample_breedSpecific_merged_germline_variants_SNP.vcf.gz \
    test_out.vcf.gz \
    10 \
    10

grep -f all_variant_list.txt /scratch/jc33471/canine_tumor/phylogenetics/vcf/project/szlab/Kun_Lin/Pan_Cancer/ValidationData/MC/store/DepthOfCoverage/004/SRR9911377_DepthofCoverage_CDS.bed | awk 'NR==FNR {{ line[$1] = $0; next }} $1 in line {{ print line[$1] }}' - all_variant_list.txt

sed '1d' /scratch/jc33471/canine_tumor/phylogenetics/vcf/project/szlab/Kun_Lin/Pan_Cancer/ValidationData/MC/store/DepthOfCoverage/007/SRR9911361_DepthofCoverage_CDS.bed | cut -f1 > test1.txt
sed '1d' /scratch/jc33471/canine_tumor/phylogenetics/vcf/project/szlab/Kun_Lin/Pan_Cancer/ValidationData/MC/store/DepthOfCoverage/004/SRR9911377_DepthofCoverage_CDS.bed | cut -f2 > test2.txt
paste -d, test1.txt test2.txt > test.txt
for variant in vcf('chr12:220157-220157'):
    print(variant.end)


#############
awk 'BEGIN{FS=":";OFS="\t"}{print $1,$2}' breedSample_all_depth_position_list.txt > breedSample_all_depth_position_list.tab
bcftools view --threads ${SLURM_NTASKS} -O z \
    -R breedSample_all_depth_position_list.tab \
    -o breedSample_all_merged_germline_variants_SNP_depth.vcf.gz \
    breedSample_all_merged_germline_variants_SNP.vcf.gz

# test nj tree here
conda activate phylogenetics

plink --threads ${SLURM_NTASKS} --memory "60000" --seed 123 --dog \
    --distance gz square '1-ibs' \
    --vcf ${run_dir}/merge_vcf/breedSample_all_merged_germline_variants_SNP_depth.vcf.gz \
    --out ${run_dir}/merge_vcf/breedSample_all_merged_germline_variants_SNP_depth

python ${project_dir}/scripts/phylogenetic/mdist2phydist.py \
    --mdist ${run_dir}/merge_vcf/breedSample_all_merged_germline_variants_SNP_depth.mdist.gz \
    --id ${run_dir}/merge_vcf/breedSample_all_merged_germline_variants_SNP_depth.mdist.id \
    --phydist ${run_dir}/merge_vcf/breedSample_all_merged_germline_variants_SNP_depth.phylip 
# neighbor < ${run_dir}/merge_vcf/all_merged_germline_variants_SNP_matrix.mdist.gz

cd ${run_dir}/merge_vcf
echo "breedSample_all_merged_germline_variants_SNP_depth.phylip" > input
echo "Y" >> input
rm out*
neighbor < input > output

Rscript --vanilla ${project_dir}/scripts/phylogenetic/phylip_nwk_vis.R \
    breedSample_all \
    outtree \
    breedSample_all_merged_germline_variants_SNP_depth.mdist.id \
    /home/jc33471/canine_tumor_wes/metadata/data_collection_old.csv \
    test.figure

zcat breedSample_all_merged_germline_variants_SNP_depth.vcf.gz |\
    awk 'BEGIN{{OFS=":"}}/^[^#]/{{print $1,$2}}' > breedSample_all_var_and_depth_position_list.txt
awk '{print $0","}' breedSample_all_var_and_depth_position_list.txt > breedSample_all_var_and_depth_position_list_comma.txt
# depth matrix based on var and depth sites

grep -f breedSample_all_var_and_depth_position_list_comma.txt breedSample_all_germline_depth_matrix_for_phylo.txt > breedSample_all_germline_depth_matrix_for_phylo_clean.txt

python /home/jc33471/canine_tumor_wes/scripts/phylogenetic/mask_vcf.py \
    breedSample_all_germline_depth_matrix_for_phylo_clean.txt \
    breedSample_all_merged_germline_variants_SNP_depth.vcf.gz \
    breedSample_all_merged_germline_variants_SNP_depth_masked.vcf.gz \
    10 \
    ${SLURM_NTASKS} \
    "60G"

plink --threads ${SLURM_NTASKS} --memory "60000" --seed 123 --dog \
    --distance gz square '1-ibs' \
    --vcf ${run_dir}/merge_vcf/breedSample_all_merged_germline_variants_SNP_depth_masked.vcf.gz \
    --out ${run_dir}/merge_vcf/breedSample_all_merged_germline_variants_SNP_depth_masked

python ${project_dir}/scripts/phylogenetic/mdist2phydist.py \
    --mdist ${run_dir}/merge_vcf/breedSample_all_merged_germline_variants_SNP_depth_masked.mdist.gz \
    --id ${run_dir}/merge_vcf/breedSample_all_merged_germline_variants_SNP_depth_masked.mdist.id \
    --phydist ${run_dir}/merge_vcf/breedSample_all_merged_germline_variants_SNP_depth_masked.phylip 
# neighbor < ${run_dir}/merge_vcf/all_merged_germline_variants_SNP_matrix.mdist.gz

cd ${run_dir}/merge_vcf
echo "breedSample_all_merged_germline_variants_SNP_depth_masked.phylip" > input
echo "Y" >> input
rm out*
neighbor < input > output

Rscript --vanilla ${project_dir}/scripts/phylogenetic/phylip_nwk_vis.R \
    breedSample_all \
    outtree \
    breedSample_all_merged_germline_variants_SNP_depth_masked.mdist.id \
    /home/jc33471/canine_tumor_wes/metadata/data_collection_old.csv \
    test_lowcov.pdf

awk 'BEGIN{OFS=","}{print $1,$2}' /scratch/jc33471/canine_tumor/phylogenetics/vcf/project/szlab/Kun_Lin/Pan_Cancer/ValidationData/MC/store/DepthOfCoverage/004/SRR9911377_DepthofCoverage_CDS.bed | grep -f breedSample_all_var_and_depth_position_list_comma.txt | cut -f1,2 > test1.txt
awk 'BEGIN{OFS=","}{print $1,$2}' /scratch/jc33471/canine_tumor/phylogenetics/vcf/project/szlab/Burair/pancancer/scratch_backup/pancancer/WES/germline/osteosarcoma/results/normal/SAMN03732738/SAMN03732738_DepthofCoverage_CDS.bed | grep -f breedSample_all_var_and_depth_position_list_comma.txt | cut -f1,2 > test2.txt
join -a1 -e- -t "," -j 1 -o 0 2.2 <(sort test1.txt) <(sort test2.txt) | tr : , | sort -t, -k1.4,1V -k2,2n | sed 's/,/:/' | less 
join -a1 -e- -t "," -j 1 -o 0 2.2 <(sort breedSample_all_depth_variant_list.txt) <(sort test2.txt) | tr : , | sort -t, -k1.4,1V -k2,2n | sed 's/,/:/' > test2_realigned.txt

join -j 1 -t "," test1.txt test2_realigned.txt | less
sample_name=`basename /scratch/jc33471/canine_tumor_0908/results/Germline/PRJDB10211/HS5/DRR345024_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf "_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf"`
        echo "${sample_name}" > /scratch/jc33471/canine_tumor_0908/results/Germline/PRJDB10211/HS5/DRR345024_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcfrename_samples_vcf.txt
        bcftools reheader -s /scratch/jc33471/canine_tumor_0908/results/Germline/PRJDB10211/HS5/DRR345024_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcfrename_samples_vcf.txt --threads 1 /scratch/jc33471/canine_tumor_0908/results/Germline/PRJDB10211/HS5/DRR345024_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf | awk '{ if($1 ~ /##/ || NF==10)print}' > /scratch/jc33471/canine_tumor_0908/results/Germline/PRJDB10211/HS5/DRR345024_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf.reheader
        bgzip --force --threads 1 /scratch/jc33471/canine_tumor_0908/results/Germline/PRJDB10211/HS5/DRR345024_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf.reheader
        bcftools index --threads 1 /scratch/jc33471/canine_tumor_0908/results/Germline/PRJDB10211/HS5/DRR345024_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf.reheader.gz
        echo "Compressed and indexed /scratch/jc33471/canine_tumor_0908/results/Germline/PRJDB10211/HS5/DRR345024_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf"
 
ll /scratch/jc33471/canine_tumor_0908/results/Germline/PRJDB16014/Beagle1/DRR483792_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf.reheader.gz