#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=wgs_master
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --tasks-per-node=8
#SBATCH --mem=60G
#SBATCH --time=500:00:00
#SBATCH --mail-user=jc33471@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/scratch/jc33471/canine_tumor/wgs_breed_prediction/wgs_master.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate wes_env

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor/wgs_breed_prediction"
breed_dir="/home/${USER}/breed_prediction"
mkdir -p ${run_dir}
cd ${run_dir}

mkdir -p vcf/ merge_vcf/ logs/

# cd $run_dir/vcf

config="/home/jc33471/canine_tumor_wes/scripts/wgs/config.json"

snakemake \
    -np \
    --jobs 200 \
    --use-conda \
    --latency-wait 60 \
    --keep-going \
    --rerun-incomplete \
    --snakefile "${project_dir}/scripts/wgs/Snakefile" \
    --configfile ${config} \
    --rerun-triggers mtime \
    --cluster-cancel 'scancel' \
    --cluster '
        sbatch \
            --partition=batch \
            --nodes=1 \
            --ntasks={threads} \
            --tasks-per-node={threads} \
            --mem={resources.mem} \
            --time=72:00:00 \
            --parsable \
            --mail-user=jc33471@uga.edu \
            --mail-type=FAIL \
            --output=logs/slurm-%j.o \
            --error=logs/slurm-%j.e'



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

bcftools index $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing.vcf.gz

# # remove useless contig names in header
# bcftools view -h $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_chr1region.vcf.gz |\
#     grep -v -E "chrUn|chrY|chrM" > $run_dir/vcf/header_modified.txt
# bcftools reheader --header $run_dir/vcf/header_modified.txt \
#     -o $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_chr1region_reheader.vcf.gz \
#     $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_chr1region.vcf.gz 

# conda activate wes_env
# # wget -O ${run_dir}/vcf/canFam4ToCanFam3.over.chain.gz "http://hgdownload.soe.ucsc.edu/goldenPath/canFam4/liftOver/canFam4ToCanFam3.over.chain.gz"
# picard -Xmx60G LiftoverVcf \
#     I=${run_dir}/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_chr1region_reheader.vcf.gz \
#     O=${run_dir}/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_chr1region_liftover.vcf.gz \
#     CHAIN=${run_dir}/vcf/canFam4ToCanFam3.over.chain.gz \
#     REJECT=${run_dir}/vcf/rejected_variants.vcf.gz \
#     WARN_ON_MISSING_CONTIG=true \
#     R="/work/szlab/Lab_shared_PanCancer/source/canFam3.fa"
# zcat Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_chr1region_liftover.vcf.gz | awk '{if($0 !~ /^#/) print $1"\t"$2}' > test_pos.txt
# grep -f test_pos.txt /work/szlab/Lab_shared_PanCancer/source/DbSNP_canFam3_version151-DogSD_Broad_March2022.vcf > overlap.vcf


# canfam4->canfam3 liftover
# wget -O ${run_dir}/vcf/canFam4ToCanFam3.over.chain.gz "http://hgdownload.soe.ucsc.edu/goldenPath/canFam4/liftOver/canFam4ToCanFam3.over.chain.gz"
grep ">" /work/szlab/Lab_shared_PanCancer/source/canFam3.fa | sed 's/>//g'  > $run_dir/vcf/canfam3_chr.txt
for chr in `cat $run_dir/vcf/canfam3_chr.txt`; do

    echo "Processing $chr ..."
    bcftools view --threads $SLURM_NTASKS \
        -r ${chr} --output-type z \
        -o $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_${chr}.vcf.gz \
        $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing.vcf.gz
    
    picard -Xmx60G LiftoverVcf \
        I=${run_dir}/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_${chr}.vcf.gz \
        O=${run_dir}/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_liftover_${chr}.vcf.gz \
        CHAIN=${run_dir}/vcf/canFam4ToCanFam3.over.chain.gz \
        REJECT=${run_dir}/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_liftover_${chr}_rejected_variants.vcf.gz \
        WARN_ON_MISSING_CONTIG=true \
        R="/work/szlab/Lab_shared_PanCancer/source/canFam3.fa"

done


bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%VAF]\n' -o $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_liftover.vaf.txt \
    $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_liftover.vcf.gz

bcftools view --threads ${SLURM_NTASKS} \
    -S $run_dir/merge_vcf/breedSample.list --force-samples \
    $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_liftover_chr1:1-10000001.vcf.gz |\
bcftools +fill-tags - -O z --threads $SLURM_NTASKS -o $run_dir/merge_vcf/test.vcf.gz -- -t FORMAT/VAF 

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%VAF]\n' -o $run_dir/merge_vcf/test.vaf_matrix.txt

bcftools query -l $run_dir/merge_vcf/test.vcf.gz | sed '1i Chromosome\nPosition\nRef\nAlt' | awk 'BEGIN{ORS="\t"}{print}' | sed 's/\t$/\n/' > $run_dir/merge_vcf/test_header.txt
cat $run_dir/merge_vcf/test_header.txt $run_dir/merge_vcf/test.vaf_matrix.txt > $run_dir/merge_vcf/test.vaf_matrix.header.txt
awk 'BEGIN{FS=",";OFS="\t"}{print $1,$3,$4,$5,$6,$7}' "/home/jc33471/canine_tumor_wes/metadata/WGS/Dog10K_breeds.csv" > $run_dir/merge_vcf/breed_prediction_metadata.txt

time python /home/jc33471/canine_tumor_wes/scripts/wgs/vaf_matrix.py \
    test.vcf.gz \
    py_vaf_matrix.header.txt.gz \
    8 \
    "60G"

awk '{if($1=="None"){split($2,locus,":"); split($3,mut,">"); print locus[1]"\t"locus[2]"\t"mut[1]"\t"mut[2]}}' test_specific.txt > test_specific_reformat.txt
zgrep -f test_specific_reformat.txt py_vaf_matrix.header.txt.gz >  py_vaf_matrix.header.breed.txt

Rscript --vanilla /home/jc33471/canine_tumor_wes/scripts/wgs/breed_specific_variants_wgs.R             /home/jc33471/canine_tumor_wes/scripts/breed_prediction/build_sample_meta_data.R             /scratch/jc33471/canine_tumor/wgs_breed_prediction/vaf_matrix/chr34:40000001-42124431.vaf_matrix.txt.gz             /scratch/jc33471/canine_tumor/wgs_breed_prediction/breed_variants/chr34:40000001-42124431/breed_unique_variants_chr34:40000001-42124431.txt             /scratch/jc33471/canine_tumor/wgs_breed_prediction/breed_variants/chr34:40000001-42124431/breed_enriched_variants_chr34:40000001-42124431.txt             /scratch/jc33471/canine_tumor/wgs_breed_prediction/breed_variants/chr34:40000001-42124431/breed_specific_variants_chr34:40000001-42124431.txt             /scratch/jc33471/canine_tumor/wgs_breed_prediction/breed_variants/breed_prediction_metadata.txt             params.breedlist &> /scratch/jc33471/canine_tumor/wgs_breed_prediction/logs/breed_specific_idenfication_chr34:40000001-42124431.log

bcftools concat --threads ${SLURM_NTASKS} --allow-overlaps \
    -o $run_dir/merge_vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_liftoverConcat.vcf.gz -O z \
    $(ls $run_dir/vcf/Dog10K_AutoAndXPAR_SNPs_filtered_added_fmissing_liftover_*.vcf.gz | grep -v "rejected_variants" | sort -V)

grep -v "NA>NA" breed_specific_variants_concat.txt |\
    sed '1d' |\
    awk 'BEGIN{FS=OFS="\t"}{split($2,locus,":"); split($3,mut,">"); print locus[1],locus[2]-1,locus[2],mut[1],mut[2],$5}' |\
    sort -t $'\t' -V -k 1.4,1  -k 2,2  > sorted_breed_specific_variants_concat.txt

awk -F '[:-]' '{print $1"\t"$2"\t"$3}' /work/szlab/Lab_shared_PanCancer/source/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list |\
    sort -t $'\t' -V -k 1.4,1  -k 2,2 -k 3,3 |\
    uniq > sorted_canfam3_cds.txt

module load BEDTools/2.30.0-GCC-12.2.0
bedtools intersect -a sorted_breed_specific_variants_concat.txt -b sorted_canfam3_cds.txt -sorted |\
    awk 'BEGIN{FS=OFS="\t"}{print "None",$1":"$3,$4">"$5,"None",$6}' 
head -n1 breed_specific_variants_concat.txt

# merge WES with WGS variants
cd /scratch/jc33471/canine_tumor/breed_prediction
mkdir -p output_include_WGS
cat /scratch/jc33471/canine_tumor/breed_prediction/output_exclude_WGS/all_breed_specific_variants.txt \
    <(sort -u /scratch/jc33471/canine_tumor/wgs_breed_prediction/breed_variants/concat/breed_specific_variants_CDS.txt | sed '1d') \
    > /scratch/jc33471/canine_tumor/breed_prediction/output_include_WGS/all_breed_specific_variants.txt

awk -F '[:-]' '{print $1"\t"$2"\t"$3}' "/work/szlab/Lab_shared_PanCancer/source/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list" | sort -t $'\t' -V -k 1.4,1  -k 2,2 -k 3,3 | uniq > /scratch/jc33471/canine_tumor/wgs_breed_prediction/breed_variants/concat/sorted_canfam3_cds_interval.txt
bedtools intersect -a <(zcat chr26:1-10000001.vaf_matrix.txt.gz | sed '1d' | awk 'BEGIN{FS=OFS="\t"}{$2=$2-1 OFS $2}1')  -b /scratch/jc33471/canine_tumor/wgs_breed_prediction/breed_variants/concat/sorted_canfam3_cds_interval.txt  > test.out

zcat chr26:1-10000001.vaf_matrix.txt.gz | sed '1d' | awk 'BEGIN{FS=OFS="\t"}{$2=$2-1 OFS $2}1' | less -S
bedtools intersect -a <( zcat /scratch/jc33471/canine_tumor/wgs_breed_prediction/vaf_matrix/chr16:50000001-59632846.vaf_matrix.txt.gz | awk 'BEGIN{FS=OFS="\t"}{$2=$2-1 OFS $2}1' )             -b /scratch/jc33471/canine_tumor/wgs_breed_prediction/breed_variants/concat/sorted_canfam3_cds_interval.txt |            cut -f2 --complement > /scratch/jc33471/canine_tumor/wgs_breed_prediction/breed_variants/chr16:50000001-59632846/chr16:50000001-59632846.cds.vaf_matrix.txt