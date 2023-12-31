#!/bin/bash
#SBATCH --job-name=004_WES         # Job name (004_WES)
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=4           # Number of CPU cores per task (4)
#SBATCH --mem=60G                   # Job memory limit (10 GB)
#SBATCH --time=100:00:00              # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=004_WES.%j.out    # Standard output log
#SBATCH --error=004_WES.%j.err     # Standard error log


################################################################################################################################
###### complete pipeline (picard, GATK, annova, Mutect, Mutect2,Strelka, Germline mutation, Depth of coverage, Callablebases) ##########
##### WGS analysis needs to change the interval region #####


Bioproject='PRJNA552905'
Normal_Run='SRR9911377'
Tumor_Run='SRR9911376'
CaseName='004'
results='/scratch/kh31516/results/'${Bioproject}'/'${CaseName}'' #
reference='/work/szlab/Lab_shared_PanCancer/source' 
dataN='/scratch/kh31516/data/'${Bioproject}'/'${CaseName}'/'${Normal_Run}
dataT='/scratch/kh31516/data/'${Bioproject}'/'${CaseName}'/'${Tumor_Run}
annovar_index='/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf'
script='/work/szlab/kh31516_Lab_Share_script'
MuTect2_source='/work/szlab/Lab_shared_PanCancer/source/Mutect2'
MuTect_out='/scratch/kh31516/store/'${Bioproject}'/Mutect/'${CaseName}''
MuTect2_out='/scratch/kh31516/store/'${Bioproject}'/Mutect2/'${CaseName}''
Germline_out='/scratch/kh31516/store/'${Bioproject}'/Germline/'${CaseName}''
DepthOfCoverage='/scratch/kh31516/store/'${Bioproject}'/DepthOfCoverage/'${CaseName}''
strelka_out='/scratch/kh31516/store/'${Bioproject}'/strelka/'${CaseName}''
###sraRunTable='/scratch/kh31516/Pan_cancer/PRJNA247493/source/PRJNA247493_SraRunTable.txt'# path


####### Download #######
module load SRA-Toolkit/2.11.1-centos_linux64

mkdir -p $dataT
mkdir -p $dataN
mkdir -p ${results}
mkdir -p ${MuTect_out}
mkdir -p ${MuTect2_out}
mkdir -p ${Germline_out}
mkdir -p ${DepthOfCoverage}
mkdir -p $strelka_out

cd $dataN

fasterq-dump --split-files ${Normal_Run}
#/scratch/kh31516/Pan_cancer/PRJNA247493/download_sra_sample.sh ${Normal_Run} $sraRunTable
gzip -f *

cd $dataT

fasterq-dump --split-files ${Tumor_Run}
#/scratch/kh31516/Pan_cancer/PRJNA247493/download_sra_sample.sh ${Tumor_Run} $sraRunTable
gzip -f *

####### BWA Mapping #######
cd ${results}
ml BWA/0.7.17-GCC-8.3.0

# Tumor
time bwa aln -t 4 ${reference}/canFam3.fa $dataT/${Tumor_Run}_1.fastq.gz > ${results}/${Tumor_Run}_1.sai
time bwa aln -t 4 ${reference}/canFam3.fa $dataT/${Tumor_Run}_2.fastq.gz > ${results}/${Tumor_Run}_2.sai

time bwa sampe ${reference}/canFam3.fa \
${results}/${Tumor_Run}_1.sai \
${results}/${Tumor_Run}_2.sai \
$dataT/${Tumor_Run}_1.fastq.gz \
$dataT/${Tumor_Run}_2.fastq.gz \
> ${results}/${Tumor_Run}.sam

# Normal
time bwa aln -t 4 ${reference}/canFam3.fa $dataN/${Normal_Run}_1.fastq.gz > ${results}/${Normal_Run}_1.sai
time bwa aln -t 4 ${reference}/canFam3.fa $dataN/${Normal_Run}_2.fastq.gz > ${results}/${Normal_Run}_2.sai

time bwa sampe ${reference}/canFam3.fa \
${results}/${Normal_Run}_1.sai \
${results}/${Normal_Run}_2.sai \
$dataN/${Normal_Run}_1.fastq.gz \
$dataN/${Normal_Run}_2.fastq.gz \
> ${results}/${Normal_Run}.sam


####### Convert sam to bam file #######

cd ${results}
module load SAMtools/1.9-GCC-8.3.0

samtools view -bS ${results}/${Tumor_Run}.sam > ${results}/${Tumor_Run}.bam
samtools view -bS ${results}/${Normal_Run}.sam > ${results}/${Normal_Run}.bam

# get header information
#grep '@SQ\|@PG' ${results}/${Normal_Run}.sam > ${results}/header_N
#grep '@SQ\|@PG' ${results}/${Tumor_Run}.sam > ${results}/header_T
samtools view -H ${results}/${Normal_Run}.bam > ${results}/header_N
samtools view -H ${results}/${Tumor_Run}.bam > ${results}/header_T

# exclude unmapped reads based on FLAG
cat ${results}/${Normal_Run}.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > ${results}/${Normal_Run}-cleaned.sam
cat ${results}/header_N ${results}/${Normal_Run}-cleaned.sam > ${results}/fooN
mv ${results}/fooN ${results}/${Normal_Run}-cleaned.sam

cat ${results}/${Tumor_Run}.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > ${results}/${Tumor_Run}-cleaned.sam
cat ${results}/header_T ${results}/${Tumor_Run}-cleaned.sam > ${results}/fooT
mv ${results}/fooT ${results}/${Tumor_Run}-cleaned.sam

####### Picard #######
# picard sort
cd ${results}
module load picard/2.21.6-Java-11

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${results}/${Normal_Run}-cleaned.sam \
O=${results}/${Normal_Run}_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${Normal_Run}

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${results}/${Tumor_Run}-cleaned.sam \
O=${results}/${Tumor_Run}_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${Tumor_Run}

# picard MarkDuplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${results}/${Normal_Run}_rg_added_sorted.bam \
O=${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=${results}/${Normal_Run}-output.metrics REMOVE_SEQUENCING_DUPLICATES=true

java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${results}/${Tumor_Run}_rg_added_sorted.bam \
O=${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=${results}/${Tumor_Run}-output.metrics REMOVE_SEQUENCING_DUPLICATES=true


####### GATK Realign #######
# Generating interval file for sort.bam
ml GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam \
-nt 4 \
-o ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam \
-nt 4 \
-o ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

#### IndelRealigner realign ######
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam \
-targetIntervals ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam.intervals \
-o ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam \
-targetIntervals ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam.intervals \
-o ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

######## Mutect #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 

cd ${MuTect_out}
module load MuTect/1.1.7-Java-1.7.0_80

time java -jar $EBROOTMUTECT/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence ${reference}/canFam3.fa \
--dbsnp ${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
--defaultBaseQualities 30 \
--intervals ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
--input_file:normal ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--input_file:tumor ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--out ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.MuTect.vcf

####### Annovar #######
# Extract PASS records from vcf
module load Perl/5.26.1-GCCcore-6.4.0

awk '$7 == "PASS" {print $0}' ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.MuTect.vcf > ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $annovar_index

# add gene name
python $script/Add_GeneName_N_Signature.py ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


### Sanger 5 steps filtering ###
# 5 Steps filtering
grep -w KEEP ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > ${MuTect_out}/${CaseName}_PASS.stat

python $script/Filter_MutectStat_5steps.py \
${MuTect_out}/${CaseName}_PASS.stat \
${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
${MuTect_out}/${CaseName}_vaf_before.txt \
${MuTect_out}/${CaseName}_vaf_after.txt \
${MuTect_out}/${CaseName}_whyout.txt \
$CaseName 
# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut > ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput $annovar_index

# add gene name
python $script/Add_GeneName_N_Signature.py ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


############################## Germline mutation preparation Start ######################################
# GATK 3 
ml GATK/3.8-1-Java-1.8.0_144

cd $Germline_out

# Variant calling
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference}/canFam3.fa -I \
${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases \
-stand_call_conf 20.0 -o ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference}/canFam3.fa -I \
${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases \
-stand_call_conf 20.0 -o ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.vcf

# Variant filtering
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${reference}/canFam3.fa \
-V ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.vcf \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${reference}/canFam3.fa \
-V ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.vcf \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf

################ Annovar ###################

cd ${Germline_out}
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf \
> ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS

awk '$7 == "PASS" {print $0}' ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf \
> ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS

# annovar input preparation
module load Perl/5.26.1-GCCcore-6.4.0

perl $annovar_index/convert2annovar.pl -format vcf4old ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS \
> ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput

perl $annovar_index/convert2annovar.pl -format vcf4old ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS \
> ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput $annovar_index

perl $annovar_index/annotate_variation.pl --buildver canFam3 ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput $annovar_index

# add gene names to annovar output

cd ${Germline_out}
python $script/Add_GeneName_N_Signature.py ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function \
${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt

python $script/Add_GeneName_N_Signature.py ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function \
${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


#### GATK callable ####
cd ${Germline_out}
ml GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
-R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-summary ${Germline_out}/${Normal_Run}_table.txt \
-o ${Germline_out}/${Normal_Run}_callable_status.bed 

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
-R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-summary ${Germline_out}/${Tumor_Run}_table.txt \
-o ${Germline_out}/${Tumor_Run}_callable_status.bed 

############################## Germline mutation preparation End ######################################


#### GATK DepthofCoverage ####
cd ${DepthOfCoverage}
ml GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
-o ${DepthOfCoverage}/${Normal_Run}_DepthofCoverage_CDS.bed 

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
-o ${DepthOfCoverage}/${Tumor_Run}_DepthofCoverage_CDS.bed 


####### Mutect2 #######

cd ${MuTect2_out}
module load SAMtools/1.9-GCC-8.3.0
module load GATK/4.1.6.0-GCCcore-8.3.0-Java-1.8

Normal_sample=$(samtools view -H ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)
Tumor_sample=$(samtools view -H ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)

## Mutect2 follows the Glioma PanCancer parameters
gatk Mutect2 \
-R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-normal $Normal_sample \
--callable-depth 8 \
--dont-use-soft-clipped-bases true \
--panel-of-normals  $MuTect2_source/pon.vcf.gz \
--initial-tumor-lod 2.0 --normal-lod 2.2 --tumor-lod-to-emit 3.0 --pcr-indel-model CONSERVATIVE \
--f1r2-tar-gz ${MuTect2_out}/${CaseName}-f1r2.tar.gz \
-L ${MuTect2_source}/Uniq-Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval.list \
-O ${MuTect2_out}/${CaseName}_MuTect2_GATK4_noDBSNP.vcf


## pass this raw data to LearnReadOrientationModel:
gatk LearnReadOrientationModel -I ${MuTect2_out}/${CaseName}-f1r2.tar.gz -O ${MuTect2_out}/read-orientation-model.tar.gz

# filtered Mutect2 files

gatk FilterMutectCalls -R ${reference}/canFam3.fa \
-ob-priors ${MuTect2_out}/read-orientation-model.tar.gz \
-V ${MuTect2_out}/${CaseName}_MuTect2_GATK4_noDBSNP.vcf \
-O ${MuTect2_out}/filtered-${CaseName}_MuTect2_GATK4_noDBSNP.vcf

cd ${MuTect2_out}

### filter Mutect2 #####

## with DbSNP
java -cp  $script/ DbSNP_filtering \
${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${MuTect2_out}/filtered-${CaseName}_MuTect2_GATK4_noDBSNP.vcf \
${MuTect2_out}/DbSNP_filtered-${CaseName}_MuTect2_GATK4.vcf


## with PON
java -cp  $script/ DbSNP_filtering \
${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${MuTect2_out}/DbSNP_filtered-${CaseName}_MuTect2_GATK4.vcf \
${MuTect2_out}/PON_DbSNP_filtered-${CaseName}_MuTect2_GATK4.vcf


### Mutect2 5 steps filtering

cd ${MuTect2_out}

ml Anaconda3/2020.02
#source activate py35

#awk '$7 == "PASS" {print $0}' ${MuTect2_out}/PON_DbSNP_filtered-${CaseName}_MuTect2_GATK4.vcf > ${MuTect2_out}/PON_DbSNP_filtered-${CaseName}_MuTect2_GATK4.vcf-PASS

# 5 Steps filtering
python $script/Mutect2_5Steps_filtering.py \
${MuTect2_out}/PON_DbSNP_filtered-${CaseName}_MuTect2_GATK4.vcf \
${MuTect2_out}/PON_DbSNP_filtered-${CaseName}_MuTect2_GATK4.vcf \
${MuTect2_out}/${CaseName}_VAF_Before.txt \
${MuTect2_out}/${CaseName}_VAF_After.txt \
${CaseName}


## Strelka
cd ${results}
module load SAMtools/1.9-GCC-8.3.0
#samtools index ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam
#samtools index ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam

${script}/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
--normalBam ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--tumorBam ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--referenceFasta ${reference}/canFam3.fa \
--runDir $strelka_out

# your result will be generated in demo_somatic/results/variants
$strelka_out/runWorkflow.py -m local -j 20


## limit strelka result into CDS region

python2 $script/Limit_vcf_to_CDS.py $strelka_out/results/variants/somatic.indels.vcf.gz $reference/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list


### Annovar for strelka
cd $strelka_out/results/

####### Annovar #######
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS > $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS

# annovar input preparation
perl $reference/annovar_CanFam3.1.99.gtf/convert2annovar.pl -format vcf4old $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS > \
$strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS-avinput

# annovar annotate
perl $reference/annovar_CanFam3.1.99.gtf/annotate_variation.pl --buildver canFam3 $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS-avinput \
$reference/annovar_CanFam3.1.99.gtf

# add gene name
python2 $script/Add_GeneName_N_Signature.py $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS-avinput.exonic_variant_function \
$reference/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


: "
####### remove unneeded files #######
rm ${results}/*_sorted_dedupped_removed.bam
rm ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS ${MuTect_out}/${CaseName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
rm ${results}/*.bai
rm ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
rm ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
rm ${results}/${Normal_Run}.sam ${results}/${Tumor_Run}.sam 
rm ${results}/header_N ${results}/header_T
rm ${results}/${Normal_Run}_rg_added_sorted.bam ${results}/${Tumor_Run}_rg_added_sorted.bam
rm ${results}/${Normal_Run}-cleaned.sam ${results}/${Tumor_Run}-cleaned.sam
rm $dataN/*.fastq.gz $dataT/*.fastq.gz
rm ${results}/*.sai
rm ${MuTect2_out}/DbSNP_filtered-${CaseName}_MuTect2_GATK4.vcf
rm ${MuTect2_out}/read-orientation-model.tar.gz
rm ${MuTect2_out}/${CaseName}-f1r2.tar.gz
"
