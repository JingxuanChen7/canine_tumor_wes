#!/bin/bash
#SBATCH --job-name=CMT-2_Summarized_BWA_WES         # Job name (Summarized_BWA_WES)
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=1           # Number of CPU cores per task (4)
#SBATCH --mem=55G                   # Job memory limit (10 GB)
#SBATCH --time=10:00:00              # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=CMT-2_Summarized_BWA_WES.%j.out    # Standard output log
#SBATCH --error=CMT-2_Summarized_BWA_WES.%j.err     # Standard error log


# WGS 100 bps 
# WES 76 bps

Sequence_length=101
sample_name='CMT-2'
source_folder='/scratch/kh31516/Original_Mammary/results/PRJNA489159/'${sample_name} #
script='/home/kh31516/kh31516_Lab_Share_script'
sam_file_output='/scratch/kh31516/Original_Mammary/results/PRJNA489159/QC/'${sample_name}
summary_output='/scratch/kh31516/Original_Mammary/results/PRJNA489159/'
DepthOfCoverage='/scratch/kh31516/Original_Mammary/results/PRJNA489159/DepthOfCoverage/'${sample_name}
Mutect='/scratch/kh31516/Original_Mammary/results/PRJNA489159/Mutect/'${sample_name}
Cancer_Type='MT'
Normal_Sample='SRR7780922'
Tumor_Sample='SRR7780923'
Bioproject='PRJNA489159'

module load SAMtools/1.9-GCC-8.3.0
ml Anaconda3/2020.02
#source activate py35

mkdir -p ${sam_file_output}
## Convert the bam file into sam file
cd ${source_folder}/

samtools view $source_folder/${Normal_Sample}.bam > ${sam_file_output}/${Normal_Sample}.sam # normal
samtools view $source_folder/${Tumor_Sample}.bam > ${sam_file_output}/${Tumor_Sample}.sam # tumor

###### Analyze the BWA mapping and CDS region mapping result ######

cd ${sam_file_output}/

java -cp ${script}/ Line_by_line_Total_dict_Get_exon_reads \
${sam_file_output}/${Normal_Sample}.sam \
${Normal_Sample} \
${sam_file_output}/Normal-${Normal_Sample}-CDS_Mapping_summary.txt \
$Sequence_length \
${Cancer_Type} \
"Normal"

java -cp ${script}/ Line_by_line_Total_dict_Get_exon_reads \
${sam_file_output}/${Tumor_Sample}.sam \
${Tumor_Sample} \
${sam_file_output}/Tumor-${Tumor_Sample}-CDS_Mapping_summary.txt \
$Sequence_length \
${Cancer_Type} \
"Tumor"

cat ${sam_file_output}/Normal-${Normal_Sample}-CDS_Mapping_summary.txt >> $summary_output/Total_WES_BWA_CDS_${Bioproject}_${Cancer_Type}.txt
cat ${sam_file_output}/Tumor-${Tumor_Sample}-CDS_Mapping_summary.txt >>  $summary_output/Total_WES_BWA_CDS_${Bioproject}_${Cancer_Type}.txt

##### Sequence Reads Mapping Quality #####

cd ${sam_file_output}/

cat ${sam_file_output}/${Normal_Sample}.sam | cut -f5 > ${sam_file_output}/${Normal_Sample}-mapping_quality

Ngt30=$(cat ${Normal_Sample}-mapping_quality | awk '{if ($1>=30) {print $1}}' | wc -l)
Ngt60=$(cat ${Normal_Sample}-mapping_quality | awk '{if ($1>=60) {print $1}}' | wc -l)
Ntotal=$(cat ${Normal_Sample}-mapping_quality| wc -l)
Nfra30=$(echo "$((Ngt30))/$((Ntotal))" | bc -l)
Nfra60=$(echo "$((Ngt60))/$((Ntotal))" | bc -l)

printf  "%s\t%4f\t%4f\t%s\t%s\n" "${Normal_Sample}" "${Nfra30}" "${Nfra60}" "${Cancer_Type}" "Normal" >> \
$summary_output/WES_Total_Mapping_quality_${Bioproject}_${Cancer_Type}.txt

cat ${sam_file_output}/${Tumor_Sample}.sam | cut -f5 > ${sam_file_output}/${Tumor_Sample}-mapping_quality

Tgt30=$(cat ${sam_file_output}/${Tumor_Sample}-mapping_quality | awk '{if ($1>=30) {print $1}}' | wc -l)
Tgt60=$(cat ${sam_file_output}/${Tumor_Sample}-mapping_quality | awk '{if ($1>=60) {print $1}}' | wc -l)
Ttotal=$(cat ${sam_file_output}/${Tumor_Sample}-mapping_quality| wc -l)
Tfra30=$(echo "$((Tgt30))/$((Ttotal))" | bc -l)
Tfra60=$(echo "$((Tgt60))/$((Ttotal))" | bc -l)

printf   "%s\t%4f\t%4f\t%s\t%s\n" "${Tumor_Sample}" "${Tfra30}" "${Tfra60}" "${Cancer_Type}" "Tumor" >> \
$summary_output/WES_Total_Mapping_quality_${Bioproject}_${Cancer_Type}.txt

###### QC of Depth of Coverage and Randomness ######

cd ${DepthOfCoverage}/

cat ${DepthOfCoverage}/${Normal_Sample}_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > \
${DepthOfCoverage}/${Normal_Sample}_DepthofCoverage_Distribution.txt
cat ${DepthOfCoverage}/${Tumor_Sample}_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > \
${DepthOfCoverage}/${Tumor_Sample}_DepthofCoverage_Distribution.txt

#cat ${Tumor_Sample}_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
#cat ${Normal_Sample}_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
python ${script}/Randomness.py ${DepthOfCoverage}/${Normal_Sample}_DepthofCoverage_Distribution.txt ${Normal_Sample} ${Cancer_Type} Normal
python ${script}/Randomness.py ${DepthOfCoverage}/${Tumor_Sample}_DepthofCoverage_Distribution.txt ${Tumor_Sample} ${Cancer_Type} Tumor

cat ${DepthOfCoverage}/${Normal_Sample}_randomness_summary.txt >> $summary_output/Total_WES_${Cancer_Type}_${Bioproject}_Randomness_Summary.txt
cat ${DepthOfCoverage}/${Tumor_Sample}_randomness_summary.txt >> $summary_output/Total_WES_${Cancer_Type}_${Bioproject}_Randomness_Summary.txt

###### Callable bases ########


reference='/work/szlab/Lab_shared_PanCancer/source'
#Yuan_script='/work/szlab/Lab_shared_PanCancer/script'

cd ${sam_file_output}

callable=$(cat $Mutect/${sample_name}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt | grep "^1"| wc -l)

printf  "%s\t%d\t%s\t%\%s\n" "${sample_name}" "${callable}" "${Cancer_Type}" "${Bioproject}" >> $summary_output/Total_WES_Callable_${Bioproject}_${Cancer_Type}.txt

: "
rm ${sam_file_output}/${Normal_Sample}.sam
rm ${sam_file_output}/${Tumor_Sample}.sam
rm ${sam_file_output}/${Normal_Sample}_BWA_summarize
rm ${sam_file_output}/${Tumor_Sample}_BWA_summarize
rm ${sam_file_output}/Normal-${Normal_Sample}-CDS_Mapping_summary.txt
rm ${sam_file_output}/Tumor-${Tumor_Sample}-CDS_Mapping_summary.txt
rm ${sam_file_output}/${Normal_Sample}-mapping_quality ${sam_file_output}/${Tumor_Sample}-mapping_quality
rm ${DepthOfCoverage}/${Normal_Sample}_DepthofCoverage_Distribution.txt
rm ${DepthOfCoverage}/${Tumor_Sample}_DepthofCoverage_Distribution.txt
rm ${DepthOfCoverage}/${Tumor_Sample}_randomness_summary.txt
rm ${DepthOfCoverage}/${Normal_Sample}_randomness_summary.txt
"
