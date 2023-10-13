#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=master_job
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --tasks-per-node=8
#SBATCH --mem=60G
#SBATCH --time=500:00:00
#SBATCH --mail-user=jc33471@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/scratch/jc33471/canine_tumor_0908/master_job.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate wes_env

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor_0908"
mkdir -p ${run_dir}
cd ${run_dir}
mkdir -p logs/ config/ data/ out/ results/

# config="/home/jc33471/canine_tumor_wes/scripts/sum_cases/config.json"
config=$run_dir/config/master_config.json

python ${project_dir}/scripts/sum_cases/make_snakemake_config.py \
    --project_dir ${project_dir} \
    --out ${config} \
    --outdir ${run_dir} \
    --metadata ${project_dir}/metadata/data_collection_new.csv\
    --threads 8 \
    --memory "60G"

snakemake \
    --jobs 200 \
    --use-conda \
    --latency-wait 60 \
    --keep-going \
    --restart-times 3 \
    --rerun-triggers mtime \
    --configfile ${config} \
    --snakefile "${project_dir}/scripts/sum_cases/Snakefile" \
    --cluster-status "${project_dir}/scripts/status-sacct-robust.sh" \
    --cluster-cancel 'scancel' \
    --cluster '
        sbatch \
            --partition=batch \
            --job-name={wildcards.Bioproject}_{wildcards.CaseName}_smk \
            --nodes=1 \
            --ntasks={threads} \
            --tasks-per-node={threads} \
            --mem={resources.mem} \
            --time=72:00:00 \
            --parsable \
            --mail-user=jc33471@uga.edu \
            --mail-type=FAIL \
            --output=logs/{wildcards.Bioproject}_{wildcards.CaseName}_log.o \
            --error=logs/{wildcards.Bioproject}_{wildcards.CaseName}_log.e'

## combine QC results

for group in `cut -d, -f4,9 /home/jc33471/canine_tumor_wes/metadata/data_collection_new.csv | grep -v "DiseaseAcronym2" | sort -u | sed 's/\"//g'`; do
    Cancer_Type=`echo ${group} | cut -d, -f1`
    Bioproject=`echo ${group} | cut -d, -f2`    
    cat ${run_dir}/results/QC/*/*/*${Cancer_Type}_${Bioproject}-CDS_Mapping_summary.txt > ${run_dir}/results/QC/Total_WES_BWA_CDS_${Bioproject}_${Cancer_Type}.txt
    cat ${run_dir}/results/QC/*/*/*${Cancer_Type}_${Bioproject}-mapping_quality_line.txt > ${run_dir}/results/QC/WES_Total_Mapping_quality_${Bioproject}_${Cancer_Type}.txt
    sed '1d' ${run_dir}/results/QC/*/*/*${Cancer_Type}_${Bioproject}_randomness_summary.txt > ${run_dir}/results/QC/Total_WES_${Bioproject}_${Cancer_Type}_Randomness_Summary.txt
done

cat ${run_dir}/results/QC/Total_WES_BWA_CDS_*.txt | grep -v "^File_name" |\
    sed '1 i\File_name\tTotal_Pairs\tUniq_mapped_rate\tUniq_Exonic_region_mapped_rate\tCancerType\tStatus\tUnmappedRate\tDuplicateMapped_rate\tOnemapped_rate\tIncorrectMapped_rate\tTotal_line\tTotal_unique\tTotal_pass\tTotal_Unmapped\tTotal_Duplicate\tTotal_Onemapped\tTotal_Incorrect' \
    > ${run_dir}/results/QC/Total_WES_BWA_CDS.combined.txt
cat ${run_dir}/results/QC/WES_Total_Mapping_quality_*.txt | grep -v "^File_name" |\
    sed '1 i\File_name\tfra30\tfra60\tCancer_Type\tStatus' > ${run_dir}/results/QC/WES_Total_Mapping_quality.combined.txt
cat ${run_dir}/results/QC/Total_WES_*_Randomness_Summary.txt | grep -v "^file_name" |\
    sed '1 i\File_name\taverage\trmse\tsumOfSqerror\tCancer_type\tStatus' > ${run_dir}/results/QC/Total_WES_Randomness_Summary.combined.txt

