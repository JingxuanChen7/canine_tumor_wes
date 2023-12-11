# Installation
- Environment for meta data collection
```bash
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/data_collection.yml --name data_collection
```
- Environments for germline and depthofcoverage pipeline
```bash
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/wes_env.yml --name wes_env
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/annovar_env.yml --name annovar_env
# mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/mutect2_env.yml --name mutect2_env
#mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/strelka_env.yml --name strelka_env
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/java17.yml --name java17
```
- Environment for breed prediction pipeline
```
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/breed_prediction.yml --name breed_prediction
```
- Environment for phylogenetic analysis
```
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/phylogenetics.yml --name phylogenetics
```
# Data collection

Bioproject   ID | Tumor Type | Total runs under Bioproject | Paired cases (matched normal   & tumor samples) | Read length | Year | Submitter
-- | -- | -- | -- | -- | -- | --
PRJEB53653 | Mammary neoplasia (MT) | 111 | 53(106); 5 unpaired | 151 | 2022 | Uppsala Uuniversity
PRJDB10211 | Histiocytic sarcoma (HS) | 10 | 5(10) | 150 | 2022 | The University of Tokyo
PRJNA786469 | oral melanoma (OM) | 145 | 72(144); 1 unpaired | 76 | 2021 | CNRS - University Rennes1
PRJNA752630 | diffuse large B-cell lymphoma   (DLBCL) | 154 | 77(154) | 151 | 2021 | University of Padua
PRJNA701141 | Osteosarcoma (OSA) | 16 | 1(2) | 151 | 2021 | University of Padua
PRJNA695534 | B-Cell Lymphoma (BCL) | 142 | 71(142) | 125 | 2021 | Cornell University
PRJNA680382 | Osteosarcoma (OSA) | 16 | 4(8) | 151 | 2020 | University of Missouri
PRJNA677995 | splenic angiosarcoma   (hemangiosarcoma) (HSA) | 83 | 4(8) | 101 | 2020 | Translational Genomics Research   Institute
PRJDB16014 | inflammatory colorectal polyp | 4 | 4 unpaired | 106 | 2023 | National Institute for Physiological Sciences
PRJNA891496 | Multiple non-tumor diseases | 49 dogs; unpaired | 49 unpaired | 76 | 2022 | Ghent University
PRJNA630029 | fearfulness? | 22 | 16 unpaired | 101 | 2020 | University of Helsinki
PRJNA616374 | Primary bladder tumor | 11 | 11 unpaired | 100x100bp | 2020 | Colorado State University
PRJNA613479 | Osteosarcoma (OSA) | 27 | 27 unpaired | 150x150bp | 2020 | Colorado State University
PRJEB57227 | Urothelial carcinoma (UC) | 829 | 87(174) | 100x100bp | 2022 | Wellcome Sanger Institute

- The above table can also be found at `metadata/dataset_sum.xlsx`.
- More detailed metatable was generated with R script `scripts/data_collection/data_collection.R`.
  - Metatable for my runs can be found at `metadata/data_collection_new.csv`.

# Previous results preparation

- Backup results are placed in `/project`, move to `/scratch` so that they can be assessed through work node. Involoved result files include:
  - `*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf`,
  - `*_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName`,
  - and `*_DepthofCoverage_CDS.bed`.

```bash
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
```
- The `backup_vcf_file.list` for previous run prior to 2021 can also be found at `metadata/backup_vcf_file.list`.

# Run WES pipeline (germline & depthofcoverage) and QC
- Shell script to submit jobs:
  - Snakemake pipeline for each individual job (case): `scripts/per_case/Snakefile`
  - Snakemake QC pipeline for each individual job (case): `scripts/qc/Snakefile`
  - Snakemake pipeline to submit all jobs: `scripts/sum_cases/Snakefile`
    - Note: the `output_report` should be edited after finishing WES analysis pipeline for QC.

```bash
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
    --latency-wait 999 \
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
```

- NOTE: The environment installation is mostly reproducible other than `MuTect`, because `MuTect/1.1.7-Java-1.7.0_80` is pre-installed by UGA GACRC. Double check with `module spider`.
  - This module was not performed in my analysis.

- After QC finished, create dot plots and update metadata to include QC results with r script `scripts/qc/combine_qc_table.R`.
  - Note: One additional step of QC is to exclude outliers in phylogenetic analysis. The list of samples should be manually provided in the r script.

## An example to test individual snakefiles
```bash
#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=test_per_case
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --tasks-per-node=8
#SBATCH --mem=60G
#SBATCH --time=500:00:00
#SBATCH --mail-user=jc33471@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/scratch/jc33471/canine_tumor/test_per_case.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate wes_env

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor"
mkdir -p ${run_dir}
cd ${run_dir}
mkdir -p logs/ config/ data/ out/ results/

#### test WES germline calling pipeline per case
config=$run_dir/config/test.json

python ${project_dir}/scripts/per_case/make_snakemake_config.py \
    --project_dir ${project_dir} \
    --out ${config} \
    --outdir ${run_dir} \
    --Bioproject "PRJEB53653" \
    --Normal_Run "ERR9923018" \
    --Tumor_Run "ERR9923074" \
    --CaseName "24" \
    --threads 8 \
    --memory "60G"

snakemake \
    --cores ${SLURM_NTASKS} \
    --use-conda \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --configfile ${config} \
    --snakefile "${project_dir}/scripts/per_case/Snakefile"


### test QC
config=$run_dir/config/test_qc.json

python ${project_dir}/scripts/qc/make_snakemake_config.py \
    --project_dir ${project_dir} \
    --out ${config} \
    --outdir ${run_dir} \
    --Bioproject "PRJEB53653" \
    --Normal_Run "ERR9923018" \
    --Tumor_Run "ERR9923074" \
    --CaseName "24" \
    --metadata ${project_dir}/metadata/data_collection_new.csv \
    --readlength ${project_dir}/metadata/data_new_readlength.csv \
    --threads 8 \
    --memory "60G"

snakemake \
    --cores ${SLURM_NTASKS} \
    --use-conda \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --configfile ${config} \
    --snakefile ${project_dir}"/scripts/qc/Snakefile"
```

# Breed prediction pipeline
## Combining individual VCF files into MAF matrix
- Complete pipeline can be found here `scripts/breed_prediction/combine_MAF.sh`
- Note: To run this script, one has to make sure individual VCF files are available on `/scratch`. Paths to current files should be carefully checked in section `create input file list` of the above script.

## Heatmap visualization
- Firstly, run Run R script `scripts/breed_prediction/breeds_joint_heatmap_beforeQC.R` interactively in R studio.  (command line interface not implemented yet)
  - This step is to generate heatmap for breed validation/prediction.
  - If running on local machine, one has to download `breed_predection_metadata.txt`, MAF matrix with low coverage sites masked, and breed specific variants output with folder struncture.
- Run R script `scripts/breed_prediction/breeds_joint_heatmap_withQC.R` interactively in R studio.


# Phylogenetic analysis
- Use the following shell script to submit the job.
  - Snakemake file for the complete phylogenetic analysis pipeline `scripts/phylogenetic/Snakefile`
  - Note: Only run this script AFTER finishing the breed prediction pipeline, because the following outputs are requried:
    - Meta table for breed prediction: `breed_prediction/this_wunpaired_meta.csv`
    - VCF file list created in breed prediction pipeline: `breed_prediction/vcf_file_list.txt`
    - Breed specific variants list: `breed_prediction/output_exclude_WGS/all_breed_specific_variants_13.txt`
  
```bash
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
    -p \
    --latency-wait 60 \
    --cores ${SLURM_NTASKS} \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --use-conda \
    --configfile ${config} \
    --snakefile "${project_dir}/scripts/phylogenetic/Snakefile"
```

# Breed predictions with Dog10K variants
- In this section, I will start from the rencent published SNP data generated by Dog10K project (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03023-7). Based on these variants inferred from WGS data (~20x), I developed a pipeline aiming to identify breed specific variants and perform breed validation.
- Generally speaking, the pipeline includes the following steps:
  - Download Dog10K variants (Autosomal and pseudo-autsomal region of ChrX).
  - Select high quality variants.
  - Liftover variants from canFam4 to canFam3.
  - Select breeds with at least 10 samples.
  - Calculate VAF and create VAF matrix.
  - Identify breed specific variants.
- The pipeline was written in Snakemake. Use the following shell script to submit the job.
  - Snakemake file for the complete phylogenetic analysis pipeline `scripts/wgs/Snakefile`
  - Config file `scripts/wgs/config.json`

```bash
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
    -p \
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
```

# Data backup log
```bash
# xfer4, finished
cd /project/szlab/Jingxuan_Chen
mkdir -p /project/szlab/Jingxuan_Chen/pancaner_germline
rsync -axv /scratch/jc33471/canine_tumor_0908/results/ /project/szlab/Jingxuan_Chen/pancaner_germline

# xfer4, finished
cd /project/szlab/Jingxuan_Chen
mkdir -p /project/szlab/Jingxuan_Chen/pancaner_germline_PRJNA525883
rsync -axv /scratch/jc33471/canine_tumor_rerun/results/ /project/szlab/Jingxuan_Chen/pancaner_germline_PRJNA525883

#xfer4
cd /project/szlab/Jingxuan_Chen
mkdir -p /project/szlab/Jingxuan_Chen/breed_prediction
rsync -axv /scratch/jc33471/canine_tumor/breed_prediction/ /project/szlab/Jingxuan_Chen/breed_prediction
mkdir -p /project/szlab/Jingxuan_Chen/phylogenetics
rsync -axv /scratch/jc33471/canine_tumor/phylogenetics/ /project/szlab/Jingxuan_Chen/phylogenetics
mkdir -p /project/szlab/Jingxuan_Chen/wgs_breed_prediction
rsync -axv /scratch/jc33471/canine_tumor/wgs_breed_prediction/ /project/szlab/Jingxuan_Chen/wgs_breed_prediction

# backup code to work dir, qouta exceeded
cd /work/szlab/
mkdir Jingxuan_Chen
cd Jingxuan_Chen/
git clone git@github.com:JingxuanChen7/canine_tumor_wes.git
```
