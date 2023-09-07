# run pipeline (germline & depthofcoverage)
- Install conda environments to run the pipeline on cluster
```
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/wes_env.yml --name wes_env
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/annovar_env.yml --name annovar_env
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/mutect2_env.yml --name mutect2_env
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/strelka_env.yml --name strelka_env
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/java17.yml --name java17
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/data_collection.yml --name data_collection
```
- Shell script to submit individual jobs:


# merge individual VCFs
- *To be tested* https://github.com/ZhaoS-Lab/breed_prediction/blob/main/combine_germline_variants.sh


# post-processing and identify breed-specific variants
- Install conda environment for breed prediction pipeline
```
mamba env create --force -f /home/jc33471/canine_tumor_wes/scripts/envs/breed_prediction.yml --name breed_prediction
```

- Submit job using the following script
```
#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=breed
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=60G
#SBATCH --time=100:00:00
#SBATCH --mail-user=jc33471@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/scratch/jc33471/canine_tumor_test/breed.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate breed_prediction

project_dir="/home/${USER}/canine_tumor_wes"
run_dir="/scratch/${USER}/canine_tumor_test"
mkdir -p ${run_dir}
cd ${run_dir}

# convert VAF at low coverage regions to NA
Rscript --vanilla ${project_dir}/scripts/breed_prediction/germline_VAF_reset_low_coverage.R \
  "/scratch/jc33471/canine_tumor_test/breed_prediction/PanCancer_disc_val_merged_germline_VAF_01_01_2021.txt.gz" \
  "/scratch/jc33471/canine_tumor_test/breed_prediction/PanCancer_disc_val_merged_germline_depths_01_01_2021.txt.gz" \
  "/scratch/jc33471/canine_tumor_test/breed_prediction/germline_VAF_matrix.reset_low_coverage.txt.gz"

# identify breed specific variants
Rscript --vanilla ${project_dir}/scripts/breed_prediction/breed_specific_variants.R \
  "/home/jc33471/canine_tumor_wes/scripts/breed_prediction/build_sample_meta_data.R" \
  "/scratch/jc33471/canine_tumor_test/breed_prediction/germline_VAF_matrix.reset_low_coverage.txt.gz" \
  "/scratch/jc33471/canine_tumor_test/breed_prediction/output_exclude_WGS/breed_unique_variants.txt" \
  "/scratch/jc33471/canine_tumor_test/breed_prediction/output_exclude_WGS/breed_enriched_variants.txt" \
  "/scratch/jc33471/canine_tumor_test/breed_prediction/output_exclude_WGS/all_breed_specific_variants.txt" \
  "/scratch/jc33471/canine_tumor_test/breed_prediction/breed_prediction_metadata.txt"
```

# heatmap visualization
- Firstly run Run R script `scripts/breed_prediction/breeds_joint_heatmap_beforeQC.R` interactively in R studio.  (command line interface not implemented yet)
  - This step is to generate heatmap for breed validation/prediction.
  - Need to mannually add info to metatable to generate the final heatmap.
- Run R script `scripts/breed_prediction/breeds_joint_heatmap_withQC.R` interactively in R studio.
