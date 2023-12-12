- This repository includes all related files (codes/figures/presentations/etc.) for the canine breed prediction project during Fall 2023 in Zhao lab.
- Folder structure explained below:
  - `scripts/`: All codes for data processing and analysis. `scripts/README.md` includes detailed documentation for each script.
  - `results/`: Figures and summarized result files. Detailed README files are included in each subfolder.
  - `metadata/`: Meta tables from SRA run selector and merged meta tables to run the pipeline. README file included.
  - `doc/`: Slides for lab meetings.

# Results backup
- Germline variants (in VCF), DepthOfCoverage (in BED), and QC results are backed up on cluster at `/project/szlab/Jingxuan_Chen`.
  - `/project/szlab/Jingxuan_Chen/pancaner_germline`: Results for new WES datasets as included in meta table `metadata/https://github.com/JingxuanChen7/canine_tumor_wes/blob/main/metadata/data_collection_new.csv`.
  - `/project/szlab/Jingxuan_Chen/pancaner_germline_PRJNA525883`: Results for re-run of PRJNA525883, which has already been included in the original publication, but indiviual results were removed.

- Result files for breed prediction are backed up on cluster at `/project/szlab/Jingxuan_Chen/breed_prediction`.

- Result files for phylogenetic analysis are backed up on cluster at `/project/szlab/Jingxuan_Chen/phylogenetics`.

- Result files for WGS related results (breed specific variant identification, merge with WES results, etc.) are backed up on cluster at `/project/szlab/Jingxuan_Chen/wgs_breed_prediction`.

- This repository is backed up on cluster at `/project/szlab/Jingxuan_Chen/breed_prediction_code`.
