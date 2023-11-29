- This folder includes all related files (codes/figures/presentations/etc.) for the canine breed prediction project during Fall 2023 in Zhao lab.
- Folder structure explained below:
  -`scripts/`: All codes for data processing and analysis. `scripts/README.md` includes detailed documentation for each script.
  - `results/`: Figures and summarized result files.
  - `metadata/`: Meta tables from SRA run selector and merged meta tables to run the pipeline.
  - `doc/`: Slides for lab meetings.

# Results backup
- Germline variants (in VCF), DepthOfCoverage (in BED), and QC results are backed up on cluster at `/project/szlab/Jingxuan_Chen`.
  - `/project/szlab/Jingxuan_Chen/pancaner_germline`: Results for new WES datasets as included in meta table `metadata/https://github.com/JingxuanChen7/canine_tumor_wes/blob/main/metadata/data_collection_new.csv`.
  - `/project/szlab/Jingxuan_Chen/pancaner_germline_PRJNA525883`: Results for re-run of PRJNA525883, which has already been included in the original publication.


# Write before starting...
- This work is based on a previous publication from Zhao lab: https://www.nature.com/articles/s41467-021-24836-9
- The overall plan is:
  1. Replicate analysis in the publication.
  2. Add more recent data into the pipeline.
## 1. Replicate analysis in Alsaihati et al. 2021
- Public codes for the paper: https://github.com/ZhaoS-Lab/breed_prediction
- Pipeline from Kunlin is checked here: https://github.com/JingxuanChen7/canine_tumor_wes/blob/main/scripts/004_MC_WES_complete_pipline.sh
## 2. Add more recent WES/WGS data
- Add public datasets published after 2020.
- Add WES/WXS/exome data. WGS could be potentially included if everything works fine.
- Note:
  - Merge samples with same sample accessions.
  - dog10K (not available now, some data in other publications): http://www.dog10kgenomes.org/
  - A search example for new datasets: https://www.ncbi.nlm.nih.gov/bioproject/?term=(dog%20exome)%20AND%20bioproject_sra%5bfilter%5d%20NOT%20bioproject_gap%5bfilter
  - Another potential data source, IDog: https://ngdc.cncb.ac.cn/idog/
  - Another potential data source, Dog Genome Project: https://www.broadinstitute.org/scientific-community/science/projects/mammals-models/dog/dog-genome-links
  - Dog SNP database: https://academic.oup.com/nar/article/51/D1/D816/6775385
  - New reference genome (keep using canFam3 for this project): https://www.nature.com/articles/s42003-021-01698-x
  - ONLY PAIRED-END data included in this study.
