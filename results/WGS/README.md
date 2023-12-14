# Description for result files under folder `results/WGS/`

## Figures
- `results/WGS/breeds_heatmap_wgsmain_CDS_beforeQC4567.png` and `results/WGS/breeds_heatmap_wgsmain_CDS_beforeQC9999.png`: The heatmap shows the clustering of 278 dogs, using variant allele frequency (VAF) values of the 6364 breed-specific germline base substitution and small indel variants. These variants were discovered with the WGS dataset generated by Dog10K project. The “Provided Breed” bar indicate the breed for each dog provided by the source study. 23 breeds with more than 10 dogs are selected for breed specific variants identification. Two figures are generated with seeds 4567 and 9999, respectively.
Figures are generated by running `scripts/wgs/breeds_joint_heatmap_beforeQC_wgs.R` interactively in RStudio.


- `results/WGS/breeds_heatmap_wgsassignment_beforeQC4567.png`: The heatmap shows the clustering of 1590 dogs, using variant allele frequency (VAF) values of the 6364 breed-specific germline base substitution and small indel variants. These variants were discovered with 278 dogs in the WGS dataset generated by Dog10K project. The “Provided Breed” bar indicate the breed for each dog provided by the source study. 23 breeds with more than 10 dogs are selected for breed specific variants identification. Breed dogs other than 23 selected breeds are labeled as "Unknown/Missing". The figure is generated with seed 4567.
Figures are generated by running `scripts/wgs/breeds_joint_heatmap_beforeQC_wgs.R` interactively in RStudio.

- `results/WGS/breeds_heatmap_mergemain_beforeQC.png` and `results/WGS/breeds_heatmap_mergeassignment_beforeQC.png`: The heatmap shows the clustering of 965 dogs (783 with breed provided, 182 without breed information), using variant allele frequency (VAF) values of the 5945 breed-specific germline base substitution and small indel variants. For dogs with paired Normal and Tumor samples, the normal one is used; Dogs only with Normal or Tumor samples are also used; Healthy dogs in Dog10K dataset are included. Breed specific variants are merged from WES and WGS indentified sets. The “Provided Breed” bar and the “Disease” bar respectively indicate the breed and tumor type of each dog provided by the source studies. 12 breeds from WES as well as 23 breeds from WGS with more than 10 dogs are selected for breed specific variants identification.
Figures are generated by running `scripts/breed_prediction/breeds_joint_heatmap_beforeQC_withWGS.R` interactively in RStudio.

## Text files

- `results/WGS/breed_specific_variants_concat_wgs.txt` and `results/WGS/breed_unique_variants_concat_wgs.txt`: List of breed specific variants and breed unique identified with WGS dataset for 23 breeds, respectively. Note that whole genome-wide variants are included in this list.

- `results/WGS/breed_specific_variants_CDS_wgs.txt`: List of breed specific variants identified with WGS dataset for 23 breeds. Note that Only variants in CDS regions are included in this list.

- `results/WGS/main_clusters_wgs.txt` and `results/WGS/assignment_clusters_wgs.txt`: Clustering results of 278 dogs in WGS dataset without or with other breeds, i.e., corresponding figures are `results/WGS/breeds_heatmap_wgsmain_CDS_beforeQC4567.png` and `results/WGS/breeds_heatmap_wgsassignment_beforeQC4567.png`, respectively. Generated by running `scripts/wgs/breeds_joint_heatmap_beforeQC_wgs.R` interactively in RStudio.

- `results/WGS/all_breed_specific_variants_merge.txt`: List of breed specific variants merged from WES and WGS indentified sets. Note that this list has NOT been filtered based on proportion of NA. It is one of the input file for `scripts/breed_prediction/breeds_joint_heatmap_beforeQC_withWGS.R`. Further filtering of variants is done in the R script. The specific list of variants visualized in heatmap `results/WGS/breeds_heatmap_mergeassignment_beforeQC.png` can be found at `results/WGS/heatmap_variants_mergeWGS.txt`.

- `results/WGS/main_clusters_merge.txt` and `results/WGS/assignment_clusters_merge.txt`: Clustering results of merged WGS and WES datasets without or with unknown breed samples. Note that only 278 dogs (for 23 breeds) in WGS dataset were included; while most 12 breed plus unknown breed dogs in WES datasets were included. Corresponding figures are `results/WGS/breeds_heatmap_mergemain_beforeQC.png` and `results/WGS/breeds_heatmap_mergeassignment_beforeQC.png`, respectively. Generated by running `scripts/breed_prediction/breeds_joint_heatmap_beforeQC_withWGS.R` interactively in RStudio.