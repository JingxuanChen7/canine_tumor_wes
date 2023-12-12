# Description for result files under folder `results/phylogenetics/`

## Figures
- All trees were computed and visualized using pipeline `scripts/phylogenetic/Snakefile`.

- `results/phylogenetics/breedSample_breedSpecific.pdf`: Figure shows a neighbor-joining tree reconstructed using 2,783 breed specific variants identified from WES datasets, across 474 dogs with breed provided. Only normal dogs are included in this analysis. The tree is unrooted. Breed information is annotated at tip points as well as the color stripe. Tree scale is annotated as x-axis. 


- `results/phylogenetics/breedPlusMissingSample_breedSpecific.pdf`: Figure shows a neighbor-joining tree reconstructed using 2,783 breed specific variants identified from WES datasets, across 656 dogs (474 dogs with breed provided and 182 dogs with no breed provided). Only normal dogs are included in this analysis. The tree is unrooted. Breed information is annotated at tip points as well as the color stripe. Tree scale is annotated as x-axis.

- `results/phylogenetics/breedSample_all.pdf`: Figure shows a neighbor-joining tree reconstructed using 321961 germline variants identified across the whole exome, across 474 dogs. Only normal samples are included in this analysis. The tree is re-rooted with Shih Tzu breed. Breed information is annotated at tip points as well as the color stripe. Tree scale is annotated as x-axis.

- `results/phylogenetics/breedSample_all_consensus.pdf`: Figure shows a consensus neighbor-joining tree reconstructed using 321961 germline variants identified across the whole exome, across 474 dogs. To determine the significance of branch placement in the cladogram, the dataset was resampled 100 times by pulling a random 10% of the variants to make 100 distance matrices. The cladograms created from each of the random variant-set matrices were combined using consense in the PHYLIP programs. Only normal samples are included in this analysis. The tree is re-rooted with Shih Tzu breed. Breed information is annotated at tip points as well as the color stripe. Branch support values greater than 30 are annotated on branches.