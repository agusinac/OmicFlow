# Automated-Omics-Analysis

## Running pipeline from 00_main.R as follows:

```
Rscript 00_main.R --metadata ../data/metadata/metadata.tsv --biom ../data/BiotaViz/absolute-table-with-taxonomy.biom --tree ../data/phylogeny/phylogenetic_tree/tree.nwk --refseq ../data/denoise/representative_sequences/unfiltered/sequences.fasta

```

## Running 05_Model-RDA.R as Rscript as follows:

```
Rscript 05_Model-RDA.R --inomics ../data/input.xlsx --osheet genera --inpheno ../data/input.xlsx --psheet pheno --contrast RANKSTAT_BBCBaCo --log 10 --output ../results --pairwise TRUE

```