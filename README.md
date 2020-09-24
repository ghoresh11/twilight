# The twilight zone of the pan-genome

Apply an analysis similar to the one in DOI...

## Requirements

- R

R packages:
- data.table
- ggplot2
- optparse



## To classify your own pan-genome:

To apply a similar analysis on the output of a pan-genome analysis using tools such as Roary (1) or Panaroo (2), you must provide:

`-p, --presence_absence`:  The `gene_presence_absence.Rtab` file, which states for each gene the presence or absence in each genome used (tab separated).

`-g, --grouping`: A tab separated file, with no header, which states the group of each genome.

The names of the genomes in the grouping file must match the names of the genomes in presence absence file.

Examples are available in the directory `test_sets/`

### Usage

`Rscript classify_genes.R -p gene_presence_absence.Rtab -g groups.tab`

### Optional params:

`-o, --out`: output directory name (default = "out").

`-m, --min_size` : ignore groups with fewer than `min_size` genomes (default = 10).

`-c, --core_threshold` : Threshold used to define a core gene within each group (default = 0.95).

`-r, --rare_threshold` : Threshold used to define a rare gene within each group (default = 0.15).

`-h, --help` : print help message and quit.


## Outputs:

A number of files are generated in the `out` directory as follows:

1. `classification.tab` = A table with the classification of all genes based on the new definitions. Including how many times each gene was observed as core, intermediate or rare and in which groups.

2. `frequencies.csv` = A table stating the precise frequency of each gene in each of the groups.

3. `genes_per_isolate.tab` = A table stating for each genome in the collection, how many genes from each of the distribution classes were present in its genome. Useful for measuring a typical genome in the collection and per group.

4. `plots` = a directory containing the same figures from the manuscript. The `typical_per_class` directory has a plot for each distribution class, and shows how many genes from that class were present in a single genome from each group. The `pca_per_class` presents a PCA analysis on the gene frequencies for each distribution class.

