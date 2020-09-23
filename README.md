# twilight
All scripts used for the analysis of the manuscript for analysing the twilight zone of the pan-genome.

## Example usage
To apply a similar analysis on the output of a pan-genome analysis using tools such as Roary (1) or Panaroo (2), you must provide:
1. The `presence_absence.Rtab` file, which states for each gene the presence or absence in each genome used.
2. A file which states, for each genome, its lineage (tab separated).

The names of the genomes in file 2 must match the names of the genomes in file 1.

## To run:

`Rscript twilight.R presence_absence.Rtab groups.tab`

### Optional params:

`--min_size` : ignore groups with fewer than `min_size` genomes.

## Output:

-- Figures

-- Classification table

-- Typical E. coli table
