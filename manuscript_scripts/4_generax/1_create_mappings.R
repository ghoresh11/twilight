library(ape)
library(reshape2)
library(harrietr)

create_mappings <- function(gene_name) {
  old_name = old_to_new$Old[old_to_new$New == gene_name]
  tree_file = paste("/Users/gh11/poppunk_pangenome/trees/", gene_name, ".treefile" ,sep = "")
  gene_tree = read.tree(tree_file)
  gene_distance <-cophenetic(gene_tree)
  gene_distance = melt_dist(gene_distance)
  
  mapping = data.frame( sapply(sapply(gene_tree$tip.label, strsplit, split = "_", fixed = T), head, n = 1), gene_tree$tip.label)
  write.table(mapping, file = paste("generax/links/",gene_name, ".link", sep = ""),sep = ":", col.names = F, row.names = F, quote = F)
}


old_to_new = read.table("/Users/gh11/Submissions/bioresource/data/scripts/change_gene_names/old_to_new.csv",
                        header = T, quote = "", stringsAsFactors = F, sep = ",")
all_genes = sapply(list.files("/Users/gh11/poppunk_pangenome/trees/", pattern = "*.treefile"), gsub, pattern = ".treefile", replacement = "")
sapply(all_genes, create_mappings)
for (g in all_genes) {
  print(g)
  create_mappings(g)
}
