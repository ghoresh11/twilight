library(ape)

freqs = read.table("/Users/gh11/poppunk_pangenome/4_pairwise_roary/231019_corrected/freqs.csv", sep = ",",
                   comment.char = "", stringsAsFactors = F, header = T, quote = "", row.names = 1)

classification = read.table("/Users/gh11/poppunk_pangenome/5_classify_genes/classification_v2.csv", sep = "\t",
                            comment.char = "", stringsAsFactors = F, header = T, quote = "")

freqs[freqs>0] = 1
freqs = freqs[-which(rowSums(freqs) == 1),]
freqs = freqs[-which(rowSums(freqs) == 47),]
freqs = freqs[-which(rowSums(freqs) == 0),]

freqs = freqs[,-which(colnames(freqs) == "X50")]

freqs = t(freqs)
rownames(freqs) = gsub(rownames(freqs), pattern = "X", replacement = "", fixed = T)

write.table(x = match(colnames(freqs),classification$gene), file = "/Users/gh11/Submissions/my_thesis/Chapter5/prep_v2/corrections/ancestral_state_recon/order.txt",
            quote = F, col.names = F, row.names = F)

# colnames(freqs) = gsub(colnames(freqs), pattern = "*", replacement = "XX", fixed = T) ## pastml fails when there are special chars in the colnames
# colnames(freqs) = gsub(colnames(freqs), pattern = "'", replacement = "ZZ", fixed = T)
# colnames(freqs) = gsub(colnames(freqs), pattern = "(", replacement = "YY", fixed = T)
# colnames(freqs) = gsub(colnames(freqs), pattern = ")", replacement = "UU", fixed = T)
# colnames(freqs) = gsub(colnames(freqs), pattern = "/", replacement = "VV", fixed = T)

tree = read.tree("/Users/gh11/Submissions/my_thesis/Chapter5/prep_v2/corrections/ancestral_state_recon/named.tree_tree_for_treeseg.nwk")
freqs = freqs[match(tree$tip.label, rownames(freqs)),]

write.table(t(freqs), 
            file = paste("/Users/gh11/Submissions/my_thesis/Chapter5/prep_v2/corrections/ancestral_state_recon/all.tab",sep = ""),
            sep = "\t", col.names = T, row.names = T, quote = F )

## save freqs into chunks of 1000 to see where it fails
for (i in seq(from = 1, to = dim(freqs)[2], by = 1000)) {
  write.table(freqs[,i:min(i+999, dim(freqs)[2])], 
              file = paste("/Users/gh11/Submissions/my_thesis/Chapter5/prep_v2/corrections/ancestral_state_recon/inputs/",i, ".tab",sep = ""),
              sep = "\t", col.names = T, row.names = T, quote = F )
  
}

## if converting to FASTA and using pyjar by Simon Harris
freqs[freqs == 0] = "C"
freqs[freqs == 1] = "A"
classification = classification[match(colnames(freqs),classification$gene),]

for (gene_class in unique(classification$fill)) {
  curr = freqs[, which(classification$fill == gene_class)]
  print(gene_class)
  print(dim(curr))
  filename = paste("/Users/gh11/Submissions/my_thesis/Chapter5/prep_v2/corrections/pyjar_asr/", gsub(x = gene_class, pattern = " ", replacement = "-"), ".fna", sep = "")
  if (file.exists(filename))
    #Delete file if it exists
    file.remove(filename)
  
  for (i in 1:dim(curr)[1]) {
    write(paste(">", rownames(curr)[i], sep = "", collapse = ""), file = filename, append = T)
    write(paste(curr[i,], sep = "", collapse = ""), file = filename, append = T)
  }
}
