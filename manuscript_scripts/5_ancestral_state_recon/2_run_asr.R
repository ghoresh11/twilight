library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)
library(RColorBrewer)
### FUNCTIONS ###

get_range <- function(val){
  if (val == 0) {
    return("absent")
  } 
  if (val < 0.15) {
    return("rare")
  } 
  if (val < 0.95) {
    return("inter")
  }
  return("core")
}


process_b2_list <- function(l, name) {
  gained_b2_df = data.frame(Gene = names(l), X = t(data.frame(l, stringsAsFactors = F)), stringsAsFactors = F)
  gained_b2_df = gained_b2_df[-which(gained_b2_df$X <0.8), ]
  gained_b2_df$Gene = gsub(gained_b2_df$Gene, pattern = ".", replacement = "*", fixed = T)
  write.table(gained_b2_df, paste("../12_B2/",name,".csv",sep = ""), col.names = T, row.names = F, quote = F, sep = ",")
}

### MAIN #####


## perfrom acctran maximum parsimony calculations on a tree
## I don't trust the implementation from pastml, I ran it myself and the results make sense

setwd("/Users/gh11/poppunk_pangenome/11_ancestral_state_recon/")

tree = read.tree(file = "/Users/gh11/poppunk_pangenome/9_gene_properties/treeseg/tree_for_treeseg.nwk")
classification = read.table(file = "/Users/gh11/poppunk_pangenome/5_classify_genes/classification_v2.csv", sep ="\t",
                            comment.char = "", stringsAsFactors = F, header = T, quote = "")

## read all the presence/absence for the genes
## for discrete analysis
presence_absence = read.table("all.tab", sep = "\t", header = T, comment.char = "", stringsAsFactors = F, quote = "", row.names = 1)

# for continious analysis
freqs = read.table("/Users/gh11/poppunk_pangenome/4_pairwise_roary/231019_corrected/freqs.csv", sep = ",", 
                   comment.char = "", stringsAsFactors = F, header = F, quote = "", row.names = 1)

classification = classification[match(rownames(presence_absence), classification$gene),]
gene_classes = unique(classification$fill)


## reorder the frequency columns
freqs = freqs[,match(tree$tip.label,freqs[1,])]
colnames(freqs) = freqs[1,]
freqs = freqs[-1,]
freqs = freqs[match(rownames(presence_absence), rownames(freqs)),]

subtrees = subtrees(tree)
all_subtrees_empty = data.frame(id = o,
                                name = c(tree$tip.label, (length(tree$tip.label) + 1):((length(tree$tip.label) - 1)*2 + 1)),
                                size = rep(0, length(o)),
                                value_uppass = rep(0, length(o)),
                                value_acctran = rep(0, length(o)),
                                value_downpass = rep(0, length(o)),
                                value_deltran = rep(0, length(o)), stringsAsFactors = F)

## Initiate output files
o = 1:(length(tree$tip.label)*2-1)
out = c()
for (gene_class in gene_classes) {
  out = c(out, rep(gene_class, length(o)))
}

all_subtrees_out = data.frame(id = rep(all_subtrees_empty$id, length(gene_classes)),
                              name = rep(all_subtrees_empty$name, length(gene_classes)),
                              class = out,
                              increase = rep(0, length(out)),
                              decrease = rep(0, length(out)),
                              core_to_inter = rep(0, length(out)),
                              core_to_rare = rep(0, length(out)),
                              core_to_absent = rep(0, length(out)),
                              inter_to_core = rep(0, length(out)),
                              inter_to_rare = rep(0, length(out)),
                              inter_to_absent = rep(0, length(out)),
                              rare_to_core = rep(0, length(out)),
                              rare_to_inter = rep(0, length(out)),
                              rare_to_absent = rep(0, length(out)),
                              absent_to_rare = rep(0, length(out)),
                              absent_to_inter = rep(0, length(out)),
                              absent_to_core = rep(0, length(out)),
                              stringsAsFactors = F)



all_genes = data.frame(gene = rownames(presence_absence),
                       gains = rep(0, dim(presence_absence)[1]),
                       losses = rep(0, dim(presence_absence)[1]), stringsAsFactors = F)


## if testing on a particular gene
i = which(rownames(presence_absence) ==  "group_1280" )

gained_on_b2 = list()
lost_on_b2 = list()

for (i in 1:dim(presence_absence)[1]) {
  print(i)
  
  # curr_vec = as.character(c(presence_absence[i,], rep("0/1", length(tree$node.label))))
  
  # 
  ## continious
  vec = freqs[i,]
  vec[vec > 0] = 1
  
  curr_vec = as.character(c(vec,rep("0/1", length(tree$node.label)) ))
  ggtree(tree) + geom_label(aes(label = curr_vec, fill = curr_vec))
  
  res = ace(as.character(vec), phy = tree, type = "discrete", model = "ER")
  
  
  new_vec = c(as.numeric(vec), round(res$lik.anc[,2], digits = 2))
  ggtree(tree) + geom_label(aes(label = new_vec, fill = new_vec)) + scale_fill_viridis_b()
  
  ## count gains and losses
  for (j in 1:length(new_vec)) {
    curr_node = c(tree$tip.label, 48:(48+length(tree$node.label)))[j]
    p = Ancestors(tree, j, type = "parent")
    if (p == 0){ next } ## the root
    curr_class = classification$fill[i]
    row = which(all_subtrees_out$id == j & all_subtrees_out$class == curr_class)
    parent_val = get_range(new_vec[p])
    curr_val = get_range(new_vec[j])
    if (curr_val != parent_val) { ## nothing happens on branch
      column = which(colnames(all_subtrees_out) == paste(parent_val, "_to_", curr_val, sep = ""))
      all_subtrees_out[row, column] = all_subtrees_out[row, column]+ 1
    }
    if (new_vec[p] < new_vec[j]) {
      all_subtrees_out$increase[row] = all_subtrees_out$increase[row] + new_vec[j] - new_vec[p] 
      all_genes$gains[i] = all_genes$gains[i] + new_vec[j] - new_vec[p] 
      
      if (j == 79 && curr_class == "Multi-cluster core") {
        gained_on_b2[[rownames(presence_absence)[i]]] = new_vec[j] - new_vec[p] 
      }
    }
    else if (new_vec[p] > new_vec[j]) {
      all_subtrees_out$decrease[row] = all_subtrees_out$decrease[row] + new_vec[p] - new_vec[j] 
      all_genes$losses[i] = all_genes$losses[i] + new_vec[p] - new_vec[j] 
      
      if (j == 79  && curr_class == "Multi-cluster core") {
        lost_on_b2[[rownames(presence_absence)[i]]] = new_vec[p] - new_vec[j] 
      }
    }
  }
}


write.table(all_subtrees_out, file = "novel_asr_branches.csv",sep = "\t", row.names = F, col.names = T, quote = F)
write.table(all_genes, file = "novel_asr_genes.csv",sep = "\t", row.names = F, col.names = T, quote = F)
process_b2_list(gained_on_b2, "gained")
process_b2_list(lost_on_b2, "lost")





# ## Build the basic tree plot
# gene_class = "Multi-cluster intermediate"
# 
# curr = all_subtrees_out[all_subtrees_out$class %in% gene_class,]
# #curr$label = clades$Name[match(curr$id, clades$Node)]
# #curr$fill = clades$Colour[match(curr$id, clades$Node)]
# 
# ggtree(tree, aes(color=curr$increase), size = 1.2) + geom_tiplab(size = 4) + geom_tippoint()+geom_treescale(x = -0.0001, y = 44, width = 0.05)+
#   scale_color_viridis_c(name = "Increase in frequency", direction = -1) + theme(legend.position="bottom")
# 
# ggtree(tree, aes(color=curr$decrease), size = 1.2) + geom_tiplab(size = 4) + geom_tippoint()+geom_treescale(x = -0.0001, y = 44, width = 0.05)+
#   scale_color_viridis_c(name = "Decrease in frequency", direction = -1) + theme(legend.position="bottom")
