library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)

## perfrom acctran maximum parsimony calculations on a tree
## I don't trust the implementation from pastml, I ran it myself and the results make sense

debug = F

setwd("/Users/gh11/poppunk_pangenome/11_ancestral_state_recon/")

tree = read.tree(file = "/Users/gh11/poppunk_pangenome/9_gene_properties/treeseg/tree_for_treeseg.nwk")
classification = read.table(file = "/Users/gh11/poppunk_pangenome/5_classify_genes/classification_v2.csv", sep ="\t",
                            comment.char = "", stringsAsFactors = F, header = T, quote = "")

## read all the presence/absence for the genes
presence_absence = read.table("all.tab", sep = "\t", header = T, comment.char = "", stringsAsFactors = F, quote = "")
classification = classification[match(rownames(presence_absence), classification$gene),]

gene_classes = unique(classification$fill)

if (debug) {
  tree = read.tree("dummy.nwk")
  curr = c(1,1,0,0,1)
  ggtree(tree) + geom_label(aes(label = c(curr, rep("?", length(curr)-1))))
}


o = 1:(length(tree$tip.label)*2-1)

subtrees = subtrees(tree)
all_subtrees_empty = data.frame(id = o,
                                name = c(tree$tip.label, (length(tree$tip.label) + 1):((length(tree$tip.label) - 1)*2 + 1)),
                                size = rep(0, length(o)),
                                value_uppass = rep(0, length(o)),
                                value_acctran = rep(0, length(o)),
                                value_downpass = rep(0, length(o)),
                                value_deltran = rep(0, length(o)), stringsAsFactors = F)

out = c()
for (gene_class in gene_classes) {
  out = c(out, rep(gene_class, length(o)))
}

all_subtrees_out = data.frame(id = rep(all_subtrees_empty$id, length(gene_classes)),
                              name = rep(all_subtrees_empty$name, length(gene_classes)),
                              class = out,
                              gains = rep(0, length(out)),
                              losses = rep(0, length(out)), stringsAsFactors = F)



all_genes = data.frame(gene = rownames(presence_absence),
                       gains = rep(0, dim(presence_absence)[1]),
                       losses = rep(0, dim(presence_absence)[1]), stringsAsFactors = F)



## do this once, get the size of all subtrees
for (i in 1:length(o)){
  if (i < length(tree$tip.label) + 1) {
    all_subtrees_empty$size[i] = 1
    next
  }
  curr_tree = subtrees[[i - length(tree$tip.label)]]
  all_subtrees_empty$size[i] = length(curr_tree$tip.label)
}
all_subtrees_empty = all_subtrees_empty[order(all_subtrees_empty$size),]


common_states <- function(vec) {
  num_1s = 0
  num_0s = 0
  for (val in vec) {
    if (val == 1) {
      num_1s = num_1s + 1
    } else if (val == 0) {
      num_0s = num_0s + 1
    } 
  }
  if (num_1s > num_0s) {
    return (1)
  } 
  if (num_0s > num_1s) {
    return (0)
  }
  return("0/1")
}


## here - add loop for all genes
for (curr_index in 1:dim(presence_absence)[1]) {
  print(curr_index)
  curr = as.numeric(presence_absence[curr_index,])
  
  ## refresh
  all_subtrees = all_subtrees_empty
  
  ## UPPASS
  for (i in 1:length(o)) {
    if (all_subtrees$size[i] == 1) {
      all_subtrees$value_uppass[i] = curr[all_subtrees$id[i]]
      next
    }
    ## get the two children of the current node
    children = Children(tree, all_subtrees$id[i])
    val1 = all_subtrees$value_uppass[all_subtrees$id == children[1]]
    val2 = all_subtrees$value_uppass[all_subtrees$id == children[2]]
    all_subtrees$value_uppass[i] = common_states(c(val1, val2))
  }
  
  all_subtrees = all_subtrees[order(all_subtrees$id),]
  ggtree(tree) + geom_label(aes(label =  all_subtrees$value_uppass, fill = all_subtrees$value_uppass)) +
    ggtitle("UPPASS")
  all_subtrees = all_subtrees[order(all_subtrees$size),]
  
  all_subtrees$value_acctran = rep("0/1", dim(all_subtrees)[1])
  for (i in 1:length(o)) {
    if (all_subtrees$size[i] == 1) {
      all_subtrees$value_acctran[i] = curr[all_subtrees$id[i]]
      next
    }
  }
  
  ggtree(tree) + geom_label(aes(label =  all_subtrees$value_acctran, fill = all_subtrees$value_acctran)) +
    scale_fill_manual(values = c("#d95f02","#d3d3d3","#7570b3"))
  
  ## ACCTRAN - From root to tip, following uppass
  all_subtrees$value_acctran = all_subtrees$value_uppass
  for (i in rev(all_subtrees$id)) {
    if (all_subtrees$value_acctran[which(all_subtrees$id == i)] != "0/1") {
      next
    }
    if (all_subtrees$size[which(all_subtrees$id == i)] == 47) { ## can't do anything about the root
      next
    }
    parent = Ancestors(tree, i, type = "parent")
    all_subtrees$value_acctran[all_subtrees$id == i] = all_subtrees$value_acctran[all_subtrees$id == parent]
  }
  
  all_subtrees = all_subtrees[order(all_subtrees$id),]
  ggtree(tree) + geom_label(aes(label =  all_subtrees$value_acctran, fill = all_subtrees$value_acctran)) 
  all_subtrees = all_subtrees[order(all_subtrees$size),]
  
  
  
  # # ## DOWNPASS
  # #DOWNPASS traverses the tree starting from the root and going down till the tips,
  # #and for each node combines the state information from its supertree and its subtree (calculated at UPPASS).
  # all_subtrees$value_downpass = c(curr, rep(NA, length(tree$tip.label) - 1))
  # downpass <- function(i) {
  #   if (all_subtrees$size[all_subtrees$id == i] == 1) { ## don't do anything for a tip
  #     return(all_subtrees)
  #   }
  #   children = Children(tree, i) ## run on children
  #   all_subtrees = downpass(children[1])
  #   all_subtrees = downpass(children[2])
  #   if (all_subtrees$size[all_subtrees$id == i] != length(tree$tip.label)) { ## not the root, choose the most parsimonous result
  #     val_a = strsplit(all_subtrees$value_uppass[all_subtrees$id == children[1]], split = "/", fixed = T)[[1]]
  #     val_b = strsplit(x = all_subtrees$value_uppass[all_subtrees$id == children[2]], split = "/", fixed = T)[[1]]
  #     intersection = intersect(val_a, val_b)
  #     print("NAES")
  #     print(val_a)
  #     print(val_b)
  #     print(intersection)
  #     if (length(intersection) == 0) {
  #       all_subtrees$value_downpass[all_subtrees$id == i] = paste(sort(union(val_a, val_b)),sep="", collapse = "/")
  #     } else {
  #       print(paste(sort(intersection),sep="", collapse = "/"))
  #       all_subtrees$value_downpass[all_subtrees$id == i] = paste(sort(intersection),sep="", collapse = "/")
  #     }
  #   }
  #   children = Children(tree, i) ## run on children
  #   
  #   return(all_subtrees)
  # }
  # 
  # all_subtrees = downpass(length(tree$tip.label) + 1)
  # 
  # all_subtrees = all_subtrees[order(all_subtrees$id),]
  # ggtree(tree) + geom_label(aes(label =  all_subtrees$value_downpass, fill = all_subtrees$value_downpass)) +
  #   ggtitle("DOWNPASS")
  # all_subtrees = all_subtrees[order(all_subtrees$size),]
  # 
  # all_subtrees$value_deltran = as.character(all_subtrees$value_downpass)
  # deltran <- function(i) {
  #   print(i)
  #   
  #   if (all_subtrees$size[all_subtrees$id == i] == 1) { ## don't do anything for a tip
  #     return(all_subtrees)
  #   }
  #   if (all_subtrees$size[all_subtrees$id == i] != 47) { ## calculate for non-root
  #     p = Ancestors(tree, i, type = "parent")
  #     val_p = strsplit(all_subtrees$value_deltran[all_subtrees$id == p], split = "/", fixed = T)[[1]]
  #     curr_val = strsplit(x = all_subtrees$value_deltran[all_subtrees$id == i], split = "/", fixed = T)[[1]]
  #     #all_subtrees$value_deltran = common_states(c(val_p, curr_val))
  #     if (length(intersect(curr_val, val_p)) == 1 ) {
  #       all_subtrees$value_deltran[all_subtrees$id == i] = intersect(curr_val, val_p)
  #     } else if (length(intersect(curr_val, val_p)) == 2 ) {
  #       all_subtrees$value_deltran[all_subtrees$id == i] = "0/1"
  #     }
  #   }
  #   children = Children(tree, i) ## run on children
  #   all_subtrees = deltran(children[1])
  #   all_subtrees = deltran(children[2])
  #   return(all_subtrees)
  # }
  # 
  # all_subtrees = deltran(48)
  # 
  # 
  # 
  # ## DELTRAN - From tip to root
  # all_subtrees$value_deltran = all_subtrees$value_downpass
  # for (i in all_subtrees$id) {
  #   if (all_subtrees$value_deltran[all_subtrees$id == i] != "0/1") {
  #     next
  #   }
  #   if (all_subtrees$size[all_subtrees$id == i] == 47) {
  #     next
  #   }
  #   parent = Ancestors(tree, i, type = "parent")
  #   all_subtrees$value_deltran[all_subtrees$id == i] = 
  #     common_states(c(all_subtrees$value_deltran[all_subtrees$id == parent], all_subtrees$value_deltran[all_subtrees$id == i]))
  # }
  
  all_subtrees = all_subtrees[order(all_subtrees$id),]
  ggtree(tree) + geom_label(aes(label =  all_subtrees$value_deltran, fill = all_subtrees$value_deltran)) +
    ggtitle("DELTRAN")
  all_subtrees = all_subtrees[order(all_subtrees$size),]
  
  ## count gains and losses
  for (i in all_subtrees$id) {
    if (all_subtrees$size[all_subtrees$id == i] == 47) { ## root
      next
    }
    p = Ancestors(tree, i, type = "parent")
    curr_class = classification$fill[curr_index]
    row = which(all_subtrees_out$id == i & all_subtrees_out$class == curr_class)
    if (all_subtrees$value_acctran[all_subtrees$id == i] == 0 && all_subtrees$value_acctran[all_subtrees$id == p] == 1) {
      all_genes$losses[curr_index] = all_genes$losses[curr_index] + 1
      all_subtrees_out$losses[row] = all_subtrees_out$losses[row] + 1
    } else if (all_subtrees$value_acctran[all_subtrees$id == i] == 1 && all_subtrees$value_acctran[all_subtrees$id == p] == 0) {
      all_genes$gains[curr_index] = all_genes$gains[curr_index] + 1
      all_subtrees_out$gains[row] = all_subtrees_out$gains[row] + 1
    }
  }
}



# all_subtrees_out = all_subtrees[order(all_subtrees$id),]
# ggtree(tree) + geom_label(aes(label =  all_subtrees$gains, fill = all_subtrees$gains)) +
#   ggtitle("ACCTRAN - gains") + scale_fill_continuous(low = "yellow", high = "red")
# 
# ggtree(tree) + geom_label(aes(label =  all_subtrees$gains, fill = all_subtrees$losses)) +
#   ggtitle("ACCTRAN - losses") + scale_fill_continuous(low = "yellow", high = "red")

all_subtrees_empty = all_subtrees_empty[order(all_subtrees_empty$id),]
all_subtrees_out$size = rep(all_subtrees_empty$size, length(gene_classes))

# write.table(all_subtrees_out, file = "all_subtrees.tab", sep = "\t", col.names = T, row.names = F, quote = F )
# write.table(all_genes, file = "all_genes.tab", sep = "\t", col.names = T, row.names = F, quote = F )
