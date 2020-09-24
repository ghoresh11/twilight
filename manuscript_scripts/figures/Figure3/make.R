library(gridExtra)
library(ggfortify)
library(ape)
library(ggtree)
library(ggpubr)


# either make 2 PCA plots or the boxplots

create_pca_plot <- function(class, comp = F){
  curr_freqs = freqs[which(rownames(freqs) %in% gene_classes$Gene[gene_classes$Specific_class == class]),]
  for_pca = t(curr_freqs)
  remove = c()
  for (i in 1:dim(for_pca)[2]) {
    if (length(unique(for_pca[,i])) == 1) { ## no variation in gene
      remove = c(remove, i)
    }
  }
  if (length(remove) > 0) {
    for_pca = for_pca[,-remove]
  }
  ### PCA plot of the clusters -> what are the relationships between the clusters based on the frequencies of all genes
  freqs.pca = prcomp(for_pca , center = T)
  summary_pca = summary(freqs.pca)
  importance = round(summary_pca$importance[2,1:2] * 100, digits = 2)
  freqs.pca = data.frame(freqs.pca$x)
  freqs.pca = cbind(freqs.pca, Cluster = as.character(o))
  freqs.pca$Cluster = factor(freqs.pca$Cluster, o)
  freqs.pca$phylogroup = graphics$Phylogroup[match(freqs.pca$Cluster, graphics$Cluster)]
  p = ggplot(freqs.pca, aes(x = PC1, y = PC2, shape = phylogroup)) + geom_point(colour = "#585858", alpha = 0.7) +
    theme_classic(base_size = 12) +
    xlab(paste("PC1 (", importance[1], "%)", sep = "")) + 
    ylab(paste("PC2 (", importance[2], "%)", sep = "")) +
    scale_shape_manual(values = c(17, 3, 19, 21, 25, 8, 11, 7),
                       labels = c("A","B1","B2","C","D","E","F",expression(italic("Shigella"))),
                       name="Phylogenetic\ngroup")
  return (p)
  
}


## read gene freqs
freqs = read.table(file = "/Users/gh11/Submissions/bioresource/data/F5_freqs.csv",
                   sep = ",", header = T, stringsAsFactors = F, comment.char = "", quote = "", row.names = 1)
gene_classes = read.table("/Users/gh11/Submissions/pan_genome/Supplementary/Table_S1.csv",
                          header = T, comment.char = "", sep = ",")
graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter4//figures/cluster_graphics.csv", sep = ",",
                      header = T, comment.char = "", stringsAsFactors = F)
o = sapply(X = colnames(freqs), FUN = gsub, pattern = "X", replacement = "")


A = create_pca_plot("Multi-lineage core")
B = create_pca_plot("Intermediate and rare", res[[2]])


### use the generax output
generax = read.table("/Users/gh11/gene_distances/generax/all_stats.csv", sep = ",",
                     comment.char = "", stringsAsFactors = F, header = T)
generax$class = gene_classes$Specific_class[match(generax$Name, gene_classes$Gene)]
generax_filtered = generax[generax$class %in% c("Multi-lineage core", "Intermediate and rare"),]
generax_acc = generax[which(!generax$class %in% c("Collection core")),]
generax_acc$class = "Complete\naccessory"
generax = rbind(generax_filtered, generax_acc)

my_comparisons = list(c("Multi-lineage core", "Intermediate and rare"))
generax$class = factor(generax$class, c("Complete\naccessory","Multi-lineage core", "Intermediate and rare"))
C = ggplot(generax, aes(x=class, y = Transfer)) + geom_boxplot(width = 0.8, outlier.size = 0.5, fill = "#eeeeee")+
  ylab("Probability of transfer") +
  xlab("Distribution class") + theme_classic(base_size = 12) + 
  ggpubr::geom_signif(comparisons = my_comparisons,step_increase=0.05, tip_length = 0, map_signif_level = T)



## trees
## make a new plot where the x axis is the size of the node (from root to tip), the y axis is the 
## relative number of gain/loss events happening on the branch and that should be coloured by occurrence class
all_subtrees = read.table("/Users/gh11/poppunk_pangenome/11_ancestral_state_recon/all_subtrees.tab", sep = "\t", comment.char = "", stringsAsFactors = F, header = T)

## get the sizes of the subtrees
tree = read.tree(file = "/Users/gh11/poppunk_pangenome/9_gene_properties/treeseg/tree_for_treeseg.nwk")
the_subtrees = subtrees(tree)
clades = read.table("/Users/gh11/poppunk_pangenome/11_ancestral_state_recon/phylogroup_clades.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)


plot_tree<-function(gene_class) {
  curr = all_subtrees[all_subtrees$class %in% gene_class,]
  curr$label = clades$Name[match(curr$id, clades$Node)]
  curr$fill = clades$Colour[match(curr$id, clades$Node)]
  curr$fill[curr$label == "U"] = NA
  curr$label[curr$label == "U"] = NA
  
  
  x = as.numeric(curr$gains)
  curr$gains = (x-min(x))/(max(x)-min(x))
  
  curr$shape = graphics$Phylogroup[match(curr$name, graphics$Cluster)]
  curr$shape[is.na(curr$shape)] = "somethign"
 # pal <- pal_material(palette = "blue-grey")(10) ## change the colour pallete for the tree
 C = ggtree(tree, aes(color=curr$gains, shape = curr$shape), size = 0.6) + geom_tippoint() + 
    geom_treescale(x = -0.0001, y = 44, width = 1.3)+
    scale_color_viridis_c(direction = -1, name = "Relative number\nof gain events") +
    #  scale_color_gradientn(colours = pal,name = "Gain events", guide = F) + 
    theme(legend.position="bottom")+
    scale_fill_manual(values = clades$Colour, guide = F) +xlim(NA, 0.3)+
   scale_shape_manual(values = c(17, 3, 19, 21, 25, 8, 11, 7,20)
                      , name="Phylogenetic\ngroup", 
                      labels = c("A","B1","B2","C","D","E","F",expression(italic("Shigella"))),
                                 guide = F)
  
  
  
  D = ggtree(tree, aes(color=curr$losses), size = 0.6) + geom_tippoint() + 
    geom_treescale(x = -0.0001, y = 44, width =1.3)+
   # scale_color_gradientn(colours = pal,name = "Loss events") +  theme(legend.position="bottom") +
    scale_fill_manual(values = clades$Colour, guide = F) +xlim(NA, 0.3)
  
  return(list(C,D))
}

D = plot_tree("Multi-cluster core")[[1]]
E =   plot_tree("Intermediate and rare")[[1]] 


legend = as_ggplot(get_legend(A +theme(legend.position = "right")))
A = A + theme(legend.position = "None")
# grid.arrange(A,B,C + coord_flip(),D,E, legend, layout_matrix= rbind(c(1,2,4,5),
#                                                                     c(6,6,4,5),
#                                                                     c(NA,3,4,5)))

tree_legend = as_ggplot(get_legend(D ))


# ## panel for the function?
# functions = read.table("/Users/gh11/poppunk_pangenome/9_gene_properties/eggnog/eggnog_results.emapper.annotations",
#                        sep = "\t",  stringsAsFactors = F, quote = "")
# classification = read.table("/Users/gh11/poppunk_pangenome/5_classify_genes/classification_v2.csv", sep ="\t", comment.char = "",stringsAsFactors = F, quote = "", header = T)
# functions$class = classification$fill[match(functions$V1, classification$gene)]
# functions = functions[which(functions$class %in% c("Multi-cluster core","Intermediate and rare")),]
# cog_descs =  read.table("/Users/gh11/poppunk_pangenome/9_gene_properties/eggnog/cog_descs.csv",
#                         sep = ",",  stringsAsFactors = F, quote = "", header = T, comment.char = "")
# functions$descs = cog_descs$Title[match(functions$V12, cog_descs$COG)]
# functions$descs[is.na(functions$descs)] = "UNASSIGNED"
# num_inter = length(which(functions$class == "Intermediate and rare"))
# num_core = length(which(functions$class == "Multi-cluster core"))
# functions = data.frame(table(functions$class, functions$descs))
# functions$Freq[functions$Var1 == "Intermediate and rare"] = functions$Freq[functions$Var1 == "Intermediate and rare"]/num_inter
# functions$Freq[functions$Var1 == "Multi-cluster core"] = functions$Freq[functions$Var1 == "Multi-cluster core"]/num_core
# ggplot(functions, aes(x = Var2, fill = Var1, y= Freq)) + geom_bar(position = "dodge", stat = "identity") 


grid.arrange(A,B+ theme(legend.position = "None"),C,D + theme(legend.position = "None"),
             E+ theme(legend.position = "None"),legend,  layout_matrix= rbind(c(6,3,4,5),
                                                                              c(6,1,4,5),
                                                                              c(NA,2,4,5)))




