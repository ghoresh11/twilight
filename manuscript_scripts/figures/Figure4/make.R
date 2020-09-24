library(harrietr)
library(reshape2)

setwd("/Users/gh11/poppunk_pangenome/10_gene_sharing/")

my_comparisons = list(c("All other\nlineages","12"), c("All other\nlineages", "40"))


######## PANEL A ########
## shared genes
edges = read.table("num_shared_gene_v2.csv", sep = "\t", comment.char = "",
                   stringsAsFactors = F, quote = "", header = T)
edges = edges[edges$Gene_class == "Intermediate and rare",]

## add size info
cluster_sizes = read.table("../2_dists_roary_analysis/cluster_sizes_updated.csv", sep = ",", comment.char = "", header = T,stringsAsFactors = F)
edges$sizeA = cluster_sizes$Size[match(edges$ClusterA, cluster_sizes$Cluster)]
edges$sizeB = cluster_sizes$Size[match(edges$ClusterB, cluster_sizes$Cluster)]
edges  = transform(edges, min = pmin(edges$sizeA , edges$sizeB))
edges  = transform(edges, max = pmax(edges$sizeA , edges$sizeB))


## add distance
tree = read.tree("/Users/gh11/gene_distances/raxml_tree_mod.nwk")
tree = root(tree,outgroup = "NC_011740")
tree = drop.tip(tree, tip =  "NC_011740")
for (i in c(21,43,49)){
  tree = drop.tip(tree, tip =  as.character(i))
}
## get the distance between every two nodes
distance_matrix <-cophenetic(tree)
distance_matrix.m = melt_dist(distance_matrix)
distance_matrix.m$label = paste(distance_matrix.m$iso1, distance_matrix.m$iso2, sep = "-")
edges$label_a = paste(edges$ClusterA, edges$ClusterB, sep = "-")
edges$label_b = paste(edges$ClusterB, edges$ClusterA, sep = "-")
edges$label = rep("", dim(edges)[1])
edges$label[which(edges$label_a %in% distance_matrix.m$label)] = edges$label_a[which(edges$label_a %in% distance_matrix.m$label)]
edges$label[which(edges$label_b %in% distance_matrix.m$label)] = edges$label_b[which(edges$label_b %in% distance_matrix.m$label)]
edges$distance =  distance_matrix.m$dist[match(edges$label, distance_matrix.m$label)]


## check if some are sharing more than expected based on size and distance
# threshold = quantile(x = distance_matrix.m$dist, probs = c(0.75))
threshold = 0.15  ## -> threshold is chosen based on figure S5B, remove within phylogroup sharing
edges = edges[-which(edges$distance < threshold),]


edges = data.frame(cluster = c(edges$ClusterA, edges$ClusterB), 
                   cluster2 = c(edges$ClusterB, edges$ClusterA),
                   count = c(edges$Count, edges$Count), stringsAsFactors = F)
edges$size = cluster_sizes$Size[match(edges$cluster, cluster_sizes$Cluster)]
edges$size2 = cluster_sizes$Size[match(edges$cluster2, cluster_sizes$Cluster)]
#edges$min = log10(pmin(edges$size, edges$size2))
edges$min = log10(edges$size)

## normalisation step

res = lm(data = edges,count~min)
edges$normalised = edges$count - res$coefficients[2]*edges$min - res$coefficients[1]
edges$normalised = edges$normalised -  min(edges$normalised)
## for supplementary -> to explain the normalisation
S6_A= ggplot(edges, aes(x = size, y = count)) + scale_x_log10() + geom_point()  +geom_smooth(method = "lm") + theme_bw(base_size = 12) +
  xlab("Lineage size") + ylab("Number of shared 'intermediate and rare' genes")
S6_B = ggplot(edges, aes(x = size, y = normalised)) + scale_x_log10() + geom_point() + geom_smooth(method = "lm") + theme_bw(base_size = 12)+
  ylab("Normalised number of shared 'intermediate and rare' genes")+ xlab("Lineage size") 
summary(lm(data = edges,count~min))
summary(lm(data = edges,normalised~min))

grid.arrange(S6_A + ggtitle("A"), S6_B + ggtitle("B"), ncol = 2)


edges$label = edges$cluster
edges$label[which(!edges$label %in% c(12,40))] = "All other\nlineages"

### for main text, after normalising for the size of the lineage, lineages 12 and 40 are very different from the rest
edges$label = factor(edges$label, c("All other\nlineages",12,40))
A = ggplot(edges, aes(x=label, y = normalised)) + geom_boxplot(width = 0.8, outlier.size = 0.5,fill = "#eeeeee") +
  ylab("Normalised number of shared\n'intermediate and rare' genes") +
  xlab("Lineage") + theme_classic(base_size = 12)+
  geom_signif(comparisons = my_comparisons,step_increase=0.05,
              map_signif_level=T, tip_length = 0) 

wilcox.test(edges$normalised[edges$label == "12"], edges$normalised[edges$label == "All other\nlineages"])
wilcox.test(edges$normalised[edges$label == "40"], edges$normalised[edges$label == "All other\nlineages"])


pairwise.wilcox.test(edges$normalised, as.character(edges$label), p.adjust.methods = "FDR")


## for supplementary -> to show the grouping doesn't change anything
edges$cluster = factor(edges$cluster, 1:51)
ggplot(edges, aes(x=cluster, y = normalised)) +  geom_boxplot() +
  ylab("Shared intermediate and rare\ngenes with other lineages") +
  xlab("Lineage") + theme_bw(base_size = 12) 



######## PANEL B ########
## Number of rare genes in a single genome
num_genes = read.table("/Users/gh11/Submissions/pan_genome/Figures/newFigure4/genes_per_ecoli.tab", sep = "\t",
                       comment.char = "", stringsAsFactors = F, header = T)

num_genes$label = num_genes$cluster
num_genes$label[which(!num_genes$label %in% c(12,40))] = "All other\nlineages"

## for main text
num_genes$label = factor(num_genes$label, c("All other\nlineages",12,40))
  B=  ggplot(num_genes, aes(x=label, y = count)) + geom_boxplot(width = 0.8, outlier.size = 0.5, fill = "#eeeeee")+
  ylab("Lineage specific rare\ngenes in one genome") +
  xlab("Lineage") + theme_classic(base_size = 12) + scale_y_continuous(limits = c(0,150)) +
  geom_signif(comparisons = my_comparisons,step_increase=0.05,
              map_signif_level=T, y_position = c(143,150), tip_length = 0) 
## It's ok to say that I removed 19 outliers which were making the plot hard to view 

wilcox.test(num_genes$count[num_genes$label == "12"], num_genes$count[num_genes$label == "All other\nlineages"])
wilcox.test(num_genes$count[num_genes$label == "40"], num_genes$count[num_genes$label == "All other\nlineages"])

## for supplementary
num_genes$cluster = factor(num_genes$cluster, 1:51)
ggplot(num_genes, aes(x=cluster, y = count)) +  geom_boxplot( outlier.size = 0.5) +
  ylab("Lineage specific rare genes in one genome") +
  xlab("Lineage") + theme_bw(base_size = 12) + scale_y_continuous(limits = c(0,150)) 
## also here I removed 19 outliers.



######## PANEL C ########
### the size of the genomes
genome_size = read.table("/Users/gh11/Submissions/bioresource/data/F1_genome_metadata.csv",
                         sep = ",", comment.char = "", stringsAsFactors = F, header = T)
genome_size = genome_size[genome_size$PopPUNK %in% num_genes$cluster,]
genome_size$label = genome_size$PopPUNK
genome_size$label[which(!genome_size$label %in% c(12,40))] = "All other\nlineages"

## for main text
genome_size$label = factor(genome_size$label, c("All other\nlineages",12,40))
C =  ggplot(genome_size, aes(x=label, y = Length_Mbp)) + geom_boxplot(width = 0.8, outlier.size = 0.5, fill = "#eeeeee")+
  ylab("Genome length (Mbp)") +
  xlab("Lineage") + theme_classic(base_size = 12) +
  geom_signif(comparisons = my_comparisons,step_increase=0.05,
              map_signif_level=T, tip_length = 0) 



## for supplementary
genome_size$PopPUNK = factor(genome_size$PopPUNK, 1:51)
ggplot(genome_size, aes(x=PopPUNK, y = Length_Mbp)) + geom_boxplot( outlier.size = 0.5)+
  ylab("Genome length (Mbp)") +
  xlab("Lineage") + theme_bw(base_size = 12) 


grid.arrange(A + ggtitle("A") ,B + ggtitle("B"),C + ggtitle("C"), ncol = 3)





