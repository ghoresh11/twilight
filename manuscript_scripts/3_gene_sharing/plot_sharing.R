library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(reshape2)
library(igraph)
library(ape)
library(ggtree)
library(stringr)
library(gridExtra)
library(scatterplot3d)
library(harrietr)

setwd("/Users/gh11/poppunk_pangenome/10_gene_sharing/")


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


cluster_graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter4/figures/cluster_graphics.csv", 
                              sep = ",", comment.char = "", stringsAsFactors = F, header = T)

colours = read.table("../5_classify_genes/colours_v2.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)

### NETWORK ANALYSIS
### FUNCTIONS ### 
get_lower_tri<-function(m){
  m[upper.tri(m)] <- NA
  return(m)
}

reorder_matrix <- function(m){
  hc <- hclust(dist(m))
  m <-m[hc$order, hc$order]
  return(m)
}

plot_matrix <- function(m, legend){
  m <- reorder_matrix(m)
  lower_tri = get_lower_tri(m)
  #  lower_tri[which(lower_tri == 0)] = NA
  melted_m <- melt(lower_tri, na.rm = TRUE)
  melted_m$Var1 = factor(melted_m$Var1, melted_m$Var1[1:length(unique(melted_m$Var1))])
  melted_m$Var2 = factor(melted_m$Var2, melted_m$Var1[1:length(unique(melted_m$Var1))])
  p = ggplot(data = melted_m, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "black")+
    scale_fill_gradient( low = "white", high = "blue", space = "Lab", name = legend) +
    theme_bw(base_size = 16)+
    coord_fixed() + xlab("Cluster") + ylab("Cluster") + theme(axis.text.x = element_text(angle = 90, vjust = 1))
  return (p)
}  

swap_values_with_cols <- function(vec){
  num_cols = length(unique(vec))
  cols = brewer.pal(12, "Paired")
  vec_to_return = as.character(vec)
  i = 1
  for (value in unique(vec)){
    vec_to_return[which(vec == value)] = cols[i]
    i = i + 1
  }
  return(vec_to_return)
}


generate_one_plot <- function(column) {
  df = pca_x
  colnames(df)[which(colnames(df) == column)] = "property"
  p = ggplot(df, aes(x = PC1, y = PC2, fill = property, color = property, Label = Label)) + 
    geom_text(aes(label = Label), size=5,
              fontface = "bold",position=position_jitter(width=0.2,height=0.2, seed = 1),
              hjust = -0.3, vjust = 0.4) + 
    geom_jitter(pch = 21, size = 2, position = position_jitter(width=0.2,height=0.2,seed = 1)) +
    theme_bw(base_size = 14) + scale_fill_brewer(palette = "Paired", name = column) + scale_color_brewer(palette = "Paired", guide = F) +
    theme(legend.position = "bottom") + guides(fill=guide_legend(ncol=2))
  pdf(paste("plots_v2/",gene_class,"_",column,"_3d.pdf",sep=""))
  s3d = scatterplot3d(x = pca_x$PC1, y = pca_x$PC2, z = pca_x$PC3, axis = T, color = swap_values_with_cols(df$property), pch = 16, grid = F, 
                      xlab = "PC1", ylab = "PC2", zlab = "PC3", main = column)
  text(s3d$xyz.convert(pca_x$PC1, pca_x$PC2, pca_x$PC3), labels=pca_x$Label, pos=1,  col= swap_values_with_cols(df$property))
  dev.off()
  return(p)
}


plot_one_genetype_against_dist <- function(gene_type, title, get_legend = F) {
  ## plot the number of genes shared compared to phylogenetic distance
  all = edges[which(edges$Gene_class == gene_type),]
  
  all$Count_normalised = (all$Count - min(all$Count)) /(max(all$Count) - min(all$Count))
  all$phylo = factor(all$phylo, c("B1","C","A","E","D","B2","F","U","Different"))
  all$name = gene_type
  text_colour = "white"
  if (gene_type %in% c("Intermediate and rare","Core, intermediate and rare", "Core and rare", "Cluster specific rare", "Cluster specific intermediate")) {
    text_colour = "black"
  }
  max_colour = colours$Colour[colours$Class == gene_type]
  print(summary(lm(data = all, distance~Count)))
  
  p = ggplot(all, aes(x = distance, y = Count, fill = phylo)) + geom_point(size = 3, pch = 21, alpha = 0.7) +
    xlab("Phylogenetic distance") + ylab("Number of shared genes") +
    theme_bw(base_size = 14)  + 
    scale_fill_manual(values = c(brewer.pal(n=7,"Set2"),"brown","#d3d3d3"),name = "", drop = F) + ggtitle(title) + 
    facet_grid(. ~ name) +
    theme(strip.background = element_rect(fill=max_colour),
          strip.text = element_text(size=12, colour=text_colour)) + theme(legend.position = "bottom") + guides(fill=guide_legend(ncol=3,byrow=F))
  if (!get_legend) {p = p + theme(legend.position = "None")}
  return(p)
}


plot_one_genetype_against_min_size <-  function(gene_type, title, get_legend = F){
  curr = edges[edges$Gene_class == gene_type,]
  curr$phylo = factor(curr$phylo, c("B1","C","A","E","D","B2","F","U","Different"))
  curr$name = gene_type
  text_colour = "white"
  if (gene_type %in% c("Intermediate and rare","Core, intermediate and rare", "Core and rare", "Cluster specific rare", "Cluster specific intermediate")) {
    text_colour = "black"
  }
  max_colour = colours$Colour[colours$Class == gene_type]
  print(summary(lm(data = curr, log10(min)~Count)))
  p = ggplot(curr, aes(x = min, y = Count, fill = phylo)) + geom_point(size = 3, pch = 21, alpha = 0.9) +
    xlab("Smaller cluster size") + ylab("Shared genes") +
    theme_bw(base_size = 14) + ggtitle(gene_type) + scale_fill_manual(values = c(brewer.pal(n=7,"Set2"),"brown","#d3d3d3"), 
                                                                      drop = F)+
    ggtitle(title)+ 
    facet_grid(. ~ name) +
    theme(strip.background = element_rect(fill=max_colour),
          strip.text = element_text(size=12, colour=text_colour))
  if (!get_legend) {p = p + theme(legend.position = "None")}
  return(p)
}

### MAIN ### 
edges = read.table("num_shared_gene_v2.csv", sep = "\t", comment.char = "",
                   stringsAsFactors = F, quote = "", header = T)
clusters_md = read.table("/Users/gh11/Submissions/my_thesis/Chapter4/figures/cluster_graphics.csv", sep = ",",
                         comment.char = "", stringsAsFactors = F, quote = "", header = T)
more_md = read.csv("/Users/gh11/Submissions/my_thesis/Chapter4/figures_thesis/ecoli_cluster_summary.csv", sep = ",",
                   stringsAsFactors = F, header = T, comment.char = "")
more_md = more_md[match(clusters_md$Cluster, more_md$Cluster),]
clusters_md = cbind(clusters_md, ST = more_md$STs, MDR = more_md$MDR, Pathotype = more_md$Pathotypes)

## Get the cluster sizes
cluster_sizes = read.table("../2_dists_roary_analysis/cluster_sizes_updated.csv", sep = ",", comment.char = "", header = T,stringsAsFactors = F)
edges$sizeA = cluster_sizes$Size[match(edges$ClusterA, cluster_sizes$Cluster)]
edges$sizeB = cluster_sizes$Size[match(edges$ClusterB, cluster_sizes$Cluster)]
edges  = transform(edges, min = pmin(edges$sizeA , edges$sizeB))
edges  = transform(edges, max = pmax(edges$sizeA , edges$sizeB))

## add phylogroup info
edges$phyloA = cluster_graphics$Phylogroup[match(edges$ClusterA, cluster_graphics$Cluster)]
edges$phyloB = cluster_graphics$Phylogroup[match(edges$ClusterB, cluster_graphics$Cluster)]
edges$phylo = rep("Different", dim(edges)[1])
edges$phylo[which(edges$phyloA == edges$phyloB)] = edges$phyloA[which(edges$phyloA == edges$phyloB)]

## add distance
edges$label_a = paste(edges$ClusterA, edges$ClusterB, sep = "-")
edges$label_b = paste(edges$ClusterB, edges$ClusterA, sep = "-")
edges$label = rep("", dim(edges)[1])
edges$label[which(edges$label_a %in% distance_matrix.m$label)] = edges$label_a[which(edges$label_a %in% distance_matrix.m$label)]
edges$label[which(edges$label_b %in% distance_matrix.m$label)] = edges$label_b[which(edges$label_b %in% distance_matrix.m$label)]
edges$distance =  distance_matrix.m$dist[match(edges$label, distance_matrix.m$label)]

## Run on one gene_type
all_classes = unique(edges$Gene_class)


# # ## Build a minimum spanning tree of the gene sharing for each gene class
# for (gene_class in all_classes){
#   curr.m = edges[which(edges$Gene_class == gene_class),]
#   g = graph.data.frame(curr.m, directed=FALSE)
#   E(g)$weight = as.numeric(E(g)$Count)
#   comms = cluster_walktrap(g, steps = 3)
# 
#   curr = get.adjacency(g, attr="Count", sparse=FALSE)
#   for_mst = -curr
#   curr_mst = ape::mst(for_mst)
#   curr[which(curr_mst == 0)] = 0
# 
#   g = graph_from_adjacency_matrix(curr, mode = "undirected", weighted = T)
#   ## to colour and layout using different attributes
#   V(g)$size = rep(10,47)
#   coords = layout_nicely(g)
# 
#   #pdf(file = paste("plots_v2/",gene_class, "_msts.pdf", sep = ""), width = 15, height = 8)
#   par(mfrow = c(1,1))
#   
#   V(g)$color = swap_values_with_cols(clusters_md[,which(colnames(clusters_md) == "Phylogroup")])[match(V(g)$name, clusters_md$Cluster)]
#   plot(g, layout = coords, size =1)
#   write.graph(graph = g, file = "/Users/gh11/Submissions/my_thesis/Chapter5/prep_v2/genes_effect_pop_core/sharing_mst_core.gml",format = "gml")
#   #dev.off()
# }

# ## plot a 3d PCA plot
# for (gene_class in all_classes){
#   print(gene_class)
#   curr.m = edges[which(edges$Gene_class == gene_class),]
#   g = graph.data.frame(curr.m, directed=F)
#   E(g)$weight = as.numeric(E(g)$Count)
#   curr = get.adjacency(g, attr="Count", sparse=FALSE)
#   pca_res = prcomp(curr, scale = TRUE, center = )   
#   pca_x = pca_res$x
#   curr_clusters_md = clusters_md[match(row.names(pca_x), clusters_md$Cluster),]
#   pca_x = cbind(pca_x, curr_clusters_md)
#   pca_x$Label = row.names(pca_x)
#   p1 = generate_one_plot("Phylogroup")
#   p2 = generate_one_plot("Pathotype")
#   p3 = generate_one_plot("MDR")
#   p4 = generate_one_plot("Virulence")
#   p = arrangeGrob(p1,p2,p3, p4, nrow = 2)
#   ggsave(plot = p, filename = paste("plots_v2/",gene_class, "_pca.pdf", sep = ""), width = 10, height = 10)
# }
# 

core = plot_one_genetype_against_dist("Multi-cluster core", "", T)
intermediate = plot_one_genetype_against_dist("Multi-cluster intermediate", "")
rare = plot_one_genetype_against_dist("Multi-cluster rare", "")
core_and_inter = plot_one_genetype_against_dist("Core and intermediate","")
core_and_rare = plot_one_genetype_against_dist("Core and rare","")
core_inter_rare = plot_one_genetype_against_dist("Core, intermediate and rare","")
inter_and_rare = plot_one_genetype_against_dist("Intermediate and rare","")


## extract the legend from the first plot
legend = as_ggplot(get_legend(core))
core = core + theme(legend.position = "None")

# ## plot all og them 
# grid.arrange(core, intermediate, rare, legend, core_and_inter,core_and_rare,
#              core_inter_rare,  inter_and_rare, nrow = 2)


core2 = plot_one_genetype_against_min_size("Multi-cluster core", "A")
intermediate2 = plot_one_genetype_against_min_size("Multi-cluster intermediate", "B")
rare2 = plot_one_genetype_against_min_size("Multi-cluster rare", "C")
core_and_inter2 = plot_one_genetype_against_min_size("Core and intermediate","D")
core_and_rare2 = plot_one_genetype_against_min_size("Core and rare","E")
core_inter_rare2 = plot_one_genetype_against_min_size("Core, intermediate and rare","F")
inter_and_rare2 = plot_one_genetype_against_min_size("Intermediate and rare","G")

# grid.arrange(core2, intermediate2, rare2,legend,
#              core_and_inter2,
#              core_inter_rare2, core_and_rare2, inter_and_rare2, nrow = 2)

## check if some are sharing more than expected based on size and distance
test = edges[which(edges$Gene_class %in% c("Intermediate and rare")),]
threshold = quantile(x = distance_matrix.m$dist, probs = c(0.75))
test = test[-which(test$distance < threshold),]
test = data.frame(cluster = c(test$ClusterA, test$ClusterB), count = c(test$Count, test$Count), stringsAsFactors = F)

# create new DF and sum up all three gene classes
test = aggregate(x = test$count, by = list(test$cluster), FUN = median) 
test$size = cluster_sizes$Size[match(test$Group.1, cluster_sizes$Cluster)]


test$logged = log10(test$size)
summary(lm(data = test, formula = x~logged))
res = residuals(lm(data = test, formula = x~logged))


test$name = rep("No", dim(test)[1])
name = test$Group.1[which(res > 120 | res < -120)]
test$name[test$Group.1 %in% name] = "Yes"

 ggplot(test, aes(x = size, y = x, label = Group.1))+ geom_point()  + scale_x_log10() +
  geom_smooth(method='lm', color = "black", se = 0.5) + theme_minimal(base_size = 12) + xlab("Lineage size") + 
  ylab("Number of 'intermediate and rare' genes shared") + ggtitle("A") 
  geom_text_repel(min.segment.length = 0, data = subset(test, name == "Yes")) 
