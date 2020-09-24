library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(data.table)
library(ggpubr)
library(dplyr)

setwd("/Users/gh11/poppunk_pangenome/5_classify_genes/")

complete_presence_absence = fread("../4_pairwise_roary/231019_corrected//complete_presence_absence.csv", sep = ",", header = T, stringsAsFactors = F)
classification = read.table("classification_v2.csv", sep = "\t", header = T, stringsAsFactors = F, comment.char = "", quote = "")



strains = colnames(complete_presence_absence)[-1]
genes = unlist(complete_presence_absence[,1][-1])
clusters = unlist(complete_presence_absence[1,])[-1]

## add the phylogroup information from the graphics file
graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter4/figures/cluster_graphics.csv", sep = ",",
                      header = T, comment.char = "")

### Section to generate the random sampling
num_clusters = length(unique(clusters))
num_classes = length(unique(classification$fill))

genes_per_ecoli = data.frame( cluster = rep("", num_classes*num_clusters-1),
                              desc = rep("", num_clusters*num_classes-1),
                              mean = rep(0, num_clusters*num_classes-1),
                              sd = rep(0, num_clusters*num_classes-1), stringsAsFactors = F)

genes_per_ecoli_2 = data.frame(cluster = character(0),
                               class = character(0),
                               count= numeric(0), stringsAsFactors = F)

row_index = 0
for (curr in unique(clusters)) {
  print(curr)
  samples = which(clusters == curr)
  samples = samples + 1
  curr_presence_absence = data.frame(complete_presence_absence[,..samples])
  curr_presence_absence = curr_presence_absence[-1,]
  for (gene_class in unique(classification$fill)) {
    print(gene_class)
    gene_class_presence_absence = curr_presence_absence[which(genes %in% classification$gene[which(classification$fill == gene_class)]),]
    num_genes = colSums(gene_class_presence_absence)
    genes_per_ecoli_2 = rbind(genes_per_ecoli_2,
                              data.frame(cluster = rep(curr, length(num_genes)),
                                         class = rep(gene_class, length(num_genes)),
                                         count = num_genes, stringsAsFactors = F))
    
    # genes_per_ecoli$cluster[row_index] = curr
    # genes_per_ecoli$desc[row_index] = gene_class
    # genes_per_ecoli$mean[row_index] = mean(num_genes)
    # genes_per_ecoli$sd[row_index] = sd(num_genes)
    # row_index = row_index+1
  }
}  
## save to a file??
genes_per_ecoli$phylogroup = graphics$Phylogroup[match(genes_per_ecoli$cluster, graphics$Cluster)]


mean_per_lineage = aggregate(x = genes_per_ecoli_2$count, by = list(genes_per_ecoli_2$cluster, genes_per_ecoli_2$class), median)

## to measure fraction in lineages 12 and 40:
lineage_12 = mean_per_lineage[mean_per_lineage$Group.1 == 40,]
lineage_12$x[lineage_12$Group.2 == "Cluster specific rare"] / sum(lineage_12$x)


mean_all = aggregate(mean_per_lineage$x, by = list(mean_per_lineage$Group.2), median)
write.table(mean_all, "/Users/gh11/poppunk_pangenome/5_classify_genes/typical_ecoli_median.csv",
            col.names = T, row.names = F, sep = "\t", quote = F)

write.table(genes_per_ecoli_2, file = "/Users/gh11/Submissions/pan_genome/Figures/newFigure4/genes_per_ecoli.tab",
            sep = "\t", col.names = T, row.names = F, quote = F)





## check if a PopPUNK cluster is an outlier
is_outlier <- function(x) {
  ret_val = rep(0, length(x))
  ret_val[which(x < quantile(x, 0.25) - 1.5 * IQR(x) )] = 1.2
  ret_val[x > quantile(x, 0.75) + 1.5 * IQR(x)] = -0.5
  return(ret_val)
}

genes_per_ecoli = genes_per_ecoli %>%
  group_by(desc, phylogroup) %>%
  mutate(vjust = is_outlier(mean))

genes_per_ecoli$outlier = genes_per_ecoli$cluster 
genes_per_ecoli$outlier[which(genes_per_ecoli$vjust == 0)] = NA

genes_per_ecoli$cluster = factor(genes_per_ecoli$cluster, as.character(unique(unlist(genes_per_ecoli$cluster))))

colours = read.table("colours_v2.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)

genes_per_ecoli$desc = factor(summarised_all$desc, colours$Class)
genes_per_ecoli$cluster = as.character(unlist(genes_per_ecoli$cluster))


# ggplot(summarised_all, aes(x = cluster, y = mean, fill = desc)) + geom_bar(stat = "identity", color= "black", size = 0.2) +
#   theme_classic(base_size = 12) + scale_fill_manual(values = colours$Colour, guide = F)+
#   facet_grid(. ~ phylogroup, scales='free',switch = "x", space = "free_x")


genes_per_ecoli$desc = factor(genes_per_ecoli$desc, colours$Class)
genes_per_ecoli$phylogroup = factor(genes_per_ecoli$phylogroup, c("B1","C","A","E","D","F","B2","U"))


## If I only want to include some of the gene classes
# summarised_filtered = summarised_all[which(summarised_all$desc %in% c("Real core","Missing in one","25-39 specific","10-24 specific","Core and specific",
#                                                               "Intermediate and specific", "Multicluster rare","Rare and specific","2-9 varied","Intermediate and rare",
#                                                               "Secondary (specific)")),]

summarised_typical = aggregate(by = list(genes_per_ecoli$desc), x = genes_per_ecoli$mean, mean)
summarised_typical$sd = aggregate(by = list(genes_per_ecoli$desc), x = genes_per_ecoli$mean, sd)$x
summarised_typical = cbind(summarised_typical, col = colours$Colour[match(summarised_typical$Group.1, colours$Class)])
summarised_typical$Group.1 = factor(summarised_typical$Group.1, rev(summarised_typical$Group.1))
summarised_typical$col = as.character(summarised_typical$col)
summarised_typical$Main = colours$Main.Class[match(summarised_typical$Group.1, colours$Class)]
sum(summarised_typical$x)
summarised_typical$min = (summarised_typical$x - summarised_typical$sd)
summarised_typical$max = summarised_typical$x + summarised_typical$sd
summarised_typical$min_percent =summarised_typical$min/sum(summarised_typical$x)*100
summarised_typical$max_percent =summarised_typical$max/sum(summarised_typical$x)*100


ggplot(summarised_typical,aes(fill = summarised_typical$Group.1, x= "", y = summarised_typical$x))+ geom_bar(stat = "identity", color = "black", lwd = 0.1) + 
  coord_polar("y", start=0) +
  scale_fill_manual(values = rev(summarised_typical$col), guide = F) +  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank()) +
  theme(axis.text.x=element_blank()) 



genes_per_ecoli$Class = colours$Main.Class[match(as.character(genes_per_ecoli$desc), colours$Class)]
plot_list = list()
plot_list[["A"]] = A + ggtitle("A")

write.table(genes_per_ecoli, file = "genes_per_ecoli.csv", sep = "\t",
            col.names = T, row.names = F,quote = F)
write.table(summarised_typical, file = "typical_ecoli.csv", sep = "\t",
            col.names = T, row.names = F,quote = F)

### DONE
i = 2
for (curr_class in colours$Class) {
  curr = summarised_all[summarised_all$desc == curr_class,]
  curr_min = max(0,summarised_typical$min[summarised_typical$Group.1 == curr_class])
  curr_max = summarised_typical$max[summarised_typical$Group.1 == curr_class]
  curr_mean = summarised_typical$x[summarised_typical$Group.1 == curr_class]
  curr$name = curr_class
  text_colour = "white"
  if (curr_class %in% c("Intermediate and rare","Core, intermediate and rare", "Core and rare", "Cluster specific rare", "Cluster specific intermediate")) {
    text_colour = "black"
  }
  max_colour = colours$Colour[colours$Class == curr_class]
  title = chartr("123456789", "ABCDEFGHI",i)
  if (i == 10) {
    title = "J"
  } else if (i == 11) { 
    title = "K"
  }else if  (i == 12) {
    title ="L"
  }
  p = ggplot(curr, aes( x = phylogroup, y = mean, color = phylogroup, label = outlier, vjust = vjust)) + geom_boxplot() +   geom_text(na.rm = TRUE)+
    # facet_grid(~desc, scales = "free_y",  space = "free",switch = "x",  labeller = label_wrap_gen(width=5)) + 
    scale_color_brewer(palette = "Dark2", guide = F) + 
    theme_bw(base_size = 12)+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
    xlab("") + ylab("Genes") +scale_y_continuous(expand = c(0.2,0,0.2,0)) +
    scale_x_discrete(expand = c(0.1,0.1)) + ggtitle(curr_class) +
    annotate("rect", xmin = 0, xmax = 9, ymin = curr_min, ymax = curr_max, 
             alpha = .2) +geom_hline(yintercept = curr_mean) + facet_grid(. ~ name) +
    theme(strip.background = element_rect(fill=max_colour),
          strip.text = element_text(size=12, colour=text_colour)) + ggtitle(title)
  plot_list[[curr_class]] = p
  i = i + 1
}




## arrange nicely on a grid

do.call("grid.arrange", c(plot_list, ncol=4))



