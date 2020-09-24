library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)

setwd("/Users/gh11/poppunk_pangenome/5_classify_genes/")

colours = read.table("colours_v2.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)
gene_classes = read.table("../4_pairwise_roary/231019_corrected//gene_classes.csv", sep = ",", comment.char = "", header = T, 
                          stringsAsFactors = F, quote = "", row.names = 1)



## find genes which are only ever core
gene_classes$gene_class = rep("Varied", dim(gene_classes)[1])
gene_classes$gene_class[which(gene_classes$core == gene_classes$total_presence)] = "Core"
gene_classes$gene_class[which(gene_classes$inter == gene_classes$total_presence)] = "Intermediate"
gene_classes$gene_class[which(gene_classes$rare == gene_classes$total_presence)] = "Rare"

gene_classes$label = paste(gene_classes$total_presence, "/47",sep ="")
gene_classes$label = factor(gene_classes$label, paste(1:47, "/47",sep=""))

gene_classes$fill = rep("Core, intermediate and rare", dim(gene_classes)[1])
gene_classes$fill[which(gene_classes$gene_class == "Core" & !(gene_classes$total_presence %in% c(1,47)))] = "Multi-cluster core"
gene_classes$fill[which(gene_classes$gene_class == "Core" & gene_classes$total_presence == 47)] = "Population core"
gene_classes$fill[which(gene_classes$gene_class == "Core" & gene_classes$total_presence == 1)] = "Cluster specific core"
gene_classes$fill[which(gene_classes$gene_class == "Intermediate" & gene_classes$total_presence == 1)] = "Cluster specific intermediate"
gene_classes$fill[which(gene_classes$gene_class == "Intermediate" & gene_classes$total_presence != 1)] = "Multi-cluster intermediate"
gene_classes$fill[which(gene_classes$gene_class == "Rare" & gene_classes$total_presence == 1)] = "Cluster specific rare"
gene_classes$fill[which(gene_classes$gene_class == "Rare" & gene_classes$total_presence != 1)] = "Multi-cluster rare"
gene_classes$fill[which(gene_classes$gene_class == "Varied" & gene_classes$rare == 0)] = "Core and intermediate"
gene_classes$fill[which(gene_classes$gene_class == "Varied" & gene_classes$inter == 0)] = "Core and rare"
gene_classes$fill[which(gene_classes$gene_class == "Varied" & gene_classes$core == 0)] = "Intermediate and rare"

gene_classes$fill = factor(gene_classes$fill, colours$Class)

B = ggplot(gene_classes, aes(x = label, fill = fill)) + geom_bar(color = "black", lwd = 0.2) +
  facet_grid(gene_class~., scales = "free") + 
  xlab("") + theme_classic(base_size = 12) + ylab("Number of genes") +
  scale_fill_manual(values = colours$Colour, guide = F)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Number of lineages in which gene is present") 

B

counts = data.frame(table(gene_classes$fill))
colours$counts = counts$Freq[match(colours$Class, counts$Var1)]

colours$Main.Class=factor(colours$Main.Class, unique(colours$Main.Class))
colours$Class=factor(colours$Class, colours$Class)

A = ggplot(colours, aes(x = Main.Class, y = counts, fill = Class)) + geom_bar(stat = "identity", color = "black", lwd = 0.2) +
  scale_fill_manual(values = colours$Colour, "Distribution class", labels = colours$New_name) + theme_classic(base_size = 12) +
  scale_y_continuous( breaks = seq(0, 35000, by = 5000), expand = c(0,0,0.1,0)) +
  ylab("Genes") + xlab("Distribution class") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

A

legend = as_ggplot(get_legend(A))
legend
A = A + theme(legend.position = "None")

## draw a plot where the X axis is the mean frequency when present and the Y axis is the number of PoPPUNK cluster, colour by category
freqs = read.table("../4_pairwise_roary/231019_corrected//freqs.csv", header = T,row.names = 1,
                   stringsAsFactors = F, comment.char = "", quote = "", sep =",")
freqs = freqs[,-which(colnames(freqs) == "X50")]
freqs = freqs[match(rownames(gene_classes), row.names(freqs)),]

mean_without_zeros <- function(x) {
  x = as.numeric(x)
  zeros = which(x == 0)
  if (length(zeros > 0 )) {
    return (mean(x[-zeros]))
  }
  return(mean(x))
}
gene_classes$means =  apply(X = freqs, 1, FUN = mean_without_zeros)

## plot showing how these gene categories are on a 2d plot

gene_classes$total_presence_jittered = jitter(gene_classes$total_presence, amount = 0.1)
D =  ggplot(gene_classes, aes(y = means, x = total_presence_jittered, fill = fill)) + geom_hex(bins = 60)+
  #geom_jitter(height = 0, alpha = 0.3, width = 0.1, pch = 21, color = "black", stroke = 0.1, size = 2) +
  scale_fill_manual(values = colours$Colour, name = "Category", guide = F) + theme_classic(base_size = 12) +
  ylab("Mean frequency\nwhen present") + xlab("Number of lineages in\nwhich gene is present") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("D")


## decided to add pieplot to here
summarised_typical = read.table("typical_ecoli_median.csv", sep = "\t", header = T, stringsAsFactors = F, comment.char = "")
summarised_typical$Group.1 = factor(summarised_typical$Group.1, rev(colours$Class))
summarised_typical$x = round(summarised_typical$x , digits = 0)
summarised_typical <- summarised_typical %>% mutate(pos = cumsum(x)- x/2)
ggplot(summarised_typical,aes(fill = Group.1, x= "", y =x))+ geom_bar(stat = "identity", color = "black", lwd = 0.1) + 
  coord_polar("y", start=0) +
  scale_fill_manual(values = rev(colours$Colour), guide = F) +  theme_minimal()+ ggtitle("C")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank()) +
  theme(axis.text.x=element_blank())
  #geom_text(aes(x= "", y=pos, label = x), size=4)


### 950W, 1000H
grid.arrange(A + ggtitle("A"),  B + ggtitle("B"),legend, C + ggtitle("C"), layout_matrix = rbind(c(1,1,2,2,2,2),
                                                            c(1,1,2,2,2,2),
                                                            c(4,4,2,2,2,2),
                                                            c(4,4,2,2,2,2),
                                                            c(3,3,NA,NA,NA,NA),
                                                            c(3,3,NA,NA,NA,NA)))


# last time:
# ## 1100 x 630
# grid.arrange(A + ggtitle("A"),B + ggtitle("B"),ggplot() + theme_void(), D, legend,layout_matrix = rbind(c(1,1,2,2,2,2,5,5),
#                                            c(1,1,2,2,2,2,5,5),
#                                            c(3,3,2,2,2,2,4,4),
#                                            c(3,3,2,2,2,2,4,4)))
## save as H315*W275
C + ggtitle(C)
 ## need to save the new classification to a file
#write.table(gene_classes, "classification_v2.csv", quote = F, row.names = T, col.names = T, sep = "\t")


