library(ape)
library(ggtree)
library(gridExtra)
library(ggplot2)
library(phytools)
library(reshape2)
library(aplot)
library(ggsci)

## make 4 random phylogenies with different number of genomes.
make_tree <- function(n, i){
  t = rtree(n=n)
  t = midpoint.root(t)
  #  t$tip.label = paste(i, "(", t$tip.label, ")",sep = "")
  plot(t)
  return(t)
}


t1 = make_tree(10,1)
t2 = make_tree(8,2)
t3 = make_tree(12,3)
t4 = make_tree(6,4)

t_all = make_tree(4,"All")
t_all$tip.label = c("Lineage 1\n(n=10)","Lineage 2\n(n=8)","Lineage 3\n(n=12)","Lineage 4\n(n=6)")

grid.arrange(ggtree(t1) + ggtitle("Lineage 1\n(n=10)"),
             ggtree(t2) + ggtitle("Lineage 2\n(n=8)"),
             ggtree(t3) + ggtitle("Lineage 3\n(n=12)"),
             ggtree(t4) + ggtitle("Lineage 4\n(n=6)"),
             ggtree(t_all)+ ggtitle("Collection tree") + geom_tiplab()+xlim(NA, 2),
             layout_matrix = rbind(c(1,2,3,4,5,5),
                                   c(NA,NA,NA,NA,5,5)))



## presence and absence in one lineage
df = data.frame(tip = t1$tip.label, core = rep(1, 10), 
                inter = sample(x = c(0,1), size = 10, replace = T),
                rare = c(rep(0,4),1,rep(0,5)), stringsAsFactors = F)
df = melt(df, id.vars = "tip")

p <- ggtree(t1) + ggtitle("Lineage 1\n(n=10)")

p2 = ggplot(df, aes(x=variable, y=tip)) + 
  geom_tile(aes(fill=value), colour = "black") + scale_fill_gradient(low= "white", high = "black") + 
  theme_tree2() + theme(text = element_text(size = 12)) +
  scale_x_discrete(labels = c("Core\ngene","Intermediate\ngene","Rare\ngene"))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 %>% insert_left(p)
## save 400*400

## frequency across all lineages
df2 = data.frame(tip = t_all$tip.label, 
                collection_core = c(1,1,1,1),
                multi_cluster_core = c(0,0,1,1),
                cluster_specific_core = c(0,0,1,0),
                multi_cluster_inter = c(0,0.5,0.5,0),
                cluster_specific_inter = c(0,0.5,0,0),
                multi_cluster_rare = c(0,0.1,0.1,0),
                cluster_specific_rare = c(0,0,0.1,0),
                core_and_intermediate = c(0.5,0.5,1,1),
                core_intermediate_rare = c(0.1,0.5,1,0),
                core_and_rare = c(0,0.1,1,0),
                intermediate_and_rare = c(0,0.1,0.1,0.5),
                stringsAsFactors = F)
df2 = melt(df2, id.vars = "tip")
df2$Class = rep("Varied",dim(df2)[1])
df2$Class[1:12] = "Core"
df2$Class[13:20] = "Intermediate"
df2$Class[21:28] = "Rare"


convert= data.frame(orig = unique(df2$variable), labels = c("Collection\ncore",
                                             "Multi-lineage\ncore",
                                             "Lineage specific\ncore",
                                             "Multi-lineage\nintermediate",
                                             "Lineage specific\nintermediate",
                                             "Multi-lineage\nrare",
                                             "Lineage specific\nrare",
                                             "Core and\nintermediate",
                                             "Core, intermediate\nand rare",
                                             "Core and rare",
                                             "Intermediate\nand rare"))
df2$variable2 = convert$labels[match(df2$variable, convert$orig)]
df2$variable2 = factor(df2$variable2, convert$labels)
  


p3 = ggtree(t_all) + ggtitle("Collection tree") +   
  geom_tiplab(size=4, align=TRUE, linesize=.5) +xlim(NA, 1.5)

p4= ggplot(df2, aes(x=variable2, y=tip)) + 
  geom_tile(aes(fill=value), na.rm = T, colour = "black") + scale_fill_gradient(low = "white", high= "black")+
  theme_tree2() + theme(text = element_text(size = 12)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   facet_grid(.~Class, switch = "x", scales = "free", space = "free") 
p4

p4 %>% insert_left(p3)
## save 1100*420
