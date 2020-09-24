library(ggplot2)

setwd("/Users/gh11/gene_distances/generax/")


order = c("Collection core", "Multi-lineage core", "Core and intermediate", "Core and rare", "Core, intermediate and rare", "Intermediate and rare",
          "Multi-lineage intermediate", "Multi-lineage rare")

stats = read.table("all_stats.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)

old_to_new = read.table("/Users/gh11/Submissions/bioresource/data/scripts/change_gene_names/old_to_new.csv", sep = ",",
                        comment.char = "", stringsAsFactors = F, header = T, quote = "")
classification = read.table("/Users/gh11/poppunk_pangenome/5_classify_genes/classification_v2.csv", sep = "\t",
                        comment.char = "", stringsAsFactors = F, header = T, quote = "")
classification$new_name = old_to_new$New[match(classification$gene, old_to_new$Old)]
colours = read.table("/Users/gh11/poppunk_pangenome/5_classify_genes/colours_v2.csv", sep = ",",
                     comment.char = "", stringsAsFactors = F, header = T)

stats$class = classification$fill[match(stats$Name, classification$new_name)]
stats$count = classification$total_presence[match(stats$Name, classification$new_name)]
stats$class = colours$New_name[match(stats$class, colours$Class)]

medians = aggregate(stats$Transfer, by = list(stats$class), FUN = median)


stats$main = colours$Main.Class[match(stats$class, colours$New_name)]

stats$class = factor(stats$class, order)
stats$main = factor(stats$main, c("Core","Varied","Intermediate","Rare"))

C = ggplot(stats, aes(x = class, y = Transfer)) + geom_boxplot(fill = "#eeeeee") +
  facet_grid(.~main, scales = "free", space = "free",switch = "x") +
  theme_bw(base_size = 12) + xlab("Distribution Class") + ylab("Probability of transfer event") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pairwise.wilcox.test(as.numeric(stats$Transfer), as.character(stats$class), p.adjust.method = "fdr")


## read in the acctran results to correlate with this
gain_loss = read.table("//Users/gh11/poppunk_pangenome/11_ancestral_state_recon/all_genes.tab", sep = "\t",
                       header = T, stringsAsFactors = F, quote = "")
gain_loss$new_name = old_to_new$New[match(gain_loss$gene, old_to_new$Old)]
gain_loss$class = classification$fill[match(gain_loss$gene, classification$gene)]
gain_loss$class = colours$New_name[match(gain_loss$class, colours$Class)]
gain_loss.m = melt(gain_loss[,c(5,2,3)], id.vars = c("class"))

gain_loss.m$main = colours$Main.Class[match(gain_loss.m$class, colours$New_name)]

medians_gains = aggregate(gain_loss$gains, by = list(gain_loss$class), FUN = median)
medians_losses = aggregate(gain_loss$losses, by = list(gain_loss$class), FUN = median)


gain_loss.m$main = factor(gain_loss.m$main, c("Core","Varied","Intermediate","Rare"))
gain_loss.m$class = factor(gain_loss.m$class, order)

D = ggplot(gain_loss.m, aes(x = class, y = value, fill = variable)) + geom_boxplot() +
  facet_grid(.~main, scales = "free", space = "free",switch = "x") +
  theme_bw(base_size = 12) + xlab("Distribution Class") + ylab("Events") +
  scale_fill_brewer(palette = "Set2", name = "Event type", labels = c("Gain","Loss")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pairwise.wilcox.test(as.numeric(gain_loss$gains), as.character(gain_loss$class), p.adjust.method = "fdr")
pairwise.wilcox.test(as.numeric(gain_loss$losses), as.character(gain_loss$class), p.adjust.method = "fdr")


C



### Supplementary
gain_loss$total_presence = classification$total_presence[match(gain_loss$new_name, classification$new_name)]
gain_loss = gain_loss[which(gain_loss$class %in% c("Multi-lineage core", "Intermediate and rare")),]
gain_loss$total_presence = factor(gain_loss$total_presence, 1:47)

S_A = ggplot(gain_loss, aes(x = total_presence, y = gains, fill = class)) + geom_boxplot() +
  theme_classic(base_size = 12) + scale_fill_manual(name = "Distribution class", values = c("#dbe3de","#d7d7f5"))+
  xlab("Number of lineages in which gene is present") + ylab("Gain events")
S_B = ggplot(gain_loss, aes(x = total_presence, y = losses, fill = class)) + geom_boxplot() +
  theme_classic(base_size = 12) + scale_fill_manual(name = "Distribution class", values = c("#dbe3de","#d7d7f5"))+
  xlab("Number of lineages in which gene is present") + ylab("Loss events")

stats$total_presence = classification$total_presence[match(stats$Name, classification$new_name)]
stats = stats[which(stats$class %in% c("Multi-lineage core", "Intermediate and rare")),]
stats$total_presence = factor(stats$total_presence, 1:47)
S_C = ggplot(stats, aes(x = total_presence, y = Transfer, fill = class)) + geom_boxplot() +
  theme_classic(base_size = 12) + scale_fill_manual(name = "Distribution class", values = c("#d7d7f5","#dbe3de"))+
  xlab("Number of lineages in which gene is present") + ylab("Probability of transfer")


grid.arrange(S_A + theme(legend.position =  "None") + ggtitle("A"),
             S_B + theme(legend.position =  "None") + ggtitle("B"),
             S_C + theme(legend.position =  "None") + ggtitle("C"), ncol = 1)

as_ggplot(get_legend(S_A))
# 
# ### analyse transfers and recievers
# transfers = read.table("all_transfers.csv", sep = ",", header = T, stringsAsFactors = F, comment.char = "")
# transfers$class = classification$fill[match(transfers$Name, classification$new_name)]
# donors = transfers[transfers$Type == "donor",]
# recievers = transfers[transfers$Type == "reciever",]
