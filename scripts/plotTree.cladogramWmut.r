#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(treeio)
library(ape)
library(ggsci)
library(ComplexHeatmap)
library(phytools)
library(phangorn)
library(grid)
library(gridExtra)
library(pheatmap)


#args = c("output/chr1-25199846-25212045-IGC/iqtree/chr1-25199846-25212045-IGC.treefile.wMut.new.nwk", "output/chr1-25199846-25212045-IGC/data/mapping_hap_SV.txt", "chr1-25199846-25212045-IGC")

#args = c("output/m1_recent_noR_noM/sim5/iqtree/sim5.treefile.wMut.new.nwk", "output/m1_recent_noR_noM/sim5/data/mapping_hap_SV.txt", "sim5")


fin_tree <- args[1]
fin_mapping_hap_SV <- args[2]
outf_prefix <- args[3]
outf_suffix <- args[4]

mapping_hap_SV <- read.delim(fin_mapping_hap_SV, sep='')

tree <- read.tree(fin_tree)

df_tree <- as.data.frame(as_tibble(tree))

num_seqs = length(tree$tip.label)


tmp_df_pop <- as.data.frame(t(apply(df_tree[1:num_seqs,], 1, function(x) unlist(strsplit(x[4], split="_")))))
tmp_df_pop$label3 <- as.numeric(as.character(tmp_df_pop[,3]))
#tmp_df_pop$label3 <- as.numeric(levels(tmp_df_pop[,3]))[tmp_df_pop[,3]]
Pop <- apply(tmp_df_pop, 1, function(x) {if (x[1] == "CMP") {return(x[1])} else{ if(x[3]==1){ return(substr(x[1],1,2))} else{return(substr(x[2],1,2))}}})

#Pop <- sapply(strsplit(df_tree[1:num_seqs,"label"], split = "_"), function(x) x[[2]])


names_Pop <- names(table(Pop))
numPop <- length(names_Pop)

mycolors_tmp <- colorRampPalette(c('#a50026','#BF9000','#92c5de','#313695'))(length(table(Pop)))
cols <- c(mycolors_tmp, "black", "#7F8080", NA,"#7F8080", NA)
names(cols) <- c(names(table(Pop)),"internal","svEvent",NA,"SV",NA)

Pop <- c(Pop, rep("internal",dim(df_tree)[1]-num_seqs))

df_tree$Pop <- Pop
df_tree$Pop <- ifelse(df_tree$Pop == "CMP", "Chimpanzee", df_tree$Pop)

df_tree$branchGroup <- ifelse(grepl("CMP",df_tree$label), 2, 1)
df_tree$branchGroup <- as.factor(df_tree$branchGroup)

getState <- function(r){
  a <- unlist(strsplit(r[4], split="_"))
  if (length(a) > 2){
    return(NA)
  }else{
    return(a[2])
  }
}

df_tree$state <- apply(df_tree, 1, function(x) getState(x))
df_tree$numState <- apply(df_tree,1, function(x) ifelse(grepl("notMut",x[7]), NA, nchar(x[7])))
df_tree$svEvent <- ifelse(df_tree$numState>1, "svEvent", NA)

df_tree$SV <- mapping_hap_SV[ match(df_tree$label, mapping_hap_SV$orig_hapID), "SV"]
df_tree$SV <- ifelse(df_tree$SV == T, "SV",NA)
df_tree$SVcode <- ifelse(df_tree$SV == T, as.character(2), NA)

p_title <- paste("SV",outf_prefix,sep="\n")

ggtree(tree, layout = "circular", branch.length = "none", aes(color=Pop), size=1.2)  %<+% df_tree + 
  geom_nodepoint(aes(color=svEvent, shape=svEvent), size=5) +
  geom_tippoint(aes(color=SV, shape=SV), size=5) + 
  scale_shape_manual(values=c(17,20,NA)) + 
  scale_color_manual(values=cols, breaks=names_Pop) +
  #scale_color_manual(values=c("#7F8080","#2081f9", "#7F8080","#A51E2B")) + 
  guides(fill="none", size="none", 
         shape=guide_legend(title=p_title, nrow = 2),
         color=guide_legend(title="Population", nrow=2)) + 
  theme(legend.position="top", legend.text = element_text(size=16), legend.title = element_text(size=16), legend.key.size = unit(1, 'cm') )

ggsave(paste(outf_prefix, ".", outf_suffix, sep=''), height = 9, width = 12, limitsize=F, useDingbats=FALSE)

