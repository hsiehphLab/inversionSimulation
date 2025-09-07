#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(ggpubr)
library(ggsci)
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(glue)
#library(tidyverse)
#library(rstatix)   
source("scripts/R_rainclouds.r")
#source("/Users/hseihph/mnt/vol26/home/hsiehph/bin/scripts/scripts_uw/R_plotting/R_rainclouds.r")
#source("/Users/hsiehph/mnt/vol26/home/hsiehph/bin/scripts/scripts_uw/R_plotting/R_rainclouds.r")

#args = c("output/chr1-25199846-25212045-IGC/data/chr1-25199846-25212045-IGC.IBS.dat", "output/chr1-25199846-25212045-IGC/data/chr1-25199846-25212045-IGC.elemDiff_pairHaps.dat", "output/chr1-25199846-25212045-IGC/plot/chr1-25199846-25212045-IGC.pairwiseIBS.pdf", "output/chr1-25199846-25212045-IGC/plot/chr1-25199846-25212045-IGC.elemDiff_pairHaps.pdf")

#ibs <- read.delim("tmp.out")
#elemDiff_pairHaps <- read.delim("chr11-50136371-INV-206505.elemDiff_pairHaps.dat", header=F, sep='')
ibs <- read.delim(args[1])
ibs$h1_SVgeno <- ifelse(ibs$h1_SVgeno=="False", 0, 1)
ibs$h2_SVgeno <- ifelse(ibs$h2_SVgeno=="False", 0, 1)
elemDiff_pairHaps <- read.delim(args[2], header=F , sep='', na.strings = c("nan"))
outPDF1 <- args[3]
#outPDF2 <- args[4]


ibs$SVgeno <- paste(ibs$h1_SVgeno, ibs$h2_SVgeno, sep="")
ibs$SVgeno <- ifelse(ibs$SVgeno == "10", "01", ibs$SVgeno)

min_ibs = round(min(ibs$IBS), digits=1) - 0.05
max_ibs = round(max(ibs$IBS), digits=1) + 0.05
tick_breaks = seq(max(0,round(min_ibs,1)), min(round(max_ibs,1),1), 0.1)

elemDiff_pairHaps <- elemDiff_pairHaps[ order(elemDiff_pairHaps$V1, elemDiff_pairHaps$V2),]

row.names(elemDiff_pairHaps) <- paste(elemDiff_pairHaps$V3, elemDiff_pairHaps$V4,sep=",")

list_sv_anno <- paste(elemDiff_pairHaps$V1, elemDiff_pairHaps$V2, sep="")
list_sv_anno <- ifelse(list_sv_anno=="10", "01", list_sv_anno)

table_list_sv_anno <- table(list_sv_anno)

start_groups <- c(1, unname(cumsum(table_list_sv_anno))+1)[-(length(table_list_sv_anno)+1)]
end_groups <- unname(cumsum(table(list_sv_anno)))

df_m3 <- data.frame(names=c("00","01","11"), labels = c("Ref/Ref","Ref/SV","SV/SV"))

curr_df_m3 <- droplevels(df_m3[ df_m3$names %in% names(table_list_sv_anno),])
m3 <- curr_df_m3[,"names"]
m3_labels <- curr_df_m3[,"labels"]


p1 <- ggplot(ibs, aes(x=SVgeno, y=IBS, fill=SVgeno)) + coord_flip() + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .6, adjust =1, trim = T) +
  geom_point(position=position_jitter(width = .15), size = 1, alpha = 0.5, shape=21) +
  geom_boxplot(aes(x=SVgeno, y=IBS), outlier.shape = NA, alpha=0.7, width=0.1, color="black") +
  scale_color_nejm() + scale_fill_nejm() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.title=element_text(size=16), axis.text=element_text(size=16)) +
  ylab("Pairwise identity by state") + xlab("SV genotype") + guides(fill=FALSE) + 
  scale_y_continuous(breaks=tick_breaks, limits = c(min_ibs,max_ibs)) +
  scale_x_discrete(labels=m3_labels)


pdf(outPDF1, height = 6, width = 9, useDingbats=FALSE)
#plot_grid(p1, p2, ncol=1, align = "v", rel_heights = c(3,1))
p1
dev.off()

q()
# heatmap element-wise difference between pairs of haplotypes
cols = pal_jama()(3)
colnames = c("Match","Mismatch","NA")
names(cols) <- c(0,1,NA)


options(expressions=5e5)

m_elemDiff_pairHaps <- as.matrix(elemDiff_pairHaps[,-(1:4)])


if(any(is.na(m_elemDiff_pairHaps))){
  l_at = c(0,1,NA)
  l_labels = colnames
}else{
  l_at = c(0,1)
  l_labels = colnames[1:2]
}

# get ready for top/column annotation
if (length(args) == 6){
  fin_tped1_in <- args[5]
  fin_tped2_out <- args[6]
  tped1_in_POS <- read.delim(fin_tped1_in, header=F, sep='')[,4]
  tped2_out_POS <- read.delim(fin_tped2_out, header=F, sep='')[,4]
  tmp <- which(tped2_out_POS %in% tped1_in_POS)
  indices_start_end <- c(tmp[1], tail(tmp,1))
  start_groups_top <- c(1, indices_start_end[1], indices_start_end[2]+1)
  end_groups_top <- c(indices_start_end[1]-1, indices_start_end[2], dim(m_elemDiff_pairHaps)[2])
}else{
  start_groups_top <- c(1, 1, 1)
  end_groups_top <- c(0, dim(m_elemDiff_pairHaps)[2], 0)
}


list_type_top <- NULL
for (i in 1:length(start_groups_top)){
  list_type_top <- c(list_type_top, rep(letters[i], end_groups_top[i] - start_groups_top[i] +1))
}

list_col_top <- c("a" = "white", "b" = "skyblue", "c" = "white")
anno_top = HeatmapAnnotation(type = list_type_top, col = list(type=list_col_top), annotation_height=0.5, show_legend=FALSE)


# Calculate the entropy of a vector of counts or proportions
# Inputs: Vector of numbers
# Output: Entropy (in bits)
entropy <- function(p) {
    # Assumes: p is a numeric vector
    if (sum(p) == 0) {
        return(0) # Case shows up when calculating conditional
        # entropies
    }
    p <- p/sum(p) # Normalize so it sums to 1
    p <- p[p > 0] # Discard zero entries (because 0 log 0 = 0)
    H = -sum(p*log(p,base=2))
    return(H)
}


for (i in 1:length(m3)){
  row_split = c(rep(1,sum(list_sv_anno==m3[i])))
  tmp_m_elemDiff_pairHaps <- m_elemDiff_pairHaps[start_groups[i]:end_groups[i],]
  anno_right = rowAnnotation(foo=anno_block(gp=gpar(fill=pal_d3()(3)[i]), labels = m3_labels[i], labels_gp = gpar(col="white", fontsize=10)))

  if(i==1){
    rowTitle = "Haplotype pairs"
    curr_anno_top = anno_top
  } else {
    rowTitle=""
    curr_anno_top = NULL
  }
 
  if(is.null(dim(tmp_m_elemDiff_pairHaps))){
    ht_list <- Heatmap(t(as.data.frame(tmp_m_elemDiff_pairHaps)), cluster_rows = F, cluster_columns = F, name="PairwiseComparison", col = cols, show_heatmap_legend = F,  show_column_names = F, show_row_names = F, column_title = "SNV", row_title = rowTitle, right_annotation = anno_right, top_annotation = curr_anno_top, height = 4)
  }else{
      df_tmp_hap = data.frame(ID = row.names(tmp_m_elemDiff_pairHaps), hap=apply(tmp_m_elemDiff_pairHaps,1,function(x) glue_collapse(x)))
      df_tmp_hap_sample200 <- df_tmp_hap %>% group_by(hap) %>% sample_n(if(n() < 500) n() else 500)
      tmp_m_elemDiff_pairHaps <- tmp_m_elemDiff_pairHaps[ row.names(tmp_m_elemDiff_pairHaps) %in% df_tmp_hap_sample200$ID,]
     uniqHap = dim(unique(tmp_m_elemDiff_pairHaps))[1]
      if(uniqHap < 6){
        ht_list <- Heatmap(tmp_m_elemDiff_pairHaps, row_order = order(tmp_m_elemDiff_pairHaps[,which.max(apply(tmp_m_elemDiff_pairHaps, 2, function(x) var(x)))]), cluster_rows = F, cluster_columns = F, name="PairwiseComparison", col =  cols, heatmap_legend_param = list(at = l_at, labels = l_labels), show_column_names = F, show_row_names = F, column_title = "SNV", row_title = rowTitle, right_annotation = anno_right, top_annotation = curr_anno_top, height = 4)
      }else{
        ht_list <- Heatmap(tmp_m_elemDiff_pairHaps, cluster_rows = T, cluster_columns = F, name="PairwiseComparison", col =  cols, heatmap_legend_param = list(at = l_at, labels = l_labels), show_column_names = F, show_row_names = F, column_title = "SNV", row_title = rowTitle, right_annotation = anno_right, top_annotation = curr_anno_top, height = 4)
      }
  }
  
  if(i==1) {all_ht_list <- ht_list} else{all_ht_list <- all_ht_list %v% ht_list}
}

pdf(outPDF2, height = 6, width = 9, useDingbats=FALSE)
all_ht_list
dev.off()
