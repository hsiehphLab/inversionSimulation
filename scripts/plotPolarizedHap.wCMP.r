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

mycolors_tmp <- colorRampPalette(c('#a50026','#ffeda0','#92c5de','#313695'))(5)
cols <- c("steelblue", mycolors_tmp, "black")
#names(cols) <- c("Chimpanzee","AFR","AMR","EAS","EUR","SAS","internal")
names(cols) <- c("b","AFR","AMR","EAS","EUR","SAS","a")
#cols <- mycolors_tmp
#names(cols) <- c("AFR","AMR","EAS","EUR","SAS")

#args = c("output/chr1-25199846-25212045-IGC/data/chr1-25199846-25212045-IGC.wAA.vcf.tped","output/chr1-25199846-25212045-IGC/data/chr1-25199846-25212045-IGC.wAA.vcf.tfam","output/chr1-25199846-25212045-IGC/data/mapping_hap_SV.txt","output/chr1-25199846-25212045-IGC/data/chr1-25199846-25212045-IGC")


fin_tped <- args[1]
fin_hapID <- args[2]
SV_mapping <- args[3]
outf_prefix <- args[4]

if (length(args) == 5){
    seq2beRM = args[5]
    l_seq2beRM <- unlist(strsplit(seq2beRM,","))
}else{
    l_seq2beRM = NA
}


orig1_tped <- read.delim(fin_tped, header=F, sep='', na.strings = ".")
orig_tped <-  t(as.matrix(orig1_tped[,-c(1,2,3,4)]))
orig_tfam <- read.delim(fin_hapID, header=F, sep='')

#orig_tfam$V7 <- apply(orig_tfam, 1, function(x) paste(unlist(strsplit(as.character(x),"_"))[2:3],collapse="_"))



## create a mapping between hapID and INV status
tmp_hap1 <- rep(paste(orig_tfam$V2, "1", sep='_'), each=2)
tmp_hap2 <- rep(paste(orig_tfam$V2, "2", sep='_'), each=2)

tmp_hap3 <- rep(paste(orig_tfam$V2, "1", sep='_'), each=2)
tmp_hap4 <- rep(paste(orig_tfam$V2, "2", sep='_'), each=2)

num_sample <- dim(orig_tfam)[1]
hapID <- rep("", num_sample)
hapID[ seq(1, num_sample*2, 2)] <- tmp_hap1[ seq(1, num_sample*2, 2)]
hapID[ seq(2, num_sample*2, 2)] <- tmp_hap2[ seq(2, num_sample*2, 2)]

orig_hapID <- rep("", num_sample)
orig_hapID[ seq(1, num_sample*2, 2)] <- tmp_hap3[ seq(1, num_sample*2, 2)]
orig_hapID[ seq(2, num_sample*2, 2)] <- tmp_hap4[ seq(2, num_sample*2, 2)]

mapping_hap_SV <- read.delim(SV_mapping, sep='')

missingRows = data.frame(hapID=orig_hapID[! orig_hapID %in% mapping_hap_SV$orig_hapID], SV=NA, orig_hapID=orig_hapID[! orig_hapID %in% mapping_hap_SV$orig_hapID])

mapping_hap_SV <- rbind(mapping_hap_SV, missingRows)

mapping_hap_SV <- mapping_hap_SV[ mapping_hap_SV$orig_hapID %in% orig_hapID,]



# reorder the rows/hapIDs of mapping_hap_SV so that it follows the input tfam
mapping_hap_SV <- mapping_hap_SV[ match(orig_hapID, mapping_hap_SV$orig_hapID), ]

row.names(orig_tped) <- orig_hapID


annoRow = data.frame(SV=ifelse(mapping_hap_SV[,2]==T, 1, 0))
row.names(annoRow) <- mapping_hap_SV$orig_hapID

# get the unwanted seq IDs
l_seq2beRM <- unlist(lapply(l_seq2beRM, function(x) grep(x, rownames(orig_tped), value = TRUE)))


#run_pheatmap <- function(TPED, list_seq2beRM, permitMissingFrac){
run_pheatmap <- function(TPED, list_seq2beRM, permitMissingFrac){
  TPED = orig_tped
  list_seq2beRM = l_seq2beRM
  permitMissingFrac = 1
  new_tped <- TPED[ !rownames(TPED) %in% list_seq2beRM, ]

  # remove rows with all zeros and ones
  new_tped <- new_tped[,!(apply(new_tped, 2, function(x) all(new_tped[ !is.na(new_tped)]==0)) | apply(new_tped, 2, function(x) all(new_tped[ !is.na(new_tped)]==1)))]

  new_tped <- new_tped[rowSums(is.na(new_tped)) <= ncol(new_tped)*permitMissingFrac, ]
  
  p2 <- pheatmap(new_tped, cluster_cols = F, cluster_rows=T, clustering_distance_rows = "binary", clustering_method = "complete", color=c("darkblue","darkorange"), labels_col = "", labels_row = NULL, annotation_row =annoRow)
  return(p2)
}


for(i in seq(1,0,-0.1)){
  tryCatch({output <-  run_pheatmap(orig_tped, l_seq2beRM, i)
  possibleError <- NA
  break
  }, error = function(e) {
    PossibleError <- e
  }) 
}  

if(!is.na(possibleError)){
  output <- pheatmap(orig_tped, cluster_cols = F, cluster_rows=F, clustering_distance_rows = "binary", clustering_method = "complete", color=c("darkblue","darkorange"), labels_col = "", labels_row = NULL, annotation_row =annoRow)
}

outPDF = sprintf("%s.polarizedHap.wCMP.pdf", outf_prefix )
ggsave(outPDF, output, width=12, height = 9, useDingbats = F)
output
dev.off()



