#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(gmodels)
library(irlba)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
#library(missMDA)

t_col <- function(color, percent = 70, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  lt.col = NULL
  for (i in 1:length(color)){
    ## Get RGB values for named color
    rgb.val <- col2rgb(color[i])
    
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
    lt.col <- c(lt.col, t.col)
  }
  ## Save the color
  invisible(lt.col)
}


##### Input data

#args = c("output/chr1-25199846-25212045-IGC/data/chr1-25199846-25212045-IGC.vcf.hap.tped", "output/chr1-25199846-25212045-IGC/data/chr1-25199846-25212045-IGC.vcf.hap.tfam", "output/chr1-25199846-25212045-IGC/data/mapping_hap_SV.txt","output/chr1-25199846-25212045-IGC/plot/chr1-25199846-25212045-IGC", "CMP_CMP_0")


fin_tped = args[1]
fin_tfam = args[2]
fin_IGC = args[3]
outf_prefix = args[4]
OutGroup = args[5]
outf_dir = args[6]

orig1_tped <- read.delim(fin_tped, header=F, sep='', na.strings = ".")
# get the haplotypes by rows
orig_tped <-  t(orig1_tped[,-c(1,2,3,4)])

orig_tfam <- read.delim(fin_tfam, header=F, sep='')
l_hapIDs_tfam <- paste(rep(orig_tfam$V1, each=2), rep(c(1,2), dim(orig_tfam)[1]), sep="_")
rownames(orig_tped) <- l_hapIDs_tfam
df_l_hapIDs_tfam <- data.frame(V1=l_hapIDs_tfam, V2=l_hapIDs_tfam)


orig_IGC <- read.delim(fin_IGC, sep='')
orig_IGC <- merge(orig_IGC, df_l_hapIDs_tfam, by.x="orig_hapID", by.y="V1",all.y = T) 
missingIGC_hapIDs <- which(is.na(orig_IGC$SV))
index_SVhap = which(orig_IGC$SV==T)

if (length(missingIGC_hapIDs) != 0){
  tped <- orig_tped[ -missingIGC_hapIDs,]
  IGC <- orig_IGC[ -missingIGC_hapIDs,]
}else{
  tped <- orig_tped
  IGC <- orig_IGC
}

# only keep seqs listed in df_pass_hapIDs (aka has seq > 50% of its length)
#tped <- tped[ rownames(tped) %in% df_pass_hapIDs$V1,]
#IGC <- IGC[ IGC$orig_hapID %in% df_pass_hapIDs$V1,]


#IGC$Pop <- apply(IGC, 1, function(x) unlist(strsplit(as.character(x[4]), split ="_"))[1])


#IGC$Pop <- factor(IGC$Pop, levels=c("AFR", "AMR", "EAS", "EUR","SAS"))
mycolors_tmp <- colorRampPalette(c('#a50026','#BF9000','#92c5de','#313695'))(5)

col <- colorRampPalette(mycolors_tmp)(length(unique(IGC$SV)))[factor(IGC$SV)]

#IGC$indivPop <- apply(IGC, 1, function(x) unlist(strsplit(as.character(x[4]), split ="_"))[2])

new_col <- t_col(col, percent = 50)
names(new_col) <- as.character(IGC$SV)
names(col) <- as.character(IGC$SV)
my_pch <- rep(20, dim(IGC)[1])
my_pch[index_SVhap] <- 17
my_pch_size <- rep(2, dim(IGC)[1])
my_pch_size[index_SVhap] <- 2.5

checkUniqEntry <- function(column){
  l_entry <- unique(column)
  tmp <- unique(l_entry)
  num_nonNA_entry <- length(tmp[!is.na(tmp)])
  return(num_nonNA_entry)
}

# subsetting SNVs that are variable
uniquelength <- sapply(as.data.frame(tped), checkUniqEntry)
new_tped <- as.data.frame(subset(tped, select=uniquelength>1))
upper_k = min(dim(new_tped)) - 1

#if (any(is.na(new_tped))){
#  nb <- estim_ncpPCA(new_tped, ncp.min=0, ncp.max=30)
#  new_tped.impute <- imputePCA(new_tped, ncp=nb$ncp)
#  tped.pca <- prcomp(new_tped.impute$completeObs)
#  tped.pca$rotation

#} else{
#  tped.pca <- prcomp_irlba(new_tped, n=upper_k, approx=FALSE)
#}

tped.pca <- prcomp_irlba(new_tped, n=upper_k, approx=FALSE)
summary_pca <- summary(tped.pca)
summary_pca$importance    
list_propVar_PC <- apply(data.frame(ID=names(summary_pca$importance[2,]), PC=round(summary_pca$importance[2,] * 100, digits=1)), 1, function(x) sprintf("%s(%s%%)",x[1],x[2]))

eigenvec <- tped.pca$x
project.pca <- eigenvec


# k-means
numK <- min(which(cumsum(summary_pca$importance[2,])>0.99)) # find num of PCs explains >85% variance
mydata <- project.pca[,1:numK]
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))

for (i in 2:(numK)) wss[i] <- sum(kmeans(mydata, centers=i,  nstart=25, iter.max=1000)$withinss)

#pdf("kmeans_screeplot.pdf",height = 4, width = 6)
print(paste(outf_dir, "/", outf_prefix, "_kmeans_screeplot.pdf", sep=""))
pdf(paste(outf_dir, "/", outf_prefix, "_kmeans_screeplot.pdf", sep=""),height = 4, width = 6, useDingbats=FALSE)
plot(1:(numK), wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
dev.off()



#### plot pairwise PCs ####
if(numK > 2){
  df_pairwisePC <- combn(1:3, m=2)
}else{
  df_pairwisePC <- combn(1:numK, m=2)
}

for (ii in 1:dim(df_pairwisePC)[2]){
  i=df_pairwisePC[1,ii]; j=df_pairwisePC[2,ii]; pdfFile=sprintf("%s/%s.PC%ivsPC%i.pdf", outf_dir, outf_prefix, i, j)
  pdf(pdfFile, height = 9, width = 9, useDingbats=FALSE)
  plot(project.pca[,i], project.pca[,j], type="n", main='', adj=0.5, xlab=list_propVar_PC[i], ylab=list_propVar_PC[j], bty='n')
  points(project.pca[,i], project.pca[,j], col=new_col, pch=my_pch, cex=my_pch_size)
  dev.off()
}



l_OutGroup <- unlist(strsplit(OutGroup,","))

#l_OutGroup <- l_OutGroup[ l_OutGroup %in% df_pass_hapIDs$V1 ]

dateFile = NULL
for (i in 1:length(l_OutGroup)){
  dateFile <- rbind(dateFile, data.frame(V1=l_OutGroup[i], V2=IGC[! IGC$orig_hapID %in% l_OutGroup, "orig_hapID"], V3=-6))
}

dateFile$V0 = paste(dateFile$V1, dateFile$V2, sep=",")

outf_dateFile = paste(outf_dir, "/../data/dateFile.txt", sep="")
write.table(dateFile[,c(4,3)], outf_dateFile, quote = F, col.names=F, row.names = F, sep='\t')



q()

### output some files for later usage

orig_tfam <- read.delim(fin_tfam, header=F, sep='')

orig_tfam$V7 <- apply(orig_tfam, 1, function(x) paste(unlist(strsplit(as.character(x),"_"))[2:3],collapse="_"))


poplabels <- data.frame(sample=orig_tfam$V1, population="POP_0", group="GRP0", sex=NA)
write.table(poplabels, "../data/poplabels.relate", quote = F, row.names = F, sep='\t')

sample_relate1 <- data.frame(ID_1=0, ID_2=0, missing=0)
sample_relate <- rbind(sample_relate1,	data.frame(ID_1=orig_tfam$V1, ID_2=orig_tfam$V1, missing=0))
write.table(sample_relate, "../data/sample.relate", quote = F, row.names = F, sep='\t')


#orig_IGC$V3 <- gsub("hap","", orig_IGC$V1)

