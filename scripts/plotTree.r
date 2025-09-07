#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(treeio)
library(ape)
library(ggsci)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#cols <- mycolors_tmp
#names(cols) <- c("AFR","AMR","EAS","EUR","SAS")

#args = c("output/chr1-25199846-25212045-IGC/data/chr1-25199846-25212045-IGC.vcf.hap.tped", "output/chr1-25199846-25212045-IGC/data/mapping_hap_SV.txt",  "output/chr1-25199846-25212045-IGC/iqtree/chr1-25199846-25212045-IGC.treefile","output/chr1-25199846-25212045-IGC/iqtree/chr1-25199846-25212045-IGC.timetree.nex","chr1-25199846-25212045-IGC")
#args = c("output/m1_recent_noR_noM/sim2/data/sim2.vcf.hap.tped","output/m1_recent_noR_noM/sim2/data/mapping_hap_SV.txt", "output/m1_recent_noR_noM/sim2/iqtree/sim2.treefile", "output/m1_recent_noR_noM/sim2/iqtree/sim2.timetree.nex","sim2")

#args = c("output/m50pc_young_highR_lowM/sim0/data/sim0.vcf.hap.tped", "output/m50pc_young_highR_lowM/sim0/data/mapping_hap_SV.txt",  "output/m50pc_young_highR_lowM/sim0/iqtree/sim0.treefile", "output/m50pc_young_highR_lowM/sim0/iqtree/sim0.timetree.nex",  "output/m50pc_young_highR_lowM/sim0/plot", "sim0")


fin_tped <- args[1]
fin_SV_mapping <- args[2]
fin_iqtree <- args[3]
fin_timetree <- args[4]
outf_dir <- args[5]
outf_prefix <- args[6]


orig1_tped <- read.delim(fin_tped, header=F, sep='')
orig_tped <-  t(orig1_tped[,-c(1,2,3,4)])
df_SV_mapping = read.delim(fin_SV_mapping, sep='')


# plot unrooted phylogeny
iqtree <- read.iqtree(fin_iqtree)

num_seqs <- length(iqtree@phylo$tip.label)

df_iqtree <- as.data.frame(as_tibble(iqtree))

tmp_df_pop <- as.data.frame(t(apply(df_iqtree[1:num_seqs,], 1, function(x) unlist(strsplit(x[4], split="_")))))
#tmp_df_pop$label3 <- as.numeric(levels(tmp_df_pop$label3))[tmp_df_pop$label3]
tmp_df_pop$label3 <- as.numeric(as.character(tmp_df_pop[,3]))
Pop <- apply(tmp_df_pop, 1, function(x) if (x[1] == "CMP") return(x[1]) else if(x[3]==1){ return(substr(x[1],1,2))} else{return(substr(x[2],1,2))})


#Pop <- sapply(strsplit(df_iqtree[1:num_seqs,"label"], split = "_"), function(x) x[[2]])
names_Pop <- names(table(Pop))
if("CMP" %in% names_Pop){
  num_Pop <- length(names_Pop) - 1 
  names_Pop <- names_Pop[-which(names_Pop=="CMP")]
  mycolors_tmp <- colorRampPalette(c('#a50026','#ffeda0','#92c5de','#313695'))(num_Pop)
  cols <- c("gray", mycolors_tmp, "black")
  Pop <- c(Pop, rep("internal",dim(df_iqtree)[1]-num_seqs))
  names(cols) <- c("Chimpanzee", names_Pop, "internal")
}else{
  num_Pop <- length(names_Pop)
  mycolors_tmp <- colorRampPalette(c('#a50026','#ffeda0','#92c5de','#313695'))(num_Pop)
  cols <- c(mycolors_tmp, "black")
  Pop <- c(Pop, rep("internal",dim(df_iqtree)[1]-num_seqs))
  names(cols) <- c(names_Pop, "internal")
}


df_iqtree$Pop <- Pop
df_iqtree$Pop <- ifelse(df_iqtree$Pop == "CMP", "Chimpanzee", df_iqtree$Pop)

df_iqtree$branchGroup <- ifelse(grepl("CMP",df_iqtree$label), 2, 1)
df_iqtree$branchGroup <- as.factor(df_iqtree$branchGroup)

df_iqtree$SV <- df_SV_mapping[ match(df_iqtree$label, df_SV_mapping$orig_hapID), "SV"]
df_iqtree$SV <- ifelse(df_iqtree$SV == T, "SV",NA)
df_iqtree[df_iqtree$label == "CMP", "SV"] = NA

iqtree2 <- iqtree
iqtree@phylo$edge.length.mod <- iqtree@phylo$edge.length
tipidx_CMP <- which(grepl("CMP",iqtree@phylo$tip.label))
lengthidx_CMP <- match(tipidx_CMP, iqtree@phylo$edge[,2] )

if(length(lengthidx_CMP) > 1 ){
  lengthidx_root = getmode(iqtree@phylo$edge[ lengthidx_CMP,1])
  lengthidx_CMP = match(lengthidx_root, iqtree@phylo$edge[ lengthidx_CMP,1])
}

scaleFactor <- iqtree@phylo$edge.length[lengthidx_CMP] / max(iqtree@phylo$edge.length[-lengthidx_CMP])

if (scaleFactor > 1){
  newlength_CMP = iqtree@phylo$edge.length[lengthidx_CMP] / scaleFactor * 1.5
}else{
 newlength_CMP = iqtree@phylo$edge.length[lengthidx_CMP]
}
  
iqtree@phylo$edge.length[lengthidx_CMP] = newlength_CMP

func_fort <- function(tr) tr %>% fortify
fort_tree <- func_fort(iqtree@phylo)
upXlim <- max(fort_tree$x)
lowXlim <- min(fort_tree$x)

ggtree(iqtree, layout = "slanted", aes(color=Pop, linetype=branchGroup), size=1.2)  %<+% df_iqtree + 
  geom_label2(aes(label=UFboot, subset=UFboot >= 90), label.size=0.1, label.padding = unit(0, "lines"), show.legend = FALSE) +
  geom_tippoint(aes(color=SV, shape=SV), size=4) + scale_shape_manual(values=c(17,NA), breaks=c("SV")) + 
  geom_treescale(offset = -2) +  scale_x_continuous(limits = c(min(0,lowXlim), upXlim)) +  
  scale_color_manual(values=cols, breaks=names_Pop) + geom_tiplab(aes(label=label), linetype=NA, size=2.5, hjust = -1) +  guides(linetype="none", size="none", color=guide_legend(title="Population", nrow = 2), shape=guide_legend(title=paste("SV",outf_prefix,sep="\n"), nrow = 2)) + theme(legend.position="top", legend.text = element_text(size=16), legend.title = element_text(size=16), legend.key.size = unit(1, 'cm') )

print(paste(outf_dir, "/", outf_prefix,".treefile.pdf", sep=''))
ggsave(paste(outf_dir, "/", outf_prefix,".treefile.pdf", sep=''), height = 9, width = 12, limitsize=F, useDingbats=FALSE)


ggtree(iqtree, branch.length = "none", layout = "circular", aes(color=Pop, linetype=branchGroup), size=1.2)  %<+% 
  df_iqtree + geom_label2(aes(label=UFboot, subset=UFboot >= 90), label.size = 0.1, label.padding = unit(0, "lines"), show.legend = FALSE) +
  geom_tippoint(aes(color=SV, shape=SV), size=4) + scale_shape_manual(values=c(17,NA), breaks=c("SV")) + 
 # geom_treescale(offset = -2) +  scale_x_continuous(limits = c(min(0,lowXlim), upXlim)) +  
  scale_color_manual(values=cols, breaks=names_Pop) + guides(linetype="none", size="none", color=guide_legend(title="Population", nrow=2), shape=guide_legend(title=paste("SV",outf_prefix,sep="\n"), nrow = 2)) + theme(legend.position="top", legend.text = element_text(size=16), legend.title = element_text(size=16), legend.key.size = unit(1, 'cm'))

ggsave(paste(outf_dir, "/", outf_prefix,".cladogram.treefile.pdf",sep=''), height = 9, width = 12, limitsize=F, useDingbats=FALSE)


##  this does not really plot timetree; just to fulfill the requirement of the code 
ggsave(paste(outf_dir, "/", outf_prefix,".timetree.iqtree.pdf",sep=''), height = 9, width = 18, limitsize=F, useDingbats=FALSE)



q()




# plot rooted/dated iqtree
mytimetree <- read.beast(fin_timetree)

num_seqs <- length(mytimetree@phylo$tip.label)

df_mytimetree <- as.data.frame(as_tibble(mytimetree))

Pop <- sapply(strsplit(df_mytimetree[1:num_seqs,"label"], split = "_"), function(x) x[[1]])
Pop <- c(Pop, rep("internal",dim(df_mytimetree)[1]-num_seqs))

df_mytimetree$Pop <- Pop
df_mytimetree$Pop <- ifelse(df_mytimetree$Pop == "CMP", "Chimpanzee", df_mytimetree$Pop)

df_mytimetree$branchGroup <- ifelse(df_mytimetree$Pop == "CMP", 2, 1)
df_mytimetree$branchGroup <- as.factor(df_mytimetree$branchGroup)
df_mytimetree$SV <- df_SV_mapping[ match(df_mytimetree$label, df_SV_mapping$orig_hapID), "SV"]
df_mytimetree[df_mytimetree$Pop == "CMP", "SV"] = FALSE
df_mytimetree$INV <- as.factor(df_mytimetree$SV)

p1 <- ggtree(mytimetree, aes(color=Pop, linetype=branchGroup), size=1.2)  %<+% df_mytimetree + 
  geom_tippoint(aes(color=Pop, shape=SV), size=4.5) + scale_shape_manual(values=c(NA,17),breaks = levels(df_mytimetree$SV)) + geom_treescale(x=0.2) +
  scale_color_manual(values=cols, breaks=c("AFR","AMR","EAS","EUR","SAS")) + 
  geom_range("CI_date", color="gray",size=2,alpha=0.7) + 
#  geom_tiplab(aes(label=mapping_hap_INV[ match(df_mytimetree$label, mapping_hap_INV$orig_hapID),"hapID"], subset=INV), linetype=NA, size=5, hjust=-0.1) +  
  guides(linetype="none", size="none", color=guide_legend(title="Population", nrow=1), shape=guide_legend(title=paste("SV",outf_prefix,sep="\n"), nrow=1)) + theme_tree2() + 
  theme(legend.position="top", legend.justification='left', legend.text = element_text(size=14), legend.title = element_text(size=12), legend.key.size = unit(1, 'cm'), axis.text.x = element_text(size=14)) + xlab("Million years") 

p2 <- revts(p1) 
p3 <- p2 + scale_x_continuous(breaks=c(-6:0), labels=abs(-6:0))


outGroupID <- which(grepl("CMP", mytimetree@phylo$tip.label))
for (g in outGroupID){
  mytimetree2 <- treeio::drop.tip(mytimetree, g)
}
num_seqs <- length(mytimetree2@phylo$tip.label)

df_mytimetree2 <- as.data.frame(as_tibble(mytimetree2))

Pop <- sapply(strsplit(df_mytimetree2[1:num_seqs,"label"], split = "_"), function(x) x[[1]])
Pop <- c(Pop, rep("internal",dim(df_mytimetree2)[1]-num_seqs))

df_mytimetree2$Pop <- Pop
df_mytimetree2$Pop <- ifelse(grepl("CMP",df_mytimetree2$Pop), "Chimpanzee", df_mytimetree2$Pop)

df_mytimetree2$branchGroup <- ifelse(grepl("CMP",df_mytimetree2$Pop), 2, 1)
df_mytimetree2$branchGroup <- as.factor(df_mytimetree2$branchGroup)
df_mytimetree2$SV <- df_SV_mapping[ match(df_mytimetree2$label, df_SV_mapping$orig_hapID), "SV"]
df_mytimetree2[ grepl("CMP",df_mytimetree2$Pop), "SV"] = FALSE
df_mytimetree2$SV <- as.factor(df_mytimetree2$SV) 


p4 <- ggtree(mytimetree2, aes(color=Pop, linetype=branchGroup), size=1.2)  %<+% df_mytimetree2 + 
  geom_tippoint(aes(color=Pop, shape=SV), size=4.5) + scale_shape_manual(values=c(NA,17),breaks = levels(df_mytimetree2$SV)) + geom_treescale(x=0.1) + 
  scale_color_manual(values=cols, breaks=c("AFR","AMR","EAS","EUR","SAS")) + 
  geom_range("CI_date", color="gray",size=2,alpha=0.7) + 
  geom_tiplab(aes(label= df_SV_mapping[ match(df_mytimetree2$label, df_SV_mapping$orig_hapID),"hapID"], subset=SV==T), linetype=NA, size=5, hjust=-0.1) +  
  guides(linetype="none", size="none", color=guide_legend(title="Population", nrow=1), shape=guide_legend(title=paste("SV",outf_prefix,sep="\n"), nrow=1)) + theme_tree2() + 
  theme(legend.position="top", legend.justification='left', legend.text = element_text(size=14), legend.title = element_text(size=12), legend.key.size = unit(1, 'cm'), axis.text.x = element_text(size=14)) + xlab("Million years")

upperTMRCA <- round(abs(min(unlist(mytimetree2@data$CI_date),na.rm = T)), digits = 1)
p5 <- revts(p4) 
breaks=seq(-upperTMRCA,0, upperTMRCA/2)
labels = abs(seq(-upperTMRCA,0,upperTMRCA/2))
for (i in 1:length(labels)) { if (i%%2==0) labels[i]=""}
p6 <- p5 + scale_x_continuous(breaks=breaks, labels=labels)

ggarrange(p3,p6,ncol=2,common.legend = T, widths = c(6,12))

ggsave(paste(outf_prefix,".timetree.iqtree.pdf",sep=''), height = 9, width = 18, limitsize=F, useDingbats=FALSE)
print (unlist(mytimetree2@data$CI_date))



q()
# date tree using APE's penalized algorithm
rt_iqtree2 <- root(iqtree2@phylo, outgroup = "CMP")
numZeroBL = table(rt_iqtree2$edge.length == 0)[2]
if (numZeroBL > 0.3 * Ntip(rt_iqtree2)) {
#  rt_iqtree2$edge.length = ifelse(rt_iqtree2$edge.length == 0, min(min(rt_iqtree2$edge.length[rt_iqtree2$edge.length!=0])/1e4, 1e-10), rt_iqtree2$edge.length)
  rt_iqtree2$edge.length = rt_iqtree2$edge.length + 1e-10
}
my_node <- getMRCA(rt_iqtree2, tip=c("CMP", rt_iqtree2$tip.label[2])) 
mycalibration <- makeChronosCalib(phy = rt_iqtree2, node = my_node, age.min=4, age.max = 6 )
count_eval = 0
stopTrying = FALSE
list_mytimetree <- vector(mode="list", length=10)
for (i in 1:10){
  convergence = FALSE
  while (convergence != TRUE){
    count_eval = count_eval + 1
    print(count_eval)
    tryCatch(
      mytimetree <- chronos(rt_iqtree2, lambda = 1, calibration = mycalibration, model="discrete", control = chronos.control(epsilon = 1e-8, eval.max = 1e5, dual.iter.max = 1000, tol=1e-8, nb.rate.cat=1))
      ,error=function(e) e)
	if(inherits(mytimetree,"error")) next
    convergence = attr(mytimetree, "convergence")
    print(convergence)
    if (convergence == TRUE){
      list_mytimetree[[i]] <- mytimetree
      break
    }
    if (count_eval > 100){
      if (i!= 1){
        stopTrying = TRUE
        break
      }
    }
  }
  if (stopTrying == TRUE) break
}

list_mytimetree <- list_mytimetree[lapply(list_mytimetree,length)>0]
list_logLik <- unlist(lapply(list_mytimetree, function(x) attr(x,"PHII")$logLik))
mytimetree <- list_mytimetree[[which.max(list_logLik)]]

## plot dated (APE) tree 
write.tree(mytimetree, paste(outf_prefix,".mytimetree.nwk",sep=''))
mytimetree <- read.iqtree(paste(outf_prefix,".mytimetree.nwk",sep=''))


num_seqs <- length(mytimetree@phylo$tip.label)

df_mytimetree <- as.data.frame(as_tibble(mytimetree))

Pop <- sapply(strsplit(df_mytimetree[1:num_seqs,"label"], split = "_"), function(x) x[[1]])
Pop <- c(Pop, rep("internal",dim(df_mytimetree)[1]-num_seqs))

df_mytimetree$Pop <- Pop
df_mytimetree$Pop <- ifelse(df_mytimetree$Pop == "CMP", "Chimpanzee", df_mytimetree$Pop)

df_mytimetree$branchGroup <- ifelse(df_mytimetree$label == "CMP", 2, 1)
df_mytimetree$branchGroup <- as.factor(df_mytimetree$branchGroup)
df_mytimetree$INV <- mapping_hap_INV[ match(df_mytimetree$label,mapping_hap_INV$orig_hapID), "INV"]
df_mytimetree[df_mytimetree$label == "CMP", "INV"] = FALSE


p1 <- ggtree(mytimetree, aes(color=Pop, linetype=branchGroup), size=1.2)  %<+% df_mytimetree + 
  geom_label2(aes(label=UFboot, subset=UFboot > 70), show.legend = FALSE) +
  geom_tippoint(aes(color=Pop, shape=INV), size=4.5) + scale_shape_manual(values=c(20,17)) + geom_treescale(x=2) +
  scale_color_manual(values=cols, breaks=c("AFR","AMR","EAS","EUR","SAS")) + 
  geom_tiplab(aes(label=mapping_hap_INV[ match(df_mytimetree$label, mapping_hap_INV$orig_hapID),"hapID"], subset=INV), linetype=NA, size=5, hjust=-0.1) +  
  guides(linetype="none", size="none", color=guide_legend(title="Population", nrow=1), shape=guide_legend(title=paste("Inversion",outf_prefix,sep="\n"), nrow=1)) + theme_tree2() + 
  theme(legend.position="top", legend.justification='left', legend.text = element_text(size=14), legend.title = element_text(size=12), legend.key.size = unit(1, 'cm'), axis.text.x = element_text(size=14)) + xlab("Million years") 

p2 <- revts(p1) 
p3 <- p2 + scale_x_continuous(breaks=c(-6:0), labels=abs(-6:0))


#mytimetree2 <- drop.tip(mytimetree@phylo, "CMP")
keeptips <- mytimetree@phylo$tip.label[ mytimetree@phylo$tip.label != "CMP"]
mytimetree2 <- tree_subset(mytimetree, keeptips, levels_back = 100, group_node = F, root_edge = T)

num_seqs <- length(mytimetree2@phylo$tip.label)

df_mytimetree2 <- as.data.frame(as_tibble(mytimetree2))

Pop <- sapply(strsplit(df_mytimetree2[1:num_seqs,"label"], split = "_"), function(x) x[[1]])
Pop <- c(Pop, rep("internal",dim(df_mytimetree2)[1]-num_seqs))

df_mytimetree2$Pop <- Pop
df_mytimetree2$Pop <- ifelse(df_mytimetree2$Pop == "CMP", "Chimpanzee", df_mytimetree2$Pop)

df_mytimetree2$branchGroup <- ifelse(df_mytimetree2$Pop == "CMP", 2, 1)
df_mytimetree2$branchGroup <- as.factor(df_mytimetree2$branchGroup)
df_mytimetree2$INV <- mapping_hap_INV[ match(df_mytimetree2$label,mapping_hap_INV$orig_hapID), "INV"]
df_mytimetree2[df_mytimetree2$Pop == "CMP", "INV"] = FALSE




p4 <- ggtree(mytimetree2, aes(color=Pop, linetype=branchGroup), size=1.2)  %<+% df_mytimetree2 + 
  geom_tippoint(aes(color=Pop, shape=INV), size=4.5) + scale_shape_manual(values=c(20,17)) + geom_treescale(x=.6) + 
  scale_color_manual(values=cols, breaks=c("AFR","AMR","EAS","EUR","SAS")) + 
  geom_tiplab(aes(label=mapping_hap_INV[ match(df_mytimetree2$label, mapping_hap_INV$orig_hapID),"hapID"], subset=INV), linetype=NA, size=5, hjust=-0.1) +  
  guides(linetype="none", size="none", color=guide_legend(title="Population", nrow=1), shape=guide_legend(title=paste("Inversion",outf_prefix,sep="\n"), nrow=1)) + theme_tree2() + 
  theme(legend.position="top", legend.justification='left', legend.text = element_text(size=14), legend.title = element_text(size=12), legend.key.size = unit(1, 'cm'), axis.text.x = element_text(size=14)) + xlab("Million years")


p5 <- revts(p4) 
breaks=seq(-6,0,0.25)
labels = abs(seq(-6,0,0.25))
for (i in 1:length(labels)) { if (i%%2==0) labels[i]=""}
p6 <- p5 + scale_x_continuous(breaks=breaks, labels=labels)

ggarrange(p3,p6,ncol=2,common.legend = T, widths = c(12,6))

ggsave(paste(outf_prefix,".timetree.ape.pdf",sep=''), height = 9, width = 18, limitsize=F, useDingbats=FALSE)

