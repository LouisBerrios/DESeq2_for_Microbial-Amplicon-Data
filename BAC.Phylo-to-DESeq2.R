#r scripts for bacteria DESeq2 analysis of the CRAP samples

library(phyloseq)
library(vegan)
library(ggplot2)
library(gam)
library(BiodiversityR)
library(microbiome)

setwd("/Users/louisberrios/Documents/Coprophile_Data/Group_Copro_CRAP")

######################################################################
########################### import data ##############################
#######################################################################

bac.otus.raw <- read.csv("otus.bac.itags.csv",header=T,row.names=1)
bac.taxonomy<-bac.otus.raw$taxonomy; write.csv(bac.taxonomy,"bac_taxonomy.csv")
bac.otus<-bac.otus.raw[,1:64]
bac.meta<-read.csv("iTags_meta.csv", stringsAsFactors=FALSE)
bac.meta2<-bac.meta[match(colnames(bac.otus),bac.meta$JGI_Sample),]
colnames(bac.otus)<-bac.meta2$Sample.ID
rownames(bac.meta2)<-bac.meta2$Sample.ID

sum(bac.otus)
bac.otus <- as(as.matrix(bac.otus), "matrix")
bac.sampleData <- sample_data(bac.meta2)

##############################################################################################
#create a filter to remove potential false positive occurrences based on relative abundance ##
##############################################################################################

bac.rel.abund<-decostand(bac.otus,"total",MARGIN=2) #standardize by total sequence no. per sample
bac.max.rel.abund<-apply(bac.rel.abund,1,"max") #calculate maximum relative abund. per species
bac.abund.cutoff.1000<-bac.max.rel.abund/1000 #create a cutoff after which entries will be ignored
bac.abund.cutoff.10000<-bac.max.rel.abund/10000 #create a cutoff after which entries will be ignored
bac.filter.1000<-decostand(bac.otus,'pa') #crete a pa dataset which will have 1's & 0's
bac.filter.10000<-decostand(bac.otus,'pa') #and can be multiplied to act as a mask

#use a for loop to apply the species specific cutoffs to the filter tables

for (i in 1:length(row.names(bac.otus))){
  
  bac.filter.1000[i,(bac.rel.abund[i,]<bac.abund.cutoff.1000[i])]<-0
  bac.filter.10000[i,(bac.rel.abund[i,]<bac.abund.cutoff.10000[i])]<-0
  
}

bac.otus.1000<-bac.otus*bac.filter.1000 #multiply by the filter
bac.otus.10000<-bac.otus*bac.filter.10000 #multiply by the filter

############## create phyloseq objects ###############

bac.OTU = otu_table(bac.otus.10000, taxa_are_rows = TRUE)

bac.taxmat <- read.csv("bac_taxonomy_table.csv",header=T,row.names=1)
bac.taxmat <- as(as.matrix(bac.taxmat),"matrix")
bac.TAX = tax_table(bac.taxmat)

#TAX <- head(TAX,n=120L)
bac.phylo = phyloseq(data=bac.OTU, sample_metadata=bac.sampleData,observation_metadata=bac.TAX)

bac.phylo = subset_samples(bac.phylo,Sample.Type == "Single") #drops the aggregated decomp samples
#save phyloseq object
saveRDS(bac.phylo, "BAC.Phylo.rds")
######################################################
## rarefy to even depth #############################
#####################################################

bac.crap.even<-rarefy_even_depth(bac.phylo,sample.size=min(sample_sums(bac.phylo)),rngseed=22)
#save phyloseq object
saveRDS(bac.crap.even, "BAC.Phylo.Rarefied.rds")

###########################################################
######## Run DESeq2 Analysis ##############################
###########################################################

library(DESeq2)

# Since DESeq2 only compares differences between 2 groups, we will first subset out data into 2 age groups
## Unique age groups = 0, 1, 2, 4, 8, 16, 24, 36
## Let's first look at differences between samples at age 0 and age 36

bac.crap.0.36 <- subset_samples(bac.crap.even, SampleAge == "0" | SampleAge == "36")

# Convert phyloseq to DESeq2 format
bac.ddsFUN.0.36 <- phyloseq_to_deseq2(bac.crap.0.36, ~ SampleAge)
# calculate size factors using edgeR
library(edgeR)
sizeFactors(bac.ddsFUN.0.36) <- calcNormFactors(counts(bac.ddsFUN.0.36))
# run DESeq function
bac.ddsFUN1 = DESeq(bac.ddsFUN.0.36, test="Wald", fitType="parametric")
# Generate sigtab output table
bac.resFUN <- results(bac.ddsFUN1, cooksCutoff = FALSE)
alpha = 0.1
bac.sigtabFUN1 <- bac.resFUN[which(bac.resFUN$padj < alpha), ]
bac.sigtabFUN1 <- cbind(as(bac.sigtabFUN1, "data.frame"), as(tax_table(bac.crap.even)[rownames(bac.sigtabFUN1), ], "matrix"))
x = tapply(bac.sigtabFUN1$log2FoldChange, bac.sigtabFUN1$Phylum, function(x) max(x))
x = sort(x, TRUE)
bac.sigtabFUN1$Phylum = factor(as.character(bac.sigtabFUN1$Phylum), levels=names(x))
x = tapply(bac.sigtabFUN1$log2FoldChange, bac.sigtabFUN1$Genus, function(x) max(x))
bac.sigtabFUN1$Genus = factor(as.character(bac.sigtabFUN1$Genus), levels=names(x))
# remove NA's
bac.sigtabFUN1.1 <- subset(bac.sigtabFUN1, Species != "NA")
# Add Title as a column
bac.sigtabFUN1.1$Title <- "Differentially Abundant Bacteria [Age 0 vs. Age 36]"
# plot using ggplot2
library(ggplot2)
ggplot(bac.sigtabFUN1.1, aes(x=Order, y=log2FoldChange)) + geom_point(size=4) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) 
+ theme_bw() + coord_flip() + xlab("") + ylab("[Log2] Fold-Change") + theme(axis.text.y = element_text(size = 8, face="bold.italic")) + 
  theme(plot.title = element_text(size = 10)) + geom_hline(yintercept=0, linetype ="dashed") + facet_wrap(~Title) + theme(strip.text = element_text(face="bold"))
# may have to do bacteria at the 'order' level at first because there are so many differentially abundant genera
# save plot
ggsave(
  "CRAP.BAC.Age.0-36.ORDER.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5.5,
  height = 4.5,
  units = c("in"),
  dpi = 300)

## Let's first look at differences between samples at age 0 and age 4

bac.crap.0.4 <- subset_samples(bac.crap.even, SampleAge == "0" | SampleAge == "4")

# Convert phyloseq to DESeq2 format
bac.ddsFUN.0.4 <- phyloseq_to_deseq2(bac.crap.0.4, ~ SampleAge)
# calculate size factors using edgeR
library(edgeR)
sizeFactors(bac.ddsFUN.0.4) <- calcNormFactors(counts(bac.ddsFUN.0.4))
# run DESeq function
bac.ddsFUN2 = DESeq(bac.ddsFUN.0.4, test="Wald", fitType="parametric")
# Generate sigtab output table
bac.resFUN2 <- results(bac.ddsFUN2, cooksCutoff = FALSE)
alpha = 0.1
bac.sigtabFUN2 <- bac.resFUN2[which(bac.resFUN2$padj < alpha), ]
bac.sigtabFUN2 <- cbind(as(bac.sigtabFUN2, "data.frame"), as(tax_table(bac.crap.even)[rownames(bac.sigtabFUN2), ], "matrix"))
x = tapply(bac.sigtabFUN2$log2FoldChange, bac.sigtabFUN2$Phylum, function(x) max(x))
x = sort(x, TRUE)
bac.sigtabFUN2$Phylum = factor(as.character(bac.sigtabFUN2$Phylum), levels=names(x))
x = tapply(bac.sigtabFUN2$log2FoldChange, bac.sigtabFUN2$Genus, function(x) max(x))
bac.sigtabFUN2$Genus = factor(as.character(bac.sigtabFUN2$Genus), levels=names(x))
# remove NA's
bac.sigtabFUN2.1 <- subset(bac.sigtabFUN2, Genus != "NA")
#remove samples that have no 'order' classification
bac.sigtabFUN2.2 <- subset(bac.sigtabFUN2.1, Genus !="")
# Add Title as a column
bac.sigtabFUN2.1$Title <- "Differentially Abundant Bacteria [Age 0 vs. Age 4]"
# plot using ggplot2
library(ggplot2)
ggplot(bac.sigtabFUN2.1, aes(x=Order, y=log2FoldChange)) + geom_point(size=4) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) 
+ theme_bw() + coord_flip() + xlab("") + ylab("[Log2] Fold-Change") + theme(axis.text.y = element_text(size = 8, face="bold.italic")) + 
  theme(plot.title = element_text(size = 10)) + geom_hline(yintercept=c(-1.5, 1.5), linetype ="dashed") + facet_wrap(~Title) + theme(strip.text = element_text(face="bold"))
# may have to do bacteria at the 'order' level at first because there are so many differentially abundant genera
# save plot
ggsave(
  "CRAP.BAC.Age.0-4.ORDER.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5.5,
  height = 4.5,
  units = c("in"),
  dpi = 300)

################ Now let's compare age 0 to age 8 ####################

bac.crap.0.8 <- subset_samples(bac.crap.even, SampleAge == "0" | SampleAge == "8")

# Convert phyloseq to DESeq2 format
bac.ddsFUN.0.8 <- phyloseq_to_deseq2(bac.crap.0.8, ~ SampleAge)
# calculate size factors using edgeR
library(edgeR)
sizeFactors(bac.ddsFUN.0.8) <- calcNormFactors(counts(bac.ddsFUN.0.8))
# run DESeq function
bac.ddsFUN3 = DESeq(bac.ddsFUN.0.8, test="Wald", fitType="parametric")
# Generate sigtab output table
bac.resFUN3 <- results(bac.ddsFUN3, cooksCutoff = FALSE)
alpha = 0.1
bac.sigtabFUN3 <- bac.resFUN3[which(bac.resFUN3$padj < alpha), ]
bac.sigtabFUN3 <- cbind(as(bac.sigtabFUN3, "data.frame"), as(tax_table(bac.crap.even)[rownames(bac.sigtabFUN3), ], "matrix"))
x = tapply(bac.sigtabFUN3$log2FoldChange, bac.sigtabFUN3$Phylum, function(x) max(x))
x = sort(x, TRUE)
bac.sigtabFUN3$Phylum = factor(as.character(bac.sigtabFUN3$Phylum), levels=names(x))
x = tapply(bac.sigtabFUN3$log2FoldChange, bac.sigtabFUN3$Genus, function(x) max(x))
bac.sigtabFUN3$Genus = factor(as.character(bac.sigtabFUN3$Genus), levels=names(x))
# remove NA's
bac.sigtabFUN3.1 <- subset(bac.sigtabFUN3, Genus != "NA")
# remove samples with no 'order' classification
bac.sigtabFUN3.2 <- subset(bac.sigtabFUN3.1, Order != "")
# Add Title as a column
bac.sigtabFUN3.2$Title <- "Differentially Abundant Bacteria [Age 0 vs. Age 8]"
# plot using ggplot2
library(ggplot2)
ggplot(bac.sigtabFUN3.1, aes(x=Order, y=log2FoldChange)) + geom_point(size=4) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) 
+ theme_bw() + coord_flip() + xlab("") + ylab("[Log2] Fold-Change") + theme(axis.text.y = element_text(size = 8, face="bold.italic")) + 
  theme(plot.title = element_text(size = 10)) + geom_hline(yintercept=c(-1.5, 1.5), linetype ="dashed") + facet_wrap(~Title) + theme(strip.text = element_text(face="bold"))
# may have to do bacteria at the 'order' level at first because there are so many differentially abundant genera
# save plot
ggsave(
  "CRAP.BAC.Age.0-8.ORDER.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5.5,
  height = 4.5,
  units = c("in"),
  dpi = 300)

