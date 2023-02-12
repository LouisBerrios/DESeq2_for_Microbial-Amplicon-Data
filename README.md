# DESeq2_for_Microbial-Amplicon-Data

This pipeline will allow you to analyze your amplicon data (in this case bacterial communities) using the `DESeq2 R package`. It assumes that OTU/ASV tables have been generated and merged with sample/meta data. Likewise, the pipeline assumes that a `phyloseq object` has been generated – which we will need to convert into a `DESeq object`.

 **For more information about the `DESeq2 R package`, [click here to read more about it](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).** 

---

For additional information on how to process sample reads and derived sequence count tables, see my other respository titled `BishopPine_Tripartite_Interactions`.

---

# **A Little Blurb about the *DESeq2* package**

Originally designed to analyze RNA-seq data, the `DESeq2 package` measures differential expression/abundances between groups, categories, and/or conditions. Since ASV/OTU tables follow a similar count data structure as count tables for RNA-seq outputs, the `DESeq2` package can also be used to detect abundance differences between groups. **Therefore, we can use `DESeq2` to determine which microbial taxa are enriched or depleted between our experimental groups (e.g., between samples at time '0' compared to those at time point '4').** With this information, we gain a better understanding of which taxa 'indicate' or 'signal' something potentially important about our focal system. 

# **Running the *DESeq2* Pipeline**
## Package Installation
To install *DESeq2*, enter the following:
```{}
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```
## Bacterial Samples
```{r}
library(DESeq2)

# Since DESeq2 only compares differences between 2 groups, we will first subset our data into 2 age groups.
# Unique age groups = 0, 1, 2, 4, 8, 16, 24, 36
# Let's first look at differences between samples at age 0 and age 36

# Use the rarefied and non-transforned dataset
bac.crap.even2 <- readRDS("BAC.Phylo.Rarefied.rds")
# Extract the sample data from the phyloseq object
bac.df.num2 <- as(sample_data(bac.crap.even2), "data.frame")
# Convert the numerical variable to a categorical variable
bac.df.num2$SampleAge <- as.factor(bac.df.num2$SampleAge)
# Re-assign the modified dataframe to the phyloseq object
bac.ps.num2 = sample_data(bac.crap.even2) <- bac.df.num2

# Save the modified phyloseq 
saveRDS(bac.crap.even2, "BAC.Phylo.RARE.NUMtoCAT.rds")

# Let's first look at differences between samples at age 0 and age 4

bac.crap.0.4 <- subset_samples(bac.crap.even2, SampleAge == "0" | SampleAge == "4")

# Convert phyloseq to DESeq2 format
bac.ddsFUN.0.4 <- phyloseq_to_deseq2(bac.crap.0.4, ~ SampleAge)
# Calculate size factors using edgeR
library(edgeR)
sizeFactors(bac.ddsFUN.0.4) <- calcNormFactors(counts(bac.ddsFUN.0.4))
# Run DESeq function
bac.ddsFUN2 = DESeq(bac.ddsFUN.0.4, test="Wald", fitType="parametric")
# Generate sigtab output table
bac.resFUN2 <- results(bac.ddsFUN2, cooksCutoff = FALSE)
alpha = 0.1
bac.sigtabFUN2 <- bac.resFUN2[which(bac.resFUN2$padj < alpha), ]
bac.sigtabFUN2 <- cbind(as(bac.sigtabFUN2, "data.frame"), as(tax_table(bac.crap.even2)[rownames(bac.sigtabFUN2), ], "matrix"))
x = tapply(bac.sigtabFUN2$log2FoldChange, bac.sigtabFUN2$Phylum, function(x) max(x))
x = sort(x, TRUE)
bac.sigtabFUN2$Phylum = factor(as.character(bac.sigtabFUN2$Phylum), levels=names(x))
x = tapply(bac.sigtabFUN2$log2FoldChange, bac.sigtabFUN2$Genus, function(x) max(x))
bac.sigtabFUN2$Genus = factor(as.character(bac.sigtabFUN2$Genus), levels=names(x))
# remove NA's
bac.sigtabFUN2.1 <- subset(bac.sigtabFUN2, Genus != "NA")
# Remove samples that have no 'order' classification
bac.sigtabFUN2.2 <- subset(bac.sigtabFUN2.1, Order !="")
# Add Title as a column
bac.sigtabFUN2.2$Title <- "Differentially Abundant Bacteria [Age 0 vs. Age 4]"
# Plot using ggplot2
library(ggplot2)
ggplot(bac.sigtabFUN2.2, aes(x=Order, y=log2FoldChange)) + geom_point(size=4) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + theme_bw() + coord_flip() + xlab("") + ylab("[Log2] Fold-Change") + theme(axis.text.y = element_text(size = 8, face="bold.italic")) + 
  theme(plot.title = element_text(size = 10)) + geom_hline(yintercept=0, linetype ="dashed") + facet_wrap(~Title) + theme(strip.text = element_text(face="bold"))
# May have to do bacteria at the 'order' level at first because there are so many differentially abundant genera
# Save plot
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
