library(DESeq2)
library(phyloseq)

####OTU Analysis####

#Import data - OTU file comes from transpose_mothur_shared_file.R
#Taxonomy comes from AddSpeciesToTaxonomy.R
#Sample data is summary of metadata with mothur group names as row names
OTU.matrix <- data.matrix(read.table(file="FF2_mothur_Transposed.shared", header=TRUE, row.names=1, sep="\t"))
tax.matrix <- as.matrix(read.table(file="FF2_RDP_Blast.taxonomy", header=TRUE, row.names=1, sep="\t"))
sample <- read.csv(file="FF2.meta.csv", header=TRUE, row.names=1, sep=",")

##Create phyloseq object##
OTU.p <- otu_table(OTU.matrix, taxa_are_rows =TRUE)
TAX.p <- tax_table(tax.matrix)
sample.p <- sample_data(sample)
physeq <- phyloseq(OTU.p,TAX.p,sample.p)

###Filter sequences###

##Create filter to have minimum 51 counts (0.001%) and 5% prevalence
physeq_51 <- filter_taxa(physeq, function(x) sum(x) >=51, prune=TRUE)

# Compute prevalence of each feature, store as data.frame
prevdf51 <- apply(X = otu_table(physeq_51),
                  MARGIN = ifelse(taxa_are_rows(physeq_51), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf51 <- data.frame(Prevalence = prevdf51,
                       TotalAbundance = taxa_sums(physeq_51),
                       tax_table(physeq_51))

prevalenceThreshold <- 0.05 * nsamples(physeq_51)
keepTaxa2 <- rownames(prevdf51)[(prevdf51$Prevalence >= prevalenceThreshold)]
physeq_51_5 <- prune_taxa(keepTaxa2, physeq_51)

###DESeq2 Analysis###

#convert to factors
sample$Treatment <- relevel(as.factor(sample$Treatment), ref = "Water")

#Run analysis looking at differences in Diet
OTU51_5.matrix <- as.data.frame(otu_table(physeq_51_5))
DE_Diet_data51_5 <- DESeqDataSetFromMatrix(countData = OTU51_5.matrix, colData = sample, design = ~Treatment)
DE_Diet51_5 <- DESeq(DE_Diet_data51_5)

alpha_cut <- 0.05
#Compare each treatment to water

treats <- factor(c("Amylopectin", "Banana", "CornStarch", "HAM2", "HAM4", "Potato", "PS_Ba", "PS_Rb", "RS3_Extracted",
                   "RS3_Whole", "Tapioca", "TigerNut"))

for (fac in  treats){
  res51_5 <- results(DE_Diet51_5, contrast = c("Treatment", fac, "Water"), cooksCutoff = FALSE)
  sigtab51_5 <- res51_5[which(res51_5$padj < alpha_cut), ]
  sigtab51_5 = cbind(as(sigtab51_5, "data.frame"), as(tax_table(physeq_51_5)[rownames(sigtab51_5), ], "matrix"))
  fn <- paste("DESEQ_", fac, ".csv", sep = "")
  write.csv(sigtab51_5,file = fn) 
}

