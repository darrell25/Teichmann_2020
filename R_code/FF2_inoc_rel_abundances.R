library(phyloseq)
library(data.table)
library(dplyr)
library(ggplot2)

OTU.df <- data.frame(read.table(file="FF2_mothur_Transposed.shared", header=TRUE, row.names=1, sep="\t"))
tax.df <- data.frame(read.table(file="FF2_RDP_Blast.taxonomy", header=TRUE, row.names=1, sep="\t"))
sample <- read.csv(file="FF2.meta.csv", header=TRUE, row.names=1, sep=",")

tax.df$Species <- make.names (ifelse(tax.df$Genus==tax.df$Species,tax.df$Genus,paste(tax.df$Genus, tax.df$Species, sep="_")),  
                              unique=FALSE)
tax.df$Species <- make.names (ifelse(grepl("_",tax.df$Species,fixed=TRUE ), tax.df$Species, paste(tax.df$Species, "unclassified", sep = "_")), 
                              unique = FALSE)

##Create phyloseq object##
OTU.p <- otu_table(OTU.df, taxa_are_rows =TRUE)
TAX.p <- tax_table(as.matrix(tax.df))
sample.p <- sample_data(sample)
physeq <- phyloseq(OTU.p,TAX.p,sample.p)

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

#Get the average relative abundance for each of the Inocula, not including the
#conditions where organisms were supplemented in
physeq_51_5 <- transform_sample_counts(physeq_51_5, function(x) x/sum(x))
physeq_51_5 <- subset_samples(physeq_51_5, !(Treatment %in% c("PS_Rb", "PS_Ba")))

physeq_T1 <- data.frame (otu_table(subset_samples(physeq_51_5, Inoculum == "T01")))
physeq_T2 <- data.frame (otu_table(subset_samples(physeq_51_5, Inoculum == "T02")))
physeq_T3 <- data.frame (otu_table(subset_samples(physeq_51_5, Inoculum == "T03")))
physeq_T4 <- data.frame (otu_table(subset_samples(physeq_51_5, Inoculum == "T04")))
physeq_T5 <- data.frame (otu_table(subset_samples(physeq_51_5, Inoculum == "T05")))
physeq_T6 <- data.frame (otu_table(subset_samples(physeq_51_5, Inoculum == "T06")))
physeq_T7 <- data.frame (otu_table(subset_samples(physeq_51_5, Inoculum == "T07")))
physeq_T8 <- data.frame (otu_table(subset_samples(physeq_51_5, Inoculum == "T08")))
physeq_T9 <- data.frame (otu_table(subset_samples(physeq_51_5, Inoculum == "T09")))
physeq_T10 <- data.frame (otu_table(subset_samples(physeq_51_5, Inoculum == "T10")))
physeq_T11 <- data.frame (otu_table(subset_samples(physeq_51_5, Inoculum == "T11")))

Inoc.pct <- as.data.table(data.frame(OTU=rownames(tax_table(physeq_51_5)), Species=tax_table(physeq_51_5)[,"Species"], 
                          T1 = round(rowMeans(physeq_T1)*100,3), 
                          T2 = round(rowMeans(physeq_T2)*100,3), T3 = round(rowMeans(physeq_T3)*100,3), 
                          T4 = round(rowMeans(physeq_T4)*100,3), T5 = round(rowMeans(physeq_T5)*100,3),
                          T6 = round(rowMeans(physeq_T6)*100,3), T7 = round(rowMeans(physeq_T7)*100,3),
                          T8 = round(rowMeans(physeq_T8)*100,3), T9 = round(rowMeans(physeq_T9)*100,3),
                          T10 = round(rowMeans(physeq_T10)*100,3), T11 = round(rowMeans(physeq_T11)*100,3)))

fwrite(Inoc.pct, file="Species_inoculum_percent_abundances.txt", sep = "\t")
