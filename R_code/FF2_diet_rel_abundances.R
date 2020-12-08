library(phyloseq)
library(data.table)
library(dplyr)
library(ggplot2)

OTU.df <- data.frame(read.table(file="FF2_mothur_Transposed.shared", header=TRUE, row.names=1, sep="\t"))
tax.df <- data.frame(read.table(file="FF2_RDP_Blast.taxonomy", header=TRUE, row.names=1, sep="\t"))
sample <- read.csv(file="FF2.meta.csv", header=TRUE, row.names=1, sep=",")

#Generate genus_species names
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

#Get the average relative abundance for each of the diet conditions for each member 
#of the taxonomic level, in this case phylum - can be changed to other
physeq_51_5 <- tax_glom(physeq_51_5, taxrank = "Phylum")
physeq_51_5 <- transform_sample_counts(physeq_51_5, function(x) x/sum(x))

physeq_wat <- data.frame (otu_table(subset_samples(physeq_51_5, Treatment == "Water")))
physeq_ap <- data.frame (otu_table(subset_samples(physeq_51_5, Treatment == "Amylopectin")))
physeq_bn <- data.frame (otu_table(subset_samples(physeq_51_5, Treatment == "Banana")))
physeq_cs <- data.frame (otu_table(subset_samples(physeq_51_5, Treatment == "CornStarch")))
physeq_ham2 <- data.frame (otu_table(subset_samples(physeq_51_5, Treatment == "HAM2")))
physeq_ham4 <- data.frame (otu_table(subset_samples(physeq_51_5, Treatment == "HAM4")))
physeq_pot <- data.frame (otu_table(subset_samples(physeq_51_5, Treatment == "Potato")))
physeq_ba <- data.frame (otu_table(subset_samples(physeq_51_5, Treatment == "PS_Ba")))
physeq_rb <- data.frame (otu_table(subset_samples(physeq_51_5, Treatment == "PS_Rb")))
physeq_r3e <- data.frame (otu_table(subset_samples(physeq_51_5, Treatment == "RS3_Extracted")))
physeq_r3w <- data.frame (otu_table(subset_samples(physeq_51_5, Treatment == "RS3_Whole")))
physeq_tap <- data.frame (otu_table(subset_samples(physeq_51_5, Treatment == "Tapioca")))
physeq_tn <- data.frame (otu_table(subset_samples(physeq_51_5, Treatment == "TigerNut")))

#update Phylum here if using other taxonomic level
Diet.pct <- as.data.table(data.frame(OTU=rownames(tax_table(physeq_51_5)), Phylum=tax_table(physeq_51_5)[,"Phylum"], 
                          Mean_Wat = round(rowMeans(physeq_wat)*100,3), 
                          Mean_Ap = round(rowMeans(physeq_ap)*100,3), Mean_Bn = round(rowMeans(physeq_bn)*100,3), 
                          Mean_CS = round(rowMeans(physeq_cs)*100,3), Mean_HAM2 = round(rowMeans(physeq_ham2)*100,3),
                          Mean_HAM4 = round(rowMeans(physeq_ham4)*100,3), Mean_Pot = round(rowMeans(physeq_pot)*100,3),
                          Mean_Ba = round(rowMeans(physeq_ba)*100,3), Mean_Rb = round(rowMeans(physeq_rb)*100,3),
                          Mean_R3E = round(rowMeans(physeq_r3e)*100,3), Mean_R3W = round(rowMeans(physeq_r3w)*100,3),
                          Mean_tap = round(rowMeans(physeq_tap)*100,3), Mean_TN = round(rowMeans(physeq_tn)*100,3)))

fwrite(Diet.pct, file="Phylum_percent_abundances.txt", sep = "\t")
