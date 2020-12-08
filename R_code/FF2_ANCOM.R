library(exactRankTests)
library(nlme)
library(dplyr)
library(ggplot2)
library(compositions)
library(phyloseq)
library(data.table)
library(tidyverse)
library(dplyr)
library(magrittr)

# ANCOM 2.1 from: https://github.com/FrederickHuangLin/ANCOM/blob/master/scripts/ancom_v2.1.R
# Place in working directory
source("ancom_v2.1.R")

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

#setup variables for ANCOM
sample_var <- c("SeqNum")
lib_cut = 1000
neg_lb=FALSE   
main_var <- c("Treatment")
treats <- c("Banana", "CornStarch", "HAM2", "HAM4", "Potato", "PS_Ba", "PS_Rb", "RS3_Extracted",
                   "RS3_Whole", "Tapioca", "TigerNut")

for (fac in treats){
  #set pair for comparison, change here and output file name
  physeq_test <- subset_samples(physeq_51_5, Treatment %in% c("Water", fac))
  
  feature_table <- as.data.frame(otu_table(physeq_test))
  meta_data <- as.data.table(data.frame(sample_data(physeq_test)), keep.rownames = "SeqNum")
  
  # Data Pre-Processing
  processed_data <- feature_table_pre_process(feature_table, meta_data, sample_var, group_var=NULL, 
                                              out_cut = 0.05, zero_cut = 0.95, lib_cut, neg_lb)
  
  #Run ANCOM - could not include random formula for inoculum or else it would crash
  #Didn't seem to be a memory issue, but not sure
  ANCOM_results <- ANCOM(processed_data$feature_table, processed_data$meta_data, 
                         struc_zero = NULL, main_var, p_adj_method = "BH", 
                         alpha = 0.05, adj_formula = NULL, rand_formula = NULL)
  
  #Process results
  results_51_5 <- data.table(ANCOM_results$out)
  results_0.6 <- results_51_5[detected_0.6==TRUE,]
  result_dat <- ANCOM_results$fig$data
  
  theresults <- filter(result_dat, taxa_id %in% results_0.6$taxa_id)
  theresults <- cbind(theresults, results_0.6[,detected_0.9:detected_0.6])
  theresults <- dplyr::rename(theresults, meanCLR=x)
  theresults <- dplyr::rename(theresults, W_score=y)
  
  spec_nam <- as.data.table(as.data.frame(tax_table(physeq_test)), keep.rownames = "OTU")
  theresults <- cbind(taxa_id=theresults$taxa_id, Species=spec_nam[spec_nam$OTU %in% theresults$taxa_id, Species], 
                      theresults[,c(2:8)])
  
  fn <- paste("ANCOM_", fac, ".sig", sep = "")
  fwrite(theresults, file = fn, sep = "\t") 
}

