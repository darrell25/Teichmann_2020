library(data.table)
library(phyloseq)
library(microbiome)
library(microbiomeSeq)
library(magrittr)
library(dplyr)
library(data.table)

OTU.df <- data.frame(read.table(file="FF2_mothur_Transposed.shared", header=TRUE, row.names=1, sep="\t"))
tax.df <- data.frame(read.table(file="FF2_RDP_Blast.taxonomy", header=TRUE, row.names=1, sep="\t"))
sample <- fread("FF2.meta.csv")
#Input list of identified butyrate producers in the first 200 OTUs
ButProds <- fread("Diet_but_producers.txt")
#Input list of identified RS degraders in the first 200 OTUs
RSDegraders <- fread("Diet_RS_degraders.txt")

#Make the species names combining the genus, species and OTU# (rowname)
tax.df$Species <- make.names (ifelse(tax.df$Genus==tax.df$Species,paste(rownames(tax.df), tax.df$Genus, sep = "_"),
                                     paste(rownames(tax.df), tax.df$Genus, tax.df$Species, sep="_")), unique=FALSE)

RSBut <- rbind(ButProds, RSDegraders)

#Create categories of high butyrate and low butyrate producing substrates
Low_but <- c("Water", "Amylopectin", "CornStarch", "PS_Ba", "PS_Rb", "RS3_Extracted")
High_but <- c("Banana", "HAM2", "HAM4", "PS", "RS3_Whole", "Tapioca", "TigerNut")

sample$TrtCat <- c()
sample[Treatment %in% Low_but, "TrtCat"] <- "Low.But.Trt"
sample[Treatment %in% High_but, "TrtCat"] <- "High.But.Trt"

sample <- data.frame(sample, row.names = "Group")

##Create phyloseq object##
OTU.p <- otu_table(OTU.df, taxa_are_rows =TRUE)
TAX.p <- tax_table(as.matrix(tax.df))
sample.p <- sample_data(sample)
physeq <- phyloseq(OTU.p,TAX.p,sample.p)

#CLR transform counts
physeq_clr <- microbiome::transform(physeq, "clr")

#extract OTU table and make the RS and butyrate producer sub tables
OTU.clr <- data.frame(otu_table(physeq_clr))

OTU.rsbut <- OTU.clr[rownames(OTU.clr) %in% RSBut$OTU,]
tax.rsbut <- tax.df[rownames(tax.df) %in% RSBut$OTU,]

OTU.rsdeg <- as.data.table(OTU.clr[rownames(OTU.clr) %in% RSDegraders$OTU,], keep.rownames = "OTU")
OTU.rsdeg <- data.table::transpose(OTU.rsdeg, keep.names = "Group", make.names = "OTU")

sample <- cbind(sample, OTU.rsdeg[,c(2:5)])
sample <- dplyr::rename(sample, Otu00004_B_adolescentis = Otu00004, 
                        Otu00013_R_bromii = Otu00013)

rownames(OTU.rsbut) <- tax.rsbut$Species
rownames(tax.rsbut) <- tax.rsbut$Species

#recreate phyloseq object
OTU.p <- otu_table(OTU.rsbut, taxa_are_rows =TRUE)
TAX.p <- tax_table(as.matrix(tax.rsbut))
sample.p <- sample_data(sample)
physeq_but <- phyloseq(OTU.p,TAX.p,sample.p)

#Look at correlations between butyrate producers, RS degraders, butyrate, pH
mycorr <- taxa.env.correlation(physeq_but, grouping_column = "TrtCat", method="spearman", pvalue.threshold=0.05,
                                      padjust.method="BH", adjustment=5, num.taxa=29, 
                               select.variables = c("Otu00004_B_adolescentis", "Otu00013_R_bromii", "pH", "Butyrate"))

tiff("Corrleations_ButProd_RSDeg.tiff", units="in", width=10, height=6, res=300)
plot_taxa_env(mycorr)
dev.off()




