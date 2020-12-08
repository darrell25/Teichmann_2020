library(data.table)
library(phyloseq)
library(dplyr)
library(magrittr)

#Import data - OTU file comes from transpose_mothur_shared_file.R
#Taxonomy comes from AddSpeciesToTaxonomy.R
#Sample data is summary of metadata with mothur group names as row names
OTU.df <- data.frame(read.table(file="FF2_mothur_Transposed.shared", header=TRUE, row.names=1, sep="\t"))
tax.dt <- fread(file = "FF2_RDP_Blast.taxonomy", header = TRUE, sep="\t")
sample <- read.csv(file="FF2.meta.csv", header=TRUE, row.names=1, sep=",")

#Create joined genus_species names
tax.dt$Species <- make.names (ifelse(tax.dt$Genus == tax.dt$Species, paste(tax.dt$Species, tax.dt$OTU, sep = "_"), 
                                     paste(tax.dt$Genus, tax.dt$Species, tax.dt$OTU, sep="_")), unique=TRUE)

tax.matrix<-as.matrix(tax.dt, rownames = "OTU")
#Create phyloseq object
OTU.p <- otu_table(OTU.df, taxa_are_rows =TRUE)
TAX.p <- tax_table(tax.matrix)
sample.p <- sample_data(sample)
physeq <- phyloseq(OTU.p,TAX.p,sample.p)

##Create filter to have minimum 51 counts (0.001%) and 5% prevalence

physeq_51 <- filter_taxa(physeq, function(x) sum(x) >=51, prune=TRUE)
prevalenceThreshold <- 0.05 * nsamples(physeq_51)

# Compute prevalence of each feature, store as data.frame
prevdf51 <- apply(X = otu_table(physeq_51),
                  MARGIN = ifelse(taxa_are_rows(physeq_51), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf51 <- data.frame(Prevalence = prevdf51,
                       TotalAbundance = taxa_sums(physeq_51),
                       tax_table(physeq_51))

keepTaxa2 <- rownames(prevdf51)[(prevdf51$Prevalence >= prevalenceThreshold)]
physeq_51_5 <- prune_taxa(keepTaxa2, physeq_51)


physeq_dom <- tax_glom(physeq_51_5, taxrank = "Domain")
physeq_phy <- tax_glom(physeq_51_5, taxrank = "Phylum")
physeq_class <- tax_glom(physeq_51_5, taxrank = "Class")
physeq_order <- tax_glom(physeq_51_5, taxrank = "Order")
physeq_fam <- tax_glom(physeq_51_5, taxrank = "Family")
physeq_gen <- tax_glom(physeq_51_5, taxrank = "Genus")

#Transform to relative abundances
physeq_dom <- transform_sample_counts(physeq_dom, function(x) x/sum(x))
physeq_phy <- transform_sample_counts(physeq_phy, function(x) x/sum(x))
physeq_class <- transform_sample_counts(physeq_class, function(x) x/sum(x))
physeq_order <- transform_sample_counts(physeq_order, function(x) x/sum(x))
physeq_fam <- transform_sample_counts(physeq_fam, function(x) x/sum(x))
physeq_gen <- transform_sample_counts(physeq_gen, function(x) x/sum(x))
physeq_sp <- transform_sample_counts(physeq_51_5, function(x) x/sum(x))

#Extract the otu and tax tables
otu_dom <- as.data.table(otu_table(physeq_dom), keep.rownames = "OTU")
tax_dom <- as.data.table(as.data.frame(tax_table(physeq_dom)), keep.rownames = "OTU")
otu_phy <- as.data.table(otu_table(physeq_phy), keep.rownames = "OTU")
tax_phy <- as.data.table(as.data.frame(tax_table(physeq_phy)), keep.rownames = "OTU")
otu_class <- as.data.table(otu_table(physeq_class), keep.rownames = "OTU")
tax_class <- as.data.table(as.data.frame(tax_table(physeq_class)), keep.rownames = "OTU")
otu_order <- as.data.table(otu_table(physeq_order), keep.rownames = "OTU")
tax_order <- as.data.table(as.data.frame(tax_table(physeq_order)), keep.rownames = "OTU")
otu_fam <- as.data.table(otu_table(physeq_fam), keep.rownames = "OTU")
tax_fam <- as.data.table(as.data.frame(tax_table(physeq_fam)), keep.rownames = "OTU")
otu_gen <- as.data.table(otu_table(physeq_gen), keep.rownames = "OTU")
tax_gen <- as.data.table(as.data.frame(tax_table(physeq_gen)), keep.rownames = "OTU")
otu_sp <- as.data.table(otu_table(physeq_sp), keep.rownames = "OTU")
tax_sp <- as.data.table(as.data.frame(tax_table(physeq_sp)), keep.rownames = "OTU")

#combine all taxonomic levels together
Lef_file <- data.table(Sample = tax_dom$Domain)
Lef_file <- cbind(Lef_file, otu_dom[,OTU:=NULL])
phyla <- data.table(Sample = paste(tax_phy$Domain, tax_phy$Phylum, sep = "|"))
phyla <- cbind(phyla, otu_phy[,OTU:=NULL])
Lef_file <- rbind(Lef_file, phyla)
phyla <- data.table(Sample = paste(tax_class$Domain, tax_class$Phylum, tax_class$Class, sep = "|"))
phyla <- cbind(phyla, otu_class[,OTU:=NULL])
Lef_file <- rbind(Lef_file, phyla)
phyla <- data.table(Sample = paste(tax_order$Domain, tax_order$Phylum, tax_order$Class, tax_order$Order, sep = "|"))
phyla <- cbind(phyla, otu_order[,OTU:=NULL])
Lef_file <- rbind(Lef_file, phyla)
phyla <- data.table(Sample = paste(tax_fam$Domain, tax_fam$Phylum, tax_fam$Class, tax_fam$Order, tax_fam$Family, sep = "|"))
phyla <- cbind(phyla, otu_fam[,OTU:=NULL])
Lef_file <- rbind(Lef_file, phyla)
phyla <- data.table(Sample = paste(tax_gen$Domain, tax_gen$Phylum, tax_gen$Class, tax_gen$Order, tax_gen$Family, tax_gen$Genus, sep = "|"))
phyla <- cbind(phyla, otu_gen[,OTU:=NULL])
Lef_file <- rbind(Lef_file, phyla)
phyla <- data.table(Sample = paste(tax_sp$Domain, tax_sp$Phylum, tax_sp$Class, tax_sp$Order, tax_sp$Family, tax_sp$Genus, tax_sp$Species, sep = "|"))
phyla <- cbind(phyla, otu_sp[,OTU:=NULL])
Lef_file <- rbind(Lef_file, phyla)

sample_tb <- as.data.table(sample, keep.rownames = "Group")
sample_t <- data.table::transpose(sample_tb, make.names = "Group", keep.names = "Sample")
Lef_file <- rbind(sample_t[Sample %in% c("Inoculum", "Treatment"),], Lef_file)

#all Treatments
fwrite(Lef_file, file= "FF2_OTU_LEfSe.txt", col.name = FALSE, sep = "\t")

#make versions with just water and one other treatment
Lef_ap <- Lef_file  %>% select_if(Lef_file[2,] == "Water"|Lef_file[2,] == "Amylopectin")
Lef_ap <- cbind(Sample=Lef_file$Sample, Lef_ap)
fwrite(Lef_ap, file= "FF2_OTU_LEfSe_Ap.txt", col.name = FALSE, sep = "\t")

Lef_bn <- Lef_file  %>% select_if(Lef_file[2,] == "Water"|Lef_file[2,] == "Banana")
Lef_bn <- cbind(Sample=Lef_file$Sample, Lef_bn)
fwrite(Lef_bn, file= "FF2_OTU_LEfSe_Bn.txt", col.name = FALSE, sep = "\t")

Lef_cs <- Lef_file  %>% select_if(Lef_file[2,] == "Water"|Lef_file[2,] == "CornStarch")
Lef_cs <- cbind(Sample=Lef_file$Sample, Lef_cs)
fwrite(Lef_cs, file= "FF2_OTU_LEfSe_CS.txt", col.name = FALSE, sep = "\t")

Lef_ham2 <- Lef_file  %>% select_if(Lef_file[2,] == "Water"|Lef_file[2,] == "HAM2")
Lef_ham2 <- cbind(Sample=Lef_file$Sample, Lef_ham2)
fwrite(Lef_ham2, file= "FF2_OTU_LEfSe_HAM2.txt", col.name = FALSE, sep = "\t")

Lef_ham4 <- Lef_file  %>% select_if(Lef_file[2,] == "Water"|Lef_file[2,] == "HAM4")
Lef_ham4 <- cbind(Sample=Lef_file$Sample, Lef_ham4)
fwrite(Lef_ham4, file= "FF2_OTU_LEfSe_HAM4.txt", col.name = FALSE, sep = "\t")

Lef_pot <- Lef_file  %>% select_if(Lef_file[2,] == "Water"|Lef_file[2,] == "Potato")
Lef_pot <- cbind(Sample=Lef_file$Sample, Lef_pot)
fwrite(Lef_pot, file= "FF2_OTU_LEfSe_Pot.txt", col.name = FALSE, sep = "\t")

Lef_ba <- Lef_file  %>% select_if(Lef_file[2,] == "Water"|Lef_file[2,] == "PS_Ba")
Lef_ba <- cbind(Sample=Lef_file$Sample, Lef_ba)
fwrite(Lef_ba, file= "FF2_OTU_LEfSe_PS_Ba.txt", col.name = FALSE, sep = "\t")

Lef_rb <- Lef_file  %>% select_if(Lef_file[2,] == "Water"|Lef_file[2,] == "PS_Rb")
Lef_rb <- cbind(Sample=Lef_file$Sample, Lef_rb)
fwrite(Lef_rb, file= "FF2_OTU_LEfSe_PS_Rb.txt", col.name = FALSE, sep = "\t")

Lef_r3e <- Lef_file  %>% select_if(Lef_file[2,] == "Water"|Lef_file[2,] == "RS3_Extracted")
Lef_r3e <- cbind(Sample=Lef_file$Sample, Lef_r3e)
fwrite(Lef_r3e, file= "FF2_OTU_LEfSe_RS3_Ex.txt", col.name = FALSE, sep = "\t")

Lef_r3w <- Lef_file  %>% select_if(Lef_file[2,] == "Water"|Lef_file[2,] == "RS3_Whole")
Lef_r3w <- cbind(Sample=Lef_file$Sample, Lef_r3w)
fwrite(Lef_r3w, file= "FF2_OTU_LEfSe_RS3_Wh.txt", col.name = FALSE, sep = "\t")

Lef_tap <- Lef_file  %>% select_if(Lef_file[2,] == "Water"|Lef_file[2,] == "Tapioca")
Lef_tap <- cbind(Sample=Lef_file$Sample, Lef_tap)
fwrite(Lef_tap, file= "FF2_OTU_LEfSe_TapRS4.txt", col.name = FALSE, sep = "\t")

Lef_tn <- Lef_file  %>% select_if(Lef_file[2,] == "Water"|Lef_file[2,] == "TigerNut")
Lef_tn <- cbind(Sample=Lef_file$Sample, Lef_tn)
fwrite(Lef_tn, file= "FF2_OTU_LEfSe_TN.txt", col.name = FALSE, sep = "\t")
