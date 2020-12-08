library(phyloseq)
library(DivNet)
library(breakaway)
library(lme4)
library(dplyr)
library(magrittr)
library(picante)
library(data.table)
library(ggsignif)
library(ggplot2)
library(multcomp)
library(tibble)
library(pairwiseAdonis)
library(DescTools)
library(mctoolsr)

####Data Import and Pre-processing####

#OTU and Tax tables come from previous run of Phyloseq::tax_glom, condensing to genus level
#Sample data is summary of metadata with mothur group names as row names
#tree file is from fast tree
OTU.matrix <- data.matrix(read.table(file="FF2_mothur_Transposed_gen.shared", header=TRUE, row.names=1, sep="\t"))
tax.matrix <- as.matrix(read.table(file="FF2_RDP_Blast_gen.taxonomy", header=TRUE, row.names=1, sep="\t"))
sample <- read.csv(file="FF2.meta.csv", header=TRUE, row.names=1, sep=",")

#Create phyloseq object
OTU.p <- otu_table(OTU.matrix, taxa_are_rows =TRUE)
TAX.p <- tax_table(tax.matrix)
sample.p <- sample_data(sample)
tree.p <- read_tree(treefile="FF2.tree")
physeq_gen <- phyloseq(OTU.p,TAX.p,sample.p, tree.p)

####Alpha-diversity testing####

##Run DivNet diversity tests##
a_div <- physeq_gen %>% divnet(ncores = 6)
estimates <- a_div$shannon %>% summary %$% estimate
ses <- a_div$shannon %>% summary %$% error
sample2 <- sample
sample2$Treatment <- relevel(as.factor(sample2$Treatment), ref="Water")

#calculate significance accounting for repeated measures
set.seed(3488)
a_div_test <- betta_random(chats = estimates,
                          ses = ses,
                          X = model.matrix(~Treatment, data = sample2),
                          groups=sample2$Inoculum)

#correct p-values via fdr
adjusted.shannon <- as.data.frame(a_div_test$table)
adjusted.shannon <- dplyr::rename(adjusted.shannon, pvalues='p-values')
adjusted.shannon$qvalues<- round(p.adjust(adjusted.shannon$pvalues, "BH"),3)
adjusted.shannon$Sig <- NA
adjusted.shannon$Sig[adjusted.shannon$qvalues<0.05] <- "*"
adjusted.shannon$Sig[adjusted.shannon$qvalues<0.01] <- "**"
adjusted.shannon$Sig[adjusted.shannon$qvalues<0.001] <- "***"

write ("\nGenus Level DivNet Treatment Shannon Test \n", file="ANOVA_results.txt", append=TRUE)
capture.output(adjusted.shannon, file="ANOVA_results.txt", append=TRUE)
adjusted.shannon
sample2$Shannon<-estimates

#Repeat procedure for simpson diversity
estimates2 <- a_div$simpson %>% summary %$% estimate
ses2 <- a_div$simpson %>% summary %$% error
a_div_test2 <- betta_random(chats = estimates2,
                           ses = ses2,
                           X = model.matrix(~Treatment, data = sample2),
                           groups=sample$Inoculum)

adjusted.simpson <- as.data.frame(a_div_test2$table)
adjusted.simpson <- dplyr::rename(adjusted.simpson, pvalues='p-values')
adjusted.simpson$qvalues<- round(p.adjust(adjusted.simpson$pvalues, "BH"),3)
adjusted.simpson$Sig <- NA
adjusted.simpson$Sig[adjusted.simposon$qvalues<0.05] <- "*"
adjusted.simpson$Sig[adjusted.simpson$qvalues<0.01] <- "**"
adjusted.simpson$Sig[adjusted.simpson$qvalues<0.001] <- "***"

write ("\nGenus Level DivNet Treatment Simmpson Test \n", file="ANOVA_results.txt", append=TRUE)
capture.output(adjusted.simpson, file="ANOVA_results.txt", append=TRUE)
adjusted.simpson
sample2$Simpson <- estimates2

##Faith's Phylogenetic Diversity##
set.seed(3488)
physeq_gen.rar <- rarefy_even_depth(physeq_gen, replace = FALSE)
OTU_gen.rar <- as.data.table(otu_table(physeq_gen.rar), keep.rownames = "OTU")
OTU_gen.rar.t <- data.table::transpose(OTU_gen.rar, keep.names = "Group", make.names = "OTU")
Faith_d.rar <- picante::pd(OTU_gen.rar.t, tree.p, include.root = FALSE)
sample2$Faith <- Faith_d.rar$PD

#Fit model with random effects
fit.faith <- lmerTest::lmer(Faith ~ Treatment + (1 | Inoculum), data = sample2)
faith_test <- glht(fit.faith, linfct = mcp(Treatment = "Dunnett"))

write("Genus Level lmer Treatment and Faith's Diversity Test \n", file= "ANOVA_results.txt", append=TRUE)
capture.output(summary(faith_test, test=adjusted("BH")), file="ANOVA_results.txt", append=TRUE)

#Create object to aid plotting
faith_sum <- summary(faith_test, test=adjusted("BH"))
adjusted.faith <- data.frame(qvalues=faith_sum$test$pvalues) 
adjusted.faith$Sig<- NA
adjusted.faith$Sig[adjusted.faith$qvalues<0.05] <- "*"
adjusted.faith$Sig[adjusted.faith$qvalues<0.01] <- "**"
adjusted.faith$Sig[adjusted.faith$qvalues<0.001] <- "***"

#use inverse simpson for actual plotting
inv_simp <- 1/estimates2
sample2$InvSimp <- inv_simp

#Create labels for plots and data frames
Treatments_abrev <- c("Water", "Ap", "Bn", "CS", "HAM2", "HAM4", "PS", "PS_Ba", "PS_Rb", "RS3_Ex", "RS3_Wh", "TapRS4", "TN")
Treatments_full <- c("Water", "Amylopectin", "Banana", "CornStarch", "HAM2", "HAM4", "Potato", "PS_Ba", "PS_Rb", "RS3_Extracted", "RS3_Whole", "Tapioca","TigerNut")
Treatments_test_full <- Treatments_full[-1]
Treatments_test_abrev <- Treatments_abrev[-1]

#Create objects for plotting significance of Shannon values
adjusted.shannon$max_shannon <- tapply(sample2$Shannon, sample2$Treatment, max)
adjusted.shannon$coord <- adjusted.shannon$max_shannon +0.25
Shannon_sig <- adjusted.shannon[-1,]
rownames(Shannon_sig) <- Treatments_test_full
Shannon_sig <- Shannon_sig[is.na(Shannon_sig$Sig)==FALSE, ]
Shan_label <- data.frame(Treatment = rownames(Shannon_sig), Shannon = Shannon_sig$coord)

#Create objects for plotting significance of Simpson values
adjusted.simpson$max_simpson <- tapply(sample2$Simpson, sample2$Treatment, max)
adjusted.simpson$coord <- adjusted.simpson$max_simpson +0.1
Simpson_sig <- adjusted.simpson[-1,]
rownames(Simpson_sig) <- Treatments_test_full
Simpson_sig <- Simpson_sig[is.na(Simpson_sig$Sig)==FALSE, ]
Simp_label <- data.frame(Treatment = rownames(Simpson_sig), Simpson = Simpson_sig$coord)

#Create objects for plotting Inverse Simpson
adjusted.invsimp <- data.frame(qvalues = adjusted.simpson$qvalues, Sig = adjusted.simpson$Sig, row.names = Treatments_full)
adjusted.invsimp$max_invsimpson <- tapply(sample2$InvSimp, sample2$Treatment, max)
adjusted.invsimp$coord <- adjusted.invsimp$max_invsimpson + 1
InvSimpson_sig <- adjusted.invsimp[-1,]
InvSimpson_sig <- InvSimpson_sig[is.na(InvSimpson_sig$Sig)==FALSE, ]
InvSimp_label <- data.frame(Treatment = rownames(InvSimpson_sig), invSimp = InvSimpson_sig$coord)

#Create objects for plotting Faith's Diversity
sample3 <- data.frame(Treatment = as.character(sample2$Treatment), Faith = sample2$Faith, row.names = rownames(sample2))
sample3 <- sample3[sample3$Treatment != "Water",]
adjusted.faith$max_faith <- tapply(sample3$Faith, sample3$Treatment, max)
adjusted.faith$coord <- adjusted.faith$max_faith +1
Faith_sig <- adjusted.faith
rownames(Faith_sig) <- Treatments_test_full
Faith_sig <- Faith_sig[is.na(Faith_sig$Sig)==FALSE, ]
Faith_label <- data.frame(Treatment = rownames(Faith_sig), Faith = Faith_sig$coord)

plot_sample <- data.frame(Group=rownames(sample2), Treatment=sample2$Treatment, 
                          Shannon=sample2$Shannon, Simpson=sample2$Simpson, invSimp=sample2$InvSimp, Faith=sample2$Faith)


tiff("Genus_Shannon_Diet.tiff", units="in", width=5, height=5, res=300)
ggplot(plot_sample, aes(Treatment, Shannon)) + 
  geom_boxplot() + 
  ggtitle("Genus-level Shannon Diversity") + 
  theme (plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(label=Treatments_abrev) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text (data = Shan_label, label = Shannon_sig$Sig)
dev.off()

tiff("Genus_Simpson_Diet.tiff", units="in", width=5, height=5, res=300)
ggplot(plot_sample, aes(Treatment, Simpson)) + 
  geom_boxplot() +
  ggtitle("Genus-level Simpson Diversity") +
  theme (plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(label=Treatments_abrev) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text (data = Simp_label, label = Simpson_sig$Sig)
dev.off()

tiff("Genus_Inverse_Simpson_Diet.tiff", units="in", width=5, height=5, res=300)
ggplot(plot_sample, aes(Treatment, invSimp)) + 
  geom_boxplot() +
  ggtitle("Genus-level Inverse Simpson Diversity") +
  theme (plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(label=Treatments_abrev) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text (data = InvSimp_label, label = InvSimpson_sig$Sig)
  ylab("Inverse Simpson")
dev.off()

tiff("Genus_Faith_Diet.tiff", units="in", width=5, height=5, res=300)
ggplot(plot_sample, aes(Treatment, Faith)) + 
  geom_boxplot() +
  ggtitle("Genus-level Faith Diversity") +
  theme (plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(label=Treatments_abrev) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text (data = Faith_label, label = Faith_sig$Sig)
dev.off()

####Beta-diversity Analysis####
sample.a <- as.data.frame(sample2)

##Bray-Curtis Distance via DivNet##

#Create distance object and perform ordination
b_div_bray <- as.dist(a_div$`bray-curtis`)
b_div_bray_pcoa <- ordinate(physeq_gen, method="PCoA", distance=b_div_bray) 


tiff("DivNet_Genus_Bray_Diet.tiff", units="in", width=5, height=5, res=300)
phyloseq::plot_ordination(physeq_gen, b_div_bray_pcoa, type ="samples", color = "Treatment") +
  ggtitle("Bray-Curtis Dissimilarity") +
  theme (plot.title = element_text(hjust = 0.5))
dev.off()

#Test for overall significance
write("Genus Level DivNet Diet and Bray-Curtis Diversity Test", 
      file= "PERMANOVA_results_overall.txt", append=TRUE)
set.seed(3488)
capture.output(adonis(formula = b_div_bray ~ Treatment, data = sample.a), 
               file="PERMANOVA_results_overall.txt", append=TRUE)

#Test for pair-wise significance
set.seed(3488)
adjusted.bray <- calc_pairwise_permanovas(b_div_bray, sample.a, "Treatment")
adjusted.bray <- adjusted.bray[adjusted.bray$X1=="Water",-c(5:6)]
adjusted.bray$qval <- p.adjust(adjusted.bray$pval, method = "BH")
adjusted.bray$Sig <- NA
adjusted.bray$Sig[adjusted.bray$qval<0.05] <- "*"
adjusted.bray$Sig[adjusted.bray$qval<0.01] <- "**"
adjusted.bray$Sig[adjusted.bray$qval<0.001] <- "***"

write("Genus Level DivNet Treatment and Bray-Curtis Diversity Test", file= "PERMANOVA_results.txt", append=TRUE)
capture.output(adjusted.bray, file="PERMANOVA_results.txt", append=TRUE)

##Aitchison Distance Analysis##

#CLR transformation of counts
physeq_gen_clr <- microbiome::transform(physeq_gen, "clr")
otu.clr <- t(as(otu_table(physeq_gen_clr), "matrix"))
#calculate aitchison distance (Euclidan distance of clr transformed data)
clr_ait <- dist(otu.clr, method='euc')
b_div_ait_pcoa <- ordinate(physeq_gen, method="PCoA", distance=clr_ait)


tiff("Genus_Aitchison_Diet.tiff", units="in", width=5, height=5, res=300)
phyloseq::plot_ordination(physeq_gen, b_div_ait_pcoa, type ="samples", color = "Treatment") +
  ggtitle("Aitchison Distance") +
  theme (plot.title = element_text(hjust = 0.5))
dev.off()

#Test for overall significance
write("\nGenus Level Diet and Aitchison Diversity Test", 
      file= "PERMANOVA_results_overall.txt", append=TRUE)
set.seed(3488)
capture.output(adonis(formula = clr_ait ~ Treatment, data = sample.a), 
               file="PERMANOVA_results_overall.txt", append=TRUE)

#Test for pair-wise significance
set.seed(3488)
adjusted.ait <- calc_pairwise_permanovas(clr_ait, sample.a, "Treatment")
adjusted.ait <- adjusted.ait[adjusted.ait$X1=="Water",-c(5:6)]
adjusted.ait$qval <- p.adjust(adjusted.ait$pval, method = "BH")
adjusted.ait$Sig <- NA
adjusted.ait$Sig[adjusted.ait$qval<0.05] <- "*"
adjusted.ait$Sig[adjusted.ait$qval<0.01] <- "**"
adjusted.ait$Sig[adjusted.ait$qval<0.001] <- "***"

write("\nGenus Level Diet and Aitchison Diversity Test", file= "PERMANOVA_results.txt", append=TRUE)
capture.output(adjusted.ait, file="PERMANOVA_results.txt", append=TRUE)

##Unifrac analysis##

dist_gen_rar_uniw <- phyloseq::UniFrac(physeq_gen.rar, weighted = TRUE)
gen_rar_uniw_pcoa <- ordinate(physeq_gen.rar, method="PCoA", distance=dist_gen_rar_uniw)

tiff("Genus_UnifracW_Diet.tiff", units="in", width=5, height=5, res=300)
phyloseq::plot_ordination(physeq_gen.rar, gen_rar_uniw_pcoa, type ="samples", color = "Treatment") +
  ggtitle("Weighted UniFrac") +
  theme (plot.title = element_text(hjust = 0.5))
dev.off()

#Test for overall significance
write("\nGenus Level Diet and Weighted Unifrac Diversity Test", 
      file= "PERMANOVA_results_overall.txt", append=TRUE)
set.seed(3488)
capture.output(adonis(formula = dist_gen_rar_uniw ~ Treatment, data = sample.a), 
               file="PERMANOVA_results_overall.txt", append=TRUE)

#Test for pair-wise significance
set.seed(3488)
adjusted.uniw <- calc_pairwise_permanovas(dist_gen_rar_uniw, sample.a, "Treatment")
adjusted.uniw <- adjusted.uniw[adjusted.uniw$X1=="Water",-c(5:6)]
adjusted.uniw$qval <- p.adjust(adjusted.uniw$pval, method = "BH")
adjusted.uniw$Sig <- NA
adjusted.uniw$Sig[adjusted.uniw$qval<0.05] <- "*"
adjusted.uniw$Sig[adjusted.uniw$qval<0.01] <- "**"
adjusted.uniw$Sig[adjusted.uniw$qval<0.001] <- "***"

write("\nGenus Level Treatment and Weighted Unifrac Diversity Test", file= "PERMANOVA_results.txt", append=TRUE)
capture.output(adjusted.uniw, file="PERMANOVA_results.txt", append=TRUE)
