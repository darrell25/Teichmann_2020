library(data.table)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(magrittr)

Treats_full <- c("Amylopectin", "Banana", "CornStarch", "HAM2", "HAM4", "PS", 
                 "PS_Ba", "PS_Rb", "RS3_Extracted", "RS3_Whole", "Tapioca", "TigerNut", "Water")
Donors <- c("T01", "T02", "T03", "T04", "T05", "T06", "T07", "T08", "T09", "T10", "T11")

#means of replicates for each individual - avoid pseudoreplication 
sample.means <- fread(file="FF2.meta.csv")

#don't think this reordering and making factors matters for this analysis, but just to be safe
sample.means <- sample.means[order(Inoculum, Treatment),]
sample.means$Treatment <- factor(sample.means$Treatment, levels = Treats_full, labels = Treats_full)
sample.means$Inoculum <- factor(sample.means$Inoculum, levels = Donors, labels = Donors)
sample.means$Treatment <- relevel(as.factor(sample.means$Treatment), ref = "Water")

TestCorrs <- c("Acetate", "Formate", "Lactate", "Propionate", "Succinate", "pH")
restab <- data.table(SCFA = TestCorrs, Rho = numeric(), pval = numeric(), qval = numeric())

#Test correlations between each fermentation parameter and butyrate
for (tests in TestCorrs){
  cols <- c("Group", "Inoculum", "Treatment", tests, "Butyrate")
  sample.means.but <- sample.means[,..cols] 
  tiff(paste(tests, "_spearman_correlation_plot2.tiff", sep = ""), units="in", width=7, height=5, res=300)
  p <- ggscatter(sample.means.but, x = tests, y = "Butyrate", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = paste (tests, "Concentration (mM)", sep = " "), ylab = "Butyrate Concentration (mM)")
  plot(p)
  dev.off()
  sp_res <- cor.test(sample.means.but[,get(tests)], sample.means.but$Butyrate, method = "spearman")
  restab[SCFA==tests, 'Rho'] <- sp_res$estimate
  restab[SCFA==tests, 'pval'] <- sp_res$p.value
}
restab$qval<-p.adjust(restab$pval, method = "BH")

fwrite(restab, file = "Butyrate_SCFA_spearman_correlations.txt", sep = "\t")
