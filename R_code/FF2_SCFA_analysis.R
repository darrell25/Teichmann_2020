library(data.table)
library(dplyr)
library(ggpubr)
library(DescTools)
library(rstatix)
library(lme4)
library(multcomp)
library(ggplot2)

#This is the example for lactate, but a simple search and replace for any other of the fermentation
#acids will generate the equivalent plot for that one (e.g. replace lactate with butyrate)

Donors <- c("T01", "T02", "T03", "T04", "T05", "T06", "T07", "T08", "T09", "T10", "T11")
Treats_full <- c("Amylopectin", "Banana", "CornStarch", "HAM2", "HAM4", "PS", 
            "PS_Ba", "PS_Rb", "RS3_Extracted", "RS3_Whole", "Tapioca", "TigerNut", "Water")
Treats_abrev <- c("Ap", "Bn", "CS", "HAM2", "HAM4", "PS", "PS_Ba", "PS_Rb", "RS3_Ex", "RS3_Wh", "Tap", "Tn", "Water")

#read in the file that has all of the replicates (each fermentation was conducted in triplicate)
sample <- fread(file="FF2_FermSum.csv")
sample <- sample[order(Inoculum, Trt),]

sample$Trt <- factor(sample$Trt, levels = Treats_abrev, labels = Treats_abrev)
sample$Inoculum <- factor(sample$Inoculum, levels = Donors, labels = Donors)
sample$Trt <- relevel(as.factor(sample$Trt), ref = "Water")

sample.but <- sample[,c("ID", "Inoculum", "Trt", "Lactate")]

#means of replicates for each individual - avoid pseudoreplication when testing overall treatment effect
sample.means <- fread(file="FF2.meta.csv")
sample.means <- sample.means[order(Inoculum, Treatment),]
sample.means$Treatment <- factor(sample.means$Treatment, levels = Treats_full, labels = Treats_full)
sample.means$Inoculum <- factor(sample.means$Inoculum, levels = Donors, labels = Donors)
sample.means$Treatment <- relevel(as.factor(sample.means$Treatment), ref = "Water")

sample.means.but <- sample.means[,c("Group", "Inoculum", "Treatment", "Lactate")] 

####Overall Treatment Effects####

##initial look at the data##

ggboxplot(sample.means.but, x = "Treatment", y = "Lactate", add = "point")
sample.means.but %>% group_by(Treatment) %>% get_summary_stats(Lactate, type = "mean_sd")

##Checking Normality Assumptions##

sample.means.but %>% group_by(Treatment) %>% shapiro_test(Lactate)
sample.means.but %>% group_by(Treatment) %>% identify_outliers(Lactate)
ggqqplot(sample.means.but, "Lactate", facet.by = "Treatment")

#normality problems in the data, log transforming
sample.means.but.log <- sample.means.but
sample.means.but.log$Lactate <- log(sample.means.but$Lactate+1) #0s in data
sample.means.but.log %>% group_by(Treatment) %>% shapiro_test(Lactate)
ggqqplot(sample.means.but.log, "Lactate", facet.by = "Treatment")
#log transformation seems to fix

##Model Building and Hypothesis Testing##

#linear mixed model of log transformed lactate values to test for differences between treatment and control
fit.but <- lmerTest::lmer(Lactate ~ Treatment + (1 | Inoculum), data = sample.means.but.log)
But_test <- glht(fit.but, linfct = mcp(Treatment = "Dunnett"))
summary(But_test, test=adjusted("BH"))
capture.output(summary(But_test, test=adjusted("BH")), file="Treatment_lactate_results.txt", append=TRUE)

####Individual Treatment Effects####
sample.but.log <- sample.but
sample.but.log$Lactate <- log(sample.but$Lactate+1) #0s in data so use log(x+1)
sample.but.log %>% group_by(Trt) %>% shapiro_test(Lactate)
ggqqplot(sample.but.log, "Lactate", facet.by = "Trt")


sample.plot <- sample.means.but
sample.plot <- sample.plot[Treatment!="Water",]
sample.plot$Water <- numeric()
sample.plot$qval <- numeric()

#identify significant differences in lactate production between treatments and control
for (inoc in Donors){
  fit.but <- aov(Lactate ~ Trt, data = sample.but.log[Inoculum==inoc,])
  But_test <- glht(fit.but, linfct = mcp(Trt = "Dunnett"))
  write(paste("Lactate Changes for Inoculum ", inoc, sep = ""), file= "Individual_lactate_results.txt", append=TRUE)
  test <-summary(But_test, test=adjusted("BH"))
  sample.plot[Inoculum==inoc,"Water"] <- sample.means[Inoculum==inoc & Treatment=="Water","Lactate"]
  sample.plot[Inoculum==inoc,"qval"]<- test$test$pvalues 
  capture.output(test, file="Individual_lactate_results.txt", append=TRUE)
}
sample.plot$Diff <- sample.plot$Lactate - sample.plot$Water
sample.plot$Positive <- sample.plot$Diff > 0 #Identify pos vs neg changes
sample.plot$Abs_Diff <- abs(sample.plot$Diff) #find absolute difference for bubble size

#Set labels based on significance levels
sample.plot$Sig[sample.plot$qval<0.05] <- "*"
sample.plot$Sig[sample.plot$qval<0.01] <- "**"
sample.plot$Sig[sample.plot$qval<0.001] <- "***"

But_sig <- sample.plot
But_sig <- But_sig[is.na(But_sig$Sig)==FALSE, ] #get the significant changes
But_label <- data.frame(Inoculum=But_sig$Inoculum, Treatment=But_sig$Treatment)
Lab_txt <- c("Amylopectin", "Banana", "CornStarch", "HAM2", "HAM4", "Potato", 
            "PS_B.adol", "PS_R.brom", "RS3_Extracted", "RS3_Whole", "Tapioca_RS4", "TigerNut")

#make the bubble plot
tiff("Lactate_bubble_plot.tiff", units="in", width=10, height=6, res=300)
ggplot(sample.plot, aes(x=Inoculum, y=Treatment)) +
  geom_point(aes(size=Abs_Diff, color=Positive)) +
  scale_size(range=range(0,15),breaks=c(10,25,50,100), guide=guide_legend(), 
             name = "Lactate (mM)") +
  scale_color_discrete(name = "Change", labels = c("Negative", "Postive")) +
  labs(title = "Change in Lactate Production by Fecal Donor",
       x = "Fecal Donor",
       y = "Starch Source") +
  theme (plot.title = element_text(size = 18, hjust = 0.5), axis.text.x = element_text(size = 12),
         axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14),
         axis.title.y = element_text(size = 14), legend.title = element_text(size = 14),
         legend.text=element_text(size = 12)) +
  guides(color = guide_legend(order = 1), size = guide_legend(order=2)) +
  scale_y_discrete(label=Lab_txt) +
  geom_text (data = But_label, label = But_sig$Sig)
dev.off()  




