library(data.table)
library(ggplot2)
library(ggsci)

#Import relative abundance files of RS degraders and top 10 but producers
ButProds.inoc <- fread("Inoculum_but_producers_10.txt")
RSDegraders.inoc <- fread("Inoculum_RS_degraders.txt")
ButProds.trt <- fread("Diet_but_producers_10.txt")
RSDegraders.trt <- fread("Diet_RS_degraders.txt")

#Make joined names with OTU and species
ButProds.inoc <- cbind("OTU" = paste(ButProds.inoc$OTU, ButProds.inoc$Species, sep = "_"), 
                       ButProds.inoc[,-c("OTU", "Species")])
RSDegraders.inoc <- cbind("OTU" = paste(RSDegraders.inoc$OTU, RSDegraders.inoc$Species, sep = "_"), 
                          RSDegraders.inoc[,-c("OTU","Species")])

ButProds.trt <- cbind("OTU" = paste(ButProds.trt$OTU, ButProds.trt$Species, sep = "_"), 
                      ButProds.trt[,-c("OTU", "Species")])
RSDegraders.trt <- cbind("OTU" = paste(RSDegraders.trt$OTU, RSDegraders.trt$Species, sep = "_"), 
                          RSDegraders.trt[,-c("OTU","Species")])

#convert to long format for plotting
but_plot.inoc <- data.table::melt(ButProds.inoc, variable.name = "Inoculum", id.vars = "OTU", value.name = "Rel_Abun")
degrader_plot.inoc <- data.table::melt(RSDegraders.inoc, variable.name = "Inoculum", id.vars = "OTU", value.name = "Rel_Abun")

but_plot.trt <- data.table::melt(ButProds.trt, variable.name = "Treatment", id.vars = "OTU", value.name = "Rel_Abun")
degrader_plot.trt <- data.table::melt(RSDegraders.trt, variable.name = "Treatment", id.vars = "OTU", value.name = "Rel_Abun")

##Stacked bar chart##

Treatment_lab <- c("Wat", "Ap", "Bn", "CS", "HAM2", "HAM4", "PS", "PS_Ba", "PS_Rb", "R3E", "R3W","Tap","TN" )
Inoc_lab <- c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11")

tiff("But_inoc_rel_OTU.tiff", units="in", width=7, height=5, res=300)
ggplot(but_plot.inoc, aes(fill=OTU, y=Rel_Abun, x=Inoculum)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Butyrate Producer Relative Abundances by Inoculum") +
  theme (plot.title = element_text(hjust = 0.5)) +
  ylab("Percent Relative Abundance") +
  scale_x_discrete(labels = Inoc_lab) +
  coord_cartesian(ylim=c(0,30)) +
  scale_fill_d3("category10")
dev.off()

tiff("But_trt_rel_OTU.tiff", units="in", width=8, height=5, res=300)
ggplot(but_plot.trt, aes(fill=OTU, y=Rel_Abun, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Butyrate Producer Relative Abundances by Treatment") +
  theme (plot.title = element_text(hjust = 0.5), 
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Percent Relative Abundance") +
  scale_x_discrete(labels = Treatment_lab) +
  coord_cartesian(ylim=c(0,25)) +
  scale_fill_d3("category10")
dev.off()

tiff("RSdeg_inoc_rel_OTU.tiff", units="in", width=7, height=5, res=300)
ggplot(degrader_plot.inoc, aes(fill=OTU, y=Rel_Abun, x=Inoculum)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("RS Degrader Relative Abundances by Inoculum") +
  theme (plot.title = element_text(hjust = 0.5)) +
  ylab("Percent Relative Abundance") +
  scale_x_discrete(labels = Inoc_lab) +
  coord_cartesian(ylim=c(0,16)) +
  scale_fill_d3("category10")
dev.off()

tiff("RSdeg_trt_rel_OTU.tiff", units="in", width=8, height=5, res=300)
ggplot(degrader_plot.trt, aes(fill=OTU, y=Rel_Abun, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("RS Degrader Relative Abundances by Treatment") +
  theme (plot.title = element_text(hjust = 0.5), 
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Percent Relative Abundance") +
  scale_x_discrete(labels = Treatment_lab) +
  coord_cartesian(ylim=c(0,25)) +
  scale_fill_d3("category10")
dev.off()