MakeBlastFasta <- function(fafilepath){
  library(Biostrings)
  library(stringr)
  
  #read in the fasta file
  repFasta <- readDNAStringSet(filepath = fafilepath, format="fasta") 
  
  #get just the OTU number from this, basically the text between the start of "Otu" and the "|"
  names(repFasta)<- substring(names(repFasta), regexpr("Otu", names(repFasta)), regexpr("\\|", names(repFasta))-1)
  
  writeXStringSet(repFasta, filepath = 'FastaforBlast.fasta',format="fasta")
}

AddSpeciesToTaxonomy <- function (taxonomyfile, blastfile){
  #input files are your mothur generated taxonomy file for Domain through genus and Blastn file in following format:
  
  library(stringr)
  library(tidyverse)
  library(data.table)
  
  Tax <- read.table(taxonomyfile, sep="\t", header = TRUE)
  Blast <- read.table(blastfile, sep="\t")
  colnames(Blast)<-c("OTU", "Pident", "TaxID", "Species")
  
  #First need to sub in placeholders for missing values from Blastn
  Blastspecies<-data.frame("OTU"=rep(0,length(Tax$OTU)), "Genus" = rep(0,length(Tax$OTU)), "Species"= rep(0,length(Tax$OTU)))
  where<-match(Blast$OTU, Tax$OTU)
  Blastspecies[where,c("OTU","Species")]<-Blast[,c("OTU","Species")]
  
  #Make mothur formated taxonomy with species, including percent identity
  Blastspecies[where,"Species"]<-paste(sapply (strsplit(Blast$Species," "), "[",1),"_", sapply (strsplit(Blast$Species," "), "[",2),"(",round(Blast$Pident),");",sep="")
  Tax$Taxonomy<- ifelse(Tax$OTU==Blastspecies$OTU,paste(Tax$Taxonomy,Blastspecies$Species, sep=""),
                        paste(Tax$Taxonomy,sapply(strsplit(Tax$Taxonomy,";"),"[",6),sep=""))
  write.table(Tax, file="TaxTable.txt",append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE)
  
  #Generate phyloseq formatted file
  #First get just the species and genus parts from Blast separately, correct if only 1 part name
  Blastspecies[where,"Species"] <- paste(ifelse(grepl(" ", Blast$Species, fixed=TRUE), sapply (strsplit(Blast$Species," "), 
                                                                                               "[",2),sapply (strsplit(Blast$Species," "), "[",1)) ,"(",round(Blast$Pident),");",sep="")
  Blastspecies[where,"Genus"] <- sapply(strsplit(Blast$Species, " "), "[",1)
  
  #Then split up the taxonomy into individual columns
  PhylTax<- data.table(Tax$OTU, str_split_fixed(Tax$Taxonomy, ";", 7), "BlastGenus" = rep("0",length(Tax$OTU)))
  names(PhylTax)<-c("OTU", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "BlastGenus")
  
  #Now replace species level with just species (no genus), removing the trailing semicolon
  PhylTax$Species<-if_else(PhylTax$OTU==Blastspecies$OTU,substring(Blastspecies$Species,1,nchar(Blastspecies$Species)-1),PhylTax$Species)
  
  #Propagate back the genus of the Blast hit if the species has > 97% identity
  PhylTax[where,"BlastGenus"] <- Blastspecies[where, "Genus"]
  PhylTax[as.numeric(substring(Species,regexpr("\\(",Species)+1, regexpr("\\)",Species)-1)) >= 97, Genus := BlastGenus]
  PhylTax<- subset(PhylTax, select = -c(BlastGenus))
  
  #Propagate forward last confident ID from RDP classifer in the event BLAST species level ID <97% identical
  PhylTax[as.numeric(substring(Species,regexpr("\\(",Species)+1, regexpr("\\)",Species)-1)) < 97, Species := Genus]
  
  #this is to remove the bracketed numbers
  #there should be a better way to do this, but this is the only one I could get to work as expected
  
  PhylTax[,Domain:=str_replace_all(Domain,"\\(.*?\\)","")]
  PhylTax[,Phylum:=str_replace_all(Phylum,"\\(.*?\\)","")]
  PhylTax[,Class:=str_replace_all(Class,"\\(.*?\\)","")]
  PhylTax[,Order:=str_replace_all(Order,"\\(.*?\\)","")]
  PhylTax[,Family:=str_replace_all(Family,"\\(.*?\\)","")]
  PhylTax[,Genus:=str_replace_all(Genus,"\\(.*?\\)","")]
  PhylTax[,Species:=str_replace_all(Species,"\\(.*?\\)","")]
  
  write.table(PhylTax, file="PhylTaxTable.txt",append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE)
  
}

TransposeTable<-function(otutable, outfile){
  #it is more complicated to install data.table, see package instructions
  library(data.table)
  
  table1<-fread(file=otutable, header=T, sep='\t')
  #note that data.table doesn't use rownames, so you cannot specify
  
  #this removes the label and numOTUs columns, if you want to keep comment the following line
  table1<- subset(table1, select = -c(label, numOtus))
  
  table1.t<-data.table::transpose(table1, keep.names = "OTU",make.names = "Group")
  
  fwrite(table1.t, file=outfile, sep="\t")
  
}