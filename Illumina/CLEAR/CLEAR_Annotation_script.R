library(dplyr)
library(tidyr)

#setwd as location of quant files
setwd("G:/My Drive/Postdoc/12. May 2023/CLEAR_output/other_samples")


#reading in quant.txt files from my Gdrive (scratch work) - I manually added header in Excel to each file
Lc1 <- read.table("Lc1_quant.txt", header=T)
Lc3 <- read.table("Lc3_quant.txt", header=T)
Lc4 <- read.table("Lc4_quant.txt", header=T)


#add annotation
setwd("G:/My Drive/Postdoc/11. April 2023/")
gff <- read.delim("Lotusjaponicus_Gifu_v1.2_allPredictedGenes.gff3", skip=10, header=F) #reading in .gff file
gff <- separate(data=gff,col=V9,into=c("isoformId","Product"),sep=";human_readable_description=") #creating separate columns
gff$isoformId <- gsub("\\;.*", "", sub("ID=", "", gff$isoformId)) #removing extra info (everything after ; and 'ID=')
#remove repeated products
gffProduct <- gff %>% group_by(isoformId) %>% summarize(Product=toString(sort(unique(Product)))) #creating new dataframe with all product info for each gene
#add new column to quant files to match to gffProduct

Lc1$Product<-gffProduct[match(Lc1$isoformName,gffProduct$isoformId),]$Product
Lc3$Product<-gffProduct[match(Lc3$isoformName,gffProduct$isoformId),]$Product
Lc4$Product<-gffProduct[match(Lc4$isoformName,gffProduct$isoformId),]$Product


#Change wd back to folder with quant output
setwd("G:/My Drive/Postdoc/12. May 2023/CLEAR_output/other_samples")

#output annotated files
write.csv(Lc1, file = "Lc1_quant_annot.csv", row.names = FALSE)
write.csv(Lc3, file = "Lc3_quant_annot.csv", row.names = FALSE)
write.csv(Lc4, file = "Lc4_quant_annot.csv", row.names = FALSE)

