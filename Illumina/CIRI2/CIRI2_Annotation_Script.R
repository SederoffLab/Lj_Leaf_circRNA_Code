#=====================================================================================
#
#  Add annotations to circRNA parent genes. Here we combine the CIRI2 outputs and the lotus annotated genome gff3
#
#=====================================================================================
library(ggplot2)
library(dplyr)
library(tidyr)

# Create merged file with script from DESeq2 code

# Set working directory to folder containing CIRI2 output .ciri text files
# change to location of .ciri files 
setwd("~/Desktop/1_yr_phd_current/Lab Notebook/inroot/Commensal Experiment/RNAseq analysis/circRNA analysis /CIRI2 results/CIRI_out/ciri as csv from excel")

# check the filename structure before creating your sample list
list.files(path="~/Desktop/1_yr_phd_current/Lab Notebook/inroot/Commensal Experiment/RNAseq analysis/circRNA analysis /CIRI2 results/CIRI_out/ciri as csv from excel")
#[1] "Lc_1.csv" "Lc_3.csv" "Lc_4.csv" 

# Create a list of filenames from above and put them in the correct order by treatment, which you want to have it in for the created file 
# change file names to match your raw counts files
files <- c("Lc_1.csv","Lc_3.csv","Lc_4.csv")

# Merge the first two files and store
file1<-read.csv(files[1], header=T)
file2<-read.csv(files[2], header=T)
out.cirifile = rbind(file1, file2) #gene_id is the parent gene
# View first 6 lines with:
head(out.cirifile)
# confirm that the correct columns (transcriptId and counts are present)

# For loop to merge the remaining files
for(i in 3:length(files))
{
  file = read.csv(files[i], header=T)
  out.cirifile <- rbind(out.cirifile, file)
}
View(out.cirifile)

#separe two gene ids in one column; also removes commas
out.cirifile<- separate(data=out.cirifile,col="gene_id",into=c("gene_id","gene_id2"),sep=",") #creating separate columns



# Prepare other dataframe, the .gff3 lotus genome annotation file used to bring in annotations 
setwd("~/Downloads") #use location of your .gff file
gff <- read.delim("Lotusgifuv1p2.gff3", skip=8, header=F) #reading in your .gff file
gff <- separate(data=gff,col=V9,into=c("GeneId","Product"),sep=";human_readable_description=") #creating separate columns
gff$GeneId <- gsub("\\..*||\\:.*||\\;.*", "", sub("ID=", "", gff$GeneId)) #removing extra text
gffProduct <- gff %>% group_by(GeneId) %>% summarize(Product=toString(sort(unique(Product)))) #creating new dataframe with all product info for each gene
#old code gff$Product[is.na(gff$Product)] <- 0 #replacing NAs in Product column to avoid data loss
#above line works but is not helpful in avoiding data loss
out.cirifile$Product<-NA
out.cirifile$Product<-gffProduct[match(out.cirifile$gene_id,gffProduct$GeneId),]$Product
out.cirifile$Product2<-NA
out.cirifile$Product2<-gffProduct[match(out.cirifile$gene_id2,gffProduct$GeneId),]$Product
#old code out.file$Product<-gff[match(out.file$GeneId,gff$GeneId),]$Product #still loses information!
View(out.cirifile) #new product column should have data but some may be missing (0's)

# save merged counts output file to see raw counts in excel
write.csv(out.cirifile, file = "~/Desktop/1_yr_phd_current/Lab Notebook/inroot/Commensal Experiment/RNAseq analysis/circRNA analysis /CIRI2 results/CIRI_out/annotciri.csv", row.names = FALSE)
