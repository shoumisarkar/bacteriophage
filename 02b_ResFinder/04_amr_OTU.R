rm(list = ls())

## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
## Extract accession numbers of reference sequences for the mapped reads
## the files have the extension ".mapstat"
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##

library(tidyverse)
library(metagenomeSeq)
mapWD <- paste("/path/to/ResFinder/map_results", sep="")
setwd(mapWD)
mfiles <- dir()

# get the mapstat files
mfiles <- mfiles[which(sapply(strsplit(mfiles,"[.]"), "[[", 2, simplify=TRUE) == "mapstat")]

dat <- NULL
i=1
for(i in 1:length(mfiles)){

  # read the mapstat files. Use  `quote = "" ` to ignore quotes
  temp <- read.table(file=mfiles[i], sep="\t", skip=6, header=T, quote = "")
  
  if(nrow(temp) == 0){
    next
  }
  
  # get the sample names
  sample_name <- sapply(strsplit(mfiles[i], "_"), function(x) paste(x[1], "_", x[2], sep=""),
                        simplify = TRUE)
  
  # get the class/DB for the AMR
  con <- file(mfiles[i], "r")
  line <- as.character(readLines(con, 3)[3])
  close(con)
  amr_class <- sapply(strsplit(line, "\t"), "[[", 2)
  
  
  # merge the data
  temp <- data.frame(sample_name=sample_name, amr_class=amr_class, temp)
  
  # combine with other data
  dat <- rbind(temp, dat)
}
rm(temp)
gc()

## collate results by AMR class
temp <- dat %>% 
  group_by(sample_name, amr_class) %>% 
  summarise(amr_counts = sum(readCount))

# count of classes of AMR present in each sample
temp2 <- temp %>% 
  group_by(sample_name) %>% 
  summarise(amr_counts = length(amr_class))

## wide format
temp <- spread(temp, key=amr_class, value=amr_counts)
temp[is.na(temp)] <- 0 # replace NA with 0

temp <- as.data.frame(temp)
rownames(temp) <- temp$sample_name
temp <- temp[, -which(colnames(temp) == "sample_name")]
temp <- t(temp)

temp <- data.frame(taxa=rownames(temp),temp)
row.names(temp) <- NULL

# saves preprocessed AMR matrix as .csv file
write.csv(temp, file = paste(mapWD, "/AMR_counts_preproc.csv", sep=""), row.names = FALSE)


## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
## Begin normalization
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
## loads the OTU counts
dat2 <- loadMeta(file=paste(mapWD, "/AMR_counts_preproc.csv", sep=""), sep=",")
dim(dat2$counts)
dim(dat2$taxa)

# creates a MRexperiment object
dat2_obj <- newMRexperiment(dat2$counts)

# computes the percentile by which to normalize the counts
p <- cumNormStatFast(dat2_obj)

# calculates the normalization factors
dat2_norm <- cumNorm(dat2_obj, p=p)
dat2_mat <- MRcounts(dat2_norm, norm = TRUE, log = TRUE)
colSums(dat2_mat)

## Transpose the data matrix
dat2 <- t(dat2_mat)

# saves preprocessed AMR matrix as .csv file
write.csv(dat2, file = paste(mapWD, "/AMR_counts_norm.csv", sep=""), row.names = TRUE)
