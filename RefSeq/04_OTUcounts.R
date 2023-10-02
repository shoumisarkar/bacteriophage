rm(list = ls())

## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
## Create the abundance table
## merge results from individual mapstat files
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
library(tidyverse)
library(data.table)
library(metagenomeSeq)

mapWD <- "path/to/RefSeq/phage/map_results"
taxonWD <- "path/to/RefSeq/phage/taxon_names"
setwd(mapWD)
mfiles <- dir()

# get the mapstat files
mfiles <- mfiles[which(sapply(strsplit(mfiles,"[.]"), "[[", 2, simplify=TRUE) == "mapstat")]

i=1
for(i in 1:length(mfiles)){
  # get the sample names
  sample_name <- sapply(strsplit(mfiles[i], "_"), function(x) paste(x[1], "_", x[2], sep=""), simplify = TRUE)

  # read the mapstat files. Use  `quote = "" ` to ignore quotes
  temp1 <- read.table(file=mfiles[i], sep="\t", skip=6, header=T, quote = "")
  
  # extract the accession numbers
  accs <- sapply(strsplit(sapply(strsplit(sapply(strsplit(temp1[,"refSequence"], " "), "[[", 1), "[|]"), "[[", 2), "_"), 
                 function(x) paste(x[1], "_", x[2], sep=""), simplify=TRUE)
  
  # remove the refSeq headers
  temp1 <- temp1[, -which(colnames(temp1) == c("refSequence"))]
  
  # merge the data
  temp1 <- data.frame(sample_name=sample_name, accs=accs, temp1)
  
  # collapse multiple accession numbers by sample ids
  temp1 <- temp1 %>% 
    group_by(sample_name, accs) %>% 
    summarise(read_counts = sum(readCount))
  
  # change to data.table
  temp1 <- as.data.table(temp1)
  
  # transform from long to wide
  temp1 <- dcast(temp1, accs ~ sample_name, fun.aggregate = NULL, value.var = "read_counts")
  setkey(temp1, accs) # set the key for DT
  
  # merge the OTU counts
  if(i > 1){
    temp <- merge(temp, temp1, all = TRUE) # merge by accession number
  } else {
    temp <- temp1
  }
}
rm(temp1)
gc()

# read the file containing accession numbers, NCBI taxon IDs and corresponding specie names
species <- read.table(paste(taxonWD, "/accs_taxID_species.txt", sep=""), header=T, sep="\t")
if(nrow(temp) != nrow(species)){
  stop("OTU table should have the same number of rows as the species identifiers")
}

## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
## Begin pre-processing
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
dat <- as.data.frame(temp)
row.names(dat) <- dat$accs
dat <- dat[, -which(colnames(dat) == "accs")]
dat[is.na(dat)] <- 0 # replace NA with 0

dim(dat)

## a) Remove OTUs with fewer than 10 reads
r1 <- which(rowSums(dat) < 10)
length(r1)
dat <- dat[-r1, ]

## b) remove OTUs which were present in fewer than 5% of samples
r2 <- which(rowSums(dat > 0) < ncol(dat)*0.5)
length(r2)
dat <- dat[-r2, ]

dat <- data.frame(taxa=rownames(dat),dat)
row.names(dat) <- NULL

# saves preprocessed OTU matrix as .csv file
write.csv(dat, file = paste(mapWD, "/OTU_counts_preproc.csv", sep=""), row.names = FALSE)


## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
## Begin normalization
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
## loads the OTU counts
dat2 <- loadMeta(file=paste(mapWD, "/OTU_counts_preproc.csv", sep=""), sep=",")
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

# saves preprocessed OTU matrix as .csv file
write.csv(dat2, file = paste(mapWD, "/OTU_counts_norm.csv", sep=""), row.names = TRUE)