rm(list = ls())

## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
## Extract accession numbers of reference sequences for the mapped reads
## the files have the extension ".mapstat"
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##

mapWD <- "/path/to/RefSeq/phage/map_results"
taxonWD <- "path/to/RefSeq/phage/taxon_names"
setwd(mapWD)
mfiles <- dir()

# get the mapstat files
mfiles <- mfiles[which(sapply(strsplit(mfiles,"[.]"), "[[", 2, simplify=TRUE) == "mapstat")]

dat <- NULL
# i=1
for(i in 1:length(mfiles)){
  # get the sample names
  sample_name <- sapply(strsplit(mfiles[i], "_"), function(x) paste(x[1], "_", x[2], sep=""), simplify = TRUE)

  # read the mapstat files. Use  `quote = "" ` to ignore quotes
  temp <- read.table(file=mfiles[i], sep="\t", skip=6, header=T, quote = "")
  
  # extract the accession numbers
  accs <- sapply(strsplit(sapply(strsplit(sapply(strsplit(temp[,"refSequence"], " "), "[[", 1), "[|]"), "[[", 2), "_"), 
                 function(x) paste(x[1], "_", x[2], sep=""), simplify=TRUE)
  
  # remove the refSeq headers
  temp <- temp[, -which(colnames(temp) == c("refSequence"))]
  
  # merge the data
  temp <- data.frame(sample_name=sample_name, accs=accs, temp)
  
  # combine with other data
  dat <- rbind(temp, dat)
}
gc()

# get unique accession numbers
uniq_accs <- data.frame(accs = unique(dat$accs))

# save the file to taxonWD
write.table(uniq_accs, file = paste(taxonWD, "/uniq_accs.txt", sep=""), sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
