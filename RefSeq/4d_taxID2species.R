## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
## Merge the accession number, NCBI taxon ids & specie names
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
rm(list = ls())

taxonWD <- "/path/to/RefSeq/phage/taxon_names"
setwd(taxonWD)

# read the file containing accession numbers and corresponding NCBI taxon IDs

### The line commented below contains the erroneous file with missing strings at line 3150.
#ncbi_taxa <- read.table("acc_NCBI_taxa.txt", header=F, sep="\t")
ncbi_taxa <- read.table("acc_NCBI_taxa_corrected.txt", header=F, sep="\t") ##This file is corrected.
colnames(ncbi_taxa) <- c("TaxID", "accession_number")

# read the file containing NCBI taxon IDs and corresponding specie names
species <- read.table("NCBI_taxa_species.txt", header=F, sep="\t")
colnames(species) <- c("TaxID", "specie_name")

# label unknown viruses
temp <- which(species$specie_name == "")
species[temp, "specie_name"] <- paste("Unknown_virus", 1:length(temp), sep="_")

# merge the accession number, NCBI taxon ids & specie names
if(length(ncbi_taxa$TaxID) == length(species$TaxID)){
  if(all(ncbi_taxa$TaxID == species$TaxID)){
    dat <- data.frame(cbind(ncbi_taxa, species))
    dat <- dat[, -which(colnames(dat) == "TaxID.1")]
  }
} else {
  
  # Note, multiple accession numbers may have the same taxon id.
  # Also, if an accession number may not correspond to any taxon id
  # This will be the case if an old NCBI Taxonomy database is used in the previous step
  dat <- NULL
  for(i in 1:nrow(ncbi_taxa)){
    tp1 <- ncbi_taxa[i, ]
    tp2 <- species[which(species$TaxID %in% tp1$TaxID), ]
    if(nrow(tp2) > 1){
      tp2 <- dplyr::distinct(tp2,TaxID,.keep_all=T)
    } else if(nrow(tp2) == 0){
      tp2 <- rbind(tp2, data.frame(TaxID=NA, specie_name="Unknown_virus"))
    }
    dat <- rbind(cbind(tp1, specie_name=tp2$specie_name), dat)
  }
  dat <- dat[order(as.numeric(rownames(dat))), ]
}

# save the file to taxonWD
write.table(dat, file = paste(taxonWD, "/accs_taxID_species.txt", sep=""), sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)


