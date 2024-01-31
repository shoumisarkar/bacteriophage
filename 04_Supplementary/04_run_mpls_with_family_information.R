#rm(list=setdiff(ls(), "boot.Y.199.phage.classes"))

library(MASS)
library(readxl)
library(writexl)
library(dplyr)
library(ade4)


##################################################################################
############################# Load phage families ################################
phage_families <- read_xlsx("/path/to/Python_output/updated_NCBI_taxa_families.xlsx")
colnames(phage_families) = c("accession_number", "TaxID", "family_name")
phage_families$family_name[which(phage_families$family_name=="Unknown")] = rep("Unclassified", times = length(which(phage_families$family_name=="Unknown")))
#phage_families = na.omit(phage_families)
phage_families = phage_families[phage_families$family_name!="",]
phage_families$family_name = as.factor(phage_families$family_name)

##################################################################################
############################# Load phage dataframe ###############################

Phages <- as.data.frame(read_xlsx("/path/to/normalized_phage_OTU_counts.xlsx"))
sampleNames = Phages$samples
Phages$samples = NULL
Phages = apply(Phages, 2, as.numeric)
Phages = data.frame(Phages)
rownames(Phages) = sampleNames

# Initialize an empty list to store dataframes
phage_family_list <- list()
# Initialize an empty list to store two-column dataframes
few_column_dfs <- list()


for(family in levels(phage_families$family_name)) {
  # Convert family to character if it is a factor
  family_name <- as.character(family)
  
  accs <- phage_families$accession_number[phage_families$family_name == family_name] # Accession numbers corresponding to the current phage family
  col_ids <- which(colnames(Phages) %in% accs)
  
  if(length(col_ids > 0))
  {
    if(length(col_ids)==1)
    {
      subdat <- data.frame(samples = rownames(Phages), Phages[, col_ids])
      colnames(subdat)[2] = colnames(Phages)[col_ids]
    }else{
      subdat <- data.frame(samples = rownames(Phages), Phages[, col_ids])
    }
    
    # Check if the dataframe has only two columns (excluding 'samples')
    if(ncol(subdat) < 5) {
      # Add it to the few_column_dfs list
      few_column_dfs[[family_name]] <- subdat
      # Remove the dataframe from phage_family_list
      phage_family_list[[family_name]] <- NULL
    } else {
      # Assign the dataframe to the list with the family name as the key
      phage_family_list[[family_name]] <- subdat
    }
  }
}

# Merge all two-column dataframes
Miscellaneous <- Reduce(function(x, y) merge(x, y, by = "samples", all = TRUE), few_column_dfs)
rownames(Miscellaneous) = rownames(Phages)

# Remove NULL entries from phage_family_list
phage_family_list <- phage_family_list[sapply(phage_family_list, function(x) !is.null(x))]

#Add in the entries for `Miscellaneous` family of phages
phage_family_list[["Misc. phage families"]] <- Miscellaneous

##################################################################################

##########################################################
## Phages (preprocessed with phage family information): ##
##########################################################

# Merge all dataframes in the list by the 'samples' column
phages_ordered_by_families <- Reduce(function(x, y) merge(x, y, by = "samples", all = TRUE), phage_family_list)
rownames(phages_ordered_by_families) = rownames(Phages)

# Create a vector with the number of columns in each component dataframe in the list
phage_block_size <- sapply(phage_family_list, ncol) - 1 #subtracting a column because `samples` form the first column in each component dataframe
phage_block_size

##########################################################
#################### Resistomes:  ########################
##########################################################

Resistomes <- as.data.frame(read_xlsx("/path/to/normalized_resistome_abundances.xlsx"))
sampleNames = Resistomes$samples
Resistomes$samples = NULL
Resistomes = data.frame(apply(Resistomes, 2, as.numeric))
rownames(Resistomes) =  sampleNames

##########################################################
#################### Confounders: ########################
##########################################################

Climate <- read_excel("/path/to/Climate/variables.xlsx")
Demographics <- read_excel("/path/to/Demographics/variables.xlsx")
Landscape <- read_excel("/path/to/Landscape/variables.xlsx")

#Character to numeric:
Landscape$Coastal_City = as.numeric(as.factor(Landscape$Coastal_City))

Climate$City = NULL
Demographics$City = NULL
Landscape$City = NULL
#Landscape$Koppen_climate_classification_coldest_to_warmest = NULL

start = nchar("gCSD16_A")
end = nchar("gCSD16_AKL")

all_predictors = phages_ordered_by_families
all_predictors$cityCode = substr(rownames(Phages), start, end)

all_predictors = merge(all_predictors, Climate, by="cityCode")
all_predictors = merge(all_predictors, Demographics, by="cityCode")
all_predictors = merge(all_predictors, Landscape, by="cityCode")

#center and scale numeric variables other than those in the Phages block.
all_predictors[,-c(1:(ncol(phages_ordered_by_families)+2))] = 
  all_predictors[,-c(1:1:(ncol(phages_ordered_by_families)+2))]%>%mutate_if(is.numeric,scale)

all_predictors$cityCode = NULL
all_predictors$samples = NULL

rownames(all_predictors) <- sampleNames
rownames(Resistomes) <- sampleNames

#####################################################################
################### Apply MBPLS #####################################
#####################################################################

#Blocks
blo = unname(c(phage_block_size,
               Climate=(ncol(Climate)-1),
               Demographics = (ncol(Demographics)-1),
               Landscape = (ncol(Landscape)-1)
)) #`samples`/`cityCodes` column is omitted, hence the -1.

tab.names = c(names(phage_family_list), c("Climate", "Demographics", "Landscape"))

#store predictor blocks in ktab format
ktabX = ktab.data.frame(df = all_predictors,
                        blocks = blo,
                        tabnames = tab.names
)

#Make the response of dudi class
dudiY = dudi.pca(df = Resistomes, center = T, scale = T, scannf = F, nf = 5)
mod4.phage.families <- mbpls(dudiY, ktabX, scale = TRUE, option = "uniform", scannf = FALSE, nf = 5)

setwd("/path/to/save/mpls/model/object")
save(mod4.phage.families, file = "mod4.phage.families.Rdata")

# To bootstrap:
start_bs = Sys.time()
set.seed(123456)
boot.Y.199.phage.families <- randboot(mod4.phage.families, nrepet = 199, optdim = 5)
end_bs = Sys.time()
elapsed_bs = end_bs - start_bs; elapsed_bs

setwd("/path/to/save/bootstrapped/mpls/objects")
save(boot.Y.phage.families, file = "boot.Y.phage.families.Rdata")
