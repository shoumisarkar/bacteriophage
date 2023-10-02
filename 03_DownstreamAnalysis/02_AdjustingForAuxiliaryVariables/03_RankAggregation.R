library(stringr)
library(readxl)
library(writexl)
library(RankAggreg)

### Get names of the resistomes: ####
Resistomes <- as.data.frame(read_xlsx("/path/to/ResFinder/map_results/AMR_counts_norm.xlsx"))
Resistome_names = colnames(Resistomes)[-c(1)]

#########################################
###### Load RF importance scores:  ######
#########################################

setwd("path/to/Results/univariate_RF_importance_scores/")

RF_imp_scores = data.frame(matrix(ncol=17,nrow=1190, 
                                   dimnames=list(NULL, Resistome_names)))
RF_var = data.frame(matrix(ncol=17,nrow=1190, 
                                dimnames=list(NULL, Resistome_names)))

for(i in 1:17)
{
  filename_RF = paste0("importances_univariate_RF_", Resistome_names[i], ".xlsx")
  RF_var[, i] = unlist(read_xlsx(filename_RF, col_names = T)$var)
  RF_imp_scores[, i] = read_xlsx(filename_RF, col_names = T)$Overall
}

###############################
######Rank Aggregation: #######
###############################

x = (t(RF_var))
w = (t(RF_imp_scores))

(CES_20_RF <- RankAggreg(x, 20, w, "CE", "Spearman", rho=.1, convIn=5))


########################
#To add species names instead of NCBI accession numbers:

### Accessing mapping of accession numbers to specie names for phages:
accs_taxID_species <- read.delim("path/to/RefSeq/phage/taxon_names/accs_taxID_species.txt", stringsAsFactors=TRUE)  

### Functions to carry out this mapping:
getSpeciesName = function(x)
{
  ind = which(accs_taxID_species$accession_number %in% x) 
  impSpecies = accs_taxID_species$specie_name[ind]
  #needs stringr package
  impSpecies = str_extract(impSpecies, '\\b[^;]+$')
  return(impSpecies)
}

#for a vector of Specie names:
getSpeciesName_vector = function(accID_vector)
{ specieNames = c()
for(i in accID_vector)
{specieNames = c(specieNames, (getSpeciesName(i)))
}
return(specieNames)
}

df_top_20_phages_RF = data.frame(accNo=CES_20_RF$top.list, 
                                 specieName=getSpeciesName_vector(CES_20_RF$top.list)
)


setwd("path/to/Results/rankAggreg_final_ordering")
write_xlsx(df_top_20_phages_RF, "Rank_aggreg_RF_20.xlsx", col_names = T)