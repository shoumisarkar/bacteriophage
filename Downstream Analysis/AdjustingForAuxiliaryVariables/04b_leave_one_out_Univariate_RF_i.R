#This script obtains the RF models, leaving out row i

#rm(list = ls())
library(readxl)
library(writexl)
library(randomForest)
library(caret)
library(stringr)

leave_out_row = commandArgs(trailingOnly=TRUE)
print(paste0("Leaving out sample # = ", leave_out_row))

str(leave_out_row)
leave_out_row = as.numeric(leave_out_row)
print(leave_out_row)

#########################################################################
####################  Load the data #####################################
#########################################################################

# Load the normalized phage data:
setwd("/path/to/Data/folder")
phage_df = data.frame(read_xlsx(path = "/path/to/RefSeq/map_results/OTU_counts_norm.xlsx", col_names = T))
rownames(phage_df) = phage_df$samples
phage_df$samples = NULL 

#Convert character entries to numeric:
phage_mat = (apply(as.matrix(phage_df), 2,          
                   function(x) as.numeric(as.character(x)))) #this conversion results in a matrix.
dimnames(phage_mat)[[1]] = rownames(phage_df) #the rownames are lost, so assign them again
phage_df = data.frame(phage_mat) #convert it back to a dataframe
##################################

#Get antibiotic class names:
Resistomes <- as.data.frame(read_xlsx("/orange/somnath.datta/CAMDA2023/CAMDA_Phages_Forensic/downstream_analysis_v2/Data/AMR_counts_norm_ordered.xlsx"))
antibiotic_class_names = colnames(Resistomes)[-c(1)]

#############################################
### Functions: ##############################
#############################################
### Accessing mapping of accession numbers to specie names for phages:
accs_taxID_species <- read.delim("/orange/somnath.datta/CAMDA2023/CAMDA_Phages_Forensic/RefSeq/phage/taxon_names_366_samples/accs_taxID_species.txt", stringsAsFactors=TRUE)  

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


#############################################


Resistomes_leave_one_out = Resistomes[-c(leave_out_row),]
phage_df_leave_one_out = phage_df[-c(leave_out_row),]


#Load the responses one by one, and fit & save the univariate RF models.
for(class_name in antibiotic_class_names)
{
  residuals_filename = paste0("lm_residuals_", class_name, ".xlsx")
  residuals_path = paste0("/path/to/Leave_one_out_residuals/leave_out_sample_", leave_out_row)
  setwd(residuals_path)
  
  residuals_for_current_antibiotic = read_xlsx(residuals_filename, col_names = T)
  
  rownames(residuals_for_current_antibiotic) = Resistomes_leave_one_out$samples #keep rownames consistent with phage_df_leave_one_out, so that we can cbind them:
  temp_df = cbind(residuals_for_current_antibiotic, phage_df_leave_one_out) 
  
  # Set a random seed
  set.seed(100)
  bestmtry <- data.frame(tuneRF(x = temp_df[,-c(1)], y = temp_df$residual,stepFactor = 1.2,
                                improve = 0.01, trace=T, plot= T, ntreeTry = 500)) 
  bestmtry = bestmtry$mtry[which((bestmtry$OOBError) == min(bestmtry$OOBError))]
  
  rf.fit = randomForest(residual ~ ., data=temp_df, mtry=bestmtry,
                        ntree = 500, importance=TRUE); rf.fit
  
  #Conditional=True, adjusts for correlations between predictors.
  i_scores <- varImp(rf.fit, conditional=TRUE)
  i_scores <- i_scores %>% tibble::rownames_to_column("var") 
  i_scores = i_scores[order(i_scores$Overall, decreasing = T),] #ordered from highest to lowest
  
  i_scores$specieName = getSpeciesName_vector(i_scores$var)

  #save the importance scores:
  setwd("/path/to/Results/leave_one_out/univariate_RF_importance_scores/")
  dir_name = paste0("leave_out_sample_", leave_out_row)
  dir.create(dir_name)
  file_path = paste0(dir_name,"/importances_univariate_RF_", class_name, ".xlsx")
  write_xlsx(i_scores, path = file_path, col_names = T)
}

print(leave_out_row)  

