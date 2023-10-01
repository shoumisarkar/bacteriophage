library(readxl)
library(writexl)
library(randomForest)
library(caret)

#########################################################################
####################  Load the data #####################################
#########################################################################

#################################
# Load the normalized phage data:
setwd("/path/to/Data/folder")
phage_df = data.frame(read_xlsx(path = "OTU_counts_norm.xlsx", col_names = T))
rownames(phage_df) = phage_df$samples
phage_df$samples = NULL 

#Convert character entries to numeric:
phage_mat = (apply(as.matrix(phage_df), 2,          
                   function(x) as.numeric(as.character(x)))) #this conversion results in a matrix.
dimnames(phage_mat)[[1]] = rownames(phage_df) #the rownames are lost, so assign them again
phage_df = data.frame(phage_mat) #convert it back to a dataframe
##################################

#Get antibiotic class names:
Resistomes <- as.data.frame(read_xlsx("/path/to/ResFinder/map_results/AMR_counts_norm.xlsx"))
antibiotic_class_names = colnames(Resistomes)[-c(1)]

#############################################
### Functions: ##############################
#############################################
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

#############################################

#Load the responses one by one, and fit & save the univariate RF models.
for(class_name in antibiotic_class_names)
{
  residuals_filename = paste0("lm_residuals_", class_name, ".xlsx")
  setwd("path/to/save/Residuals/")
  residuals_for_current_antibiotic = read_xlsx(residuals_filename, col_names = T)
  
  rownames(residuals_for_current_antibiotic) = Resistomes$samples #keep rownames consistent with phage_df, so that we can cbind them:
  temp_df = cbind(residuals_for_current_antibiotic, phage_df) 
  
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
  #i_scores$fullTaxonomy = getFullTaxonomy(i_scores$var)
  
  
  #save the model: 
  model_name = paste0("univariate_RF_", class_name, ".Rdata")
  setwd("path/to/Saved objects/Univariate_RF_models/")
  save(rf.fit, file = model_name)
  
  #save the importance scores:
  setwd("path/to/Results/univariate_RF_importance_scores/")
  file_path = paste0("importances_univariate_RF_", class_name, ".xlsx")
  write_xlsx(i_scores, path = file_path, col_names = T)
}
