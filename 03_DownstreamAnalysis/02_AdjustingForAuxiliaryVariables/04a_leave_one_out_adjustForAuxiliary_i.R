#This script obtains the adjusted residuals for the resistomes, leaving out row i

library(writexl)
library(readxl)
library(dplyr)

#########################################################################
####################  Restructure the data ##############################
#########################################################################

#Resistomes:
Resistomes <- as.data.frame(read_xlsx("/path/to/ResFinder/map_results/AMR_counts_norm.xlsx"))
sampleNames = Resistomes$samples
Resistomes$samples = NULL
Resistomes = data.frame(apply(Resistomes, 2, as.numeric))
rownames(Resistomes) =  sampleNames

#Phages
Phages <- as.data.frame(read_xlsx("/path/to/RefSeq/map_results/OTU_counts_norm.xlsx"))
sampleNames = Phages$samples
Phages$samples = NULL
Phages = apply(Phages, 2, as.numeric)
Phages = data.frame(Phages)
rownames(Phages) = sampleNames

#Confounders:
Climate <- read_excel("path/to/data/Confounders.xlsx", sheet = "Climate")
Demographics <- read_excel("path/to/data/Confounders.xlsx", sheet = "Demographics")
Landscape <- read_excel("path/to/data/Confounders.xlsx", sheet = "Landscape")

#Character to factor:
Landscape$Coastal_City = as.factor(Landscape$Coastal_City)
#Factor to numeric:
Landscape$Coastal_City = as.numeric(Landscape$Coastal_City)

Climate$City = NULL
Demographics$City = NULL
Landscape$City = NULL

start = nchar("gCSD16_A")
end = nchar("gCSD16_AKL")

all_predictors = Phages
all_predictors$cityCode = substr(rownames(Phages), start, end)

all_predictors = merge(all_predictors, Climate, by="cityCode")
all_predictors = merge(all_predictors, Demographics, by="cityCode")
all_predictors = merge(all_predictors, Landscape, by="cityCode")

colnames(all_predictors)[1]

#center and scale numeric variables other than those in the Phages block.
all_predictors[,-c(1:1190)] = all_predictors[,-c(1:1190)]%>%mutate_if(is.numeric,scale)

scaled_Climate = all_predictors[,colnames(Climate)[-c(1)]]
scaled_Demographics = all_predictors[,colnames(Demographics)[-c(1)]]
scaled_Landscape = all_predictors[,colnames(Landscape)[-c(1)]]
#these, combined with Phages, are the predictors.
#Resistomes are the (multivariate) responses.

all_predictors$cityCode = NULL
rownames(all_predictors) <- sampleNames
rownames(Resistomes) <- sampleNames

confounders = as.data.frame(all_predictors[,-c(1:1190)])

########################################################################
##################### Leave-one-sample-out: ############################
########################################################################

for(leave_out_row in 1:366)
{
  Resistomes_leave_one_out = Resistomes[-c(leave_out_row),]
  confounders_leave_one_out = confounders[-c(leave_out_row),]
  
  setwd("path/to/save/Leave_one_out_residuals/")
  folder_name = paste0("leave_out_sample_", leave_out_row)
  dir.create(folder_name)
  folder_path = paste0("path/to/save/Leave_one_out_residuals/", folder_name)
  setwd(folder_path)
  
  #################################################################
  
  for(i in 1:17)
  { 
    #############################
    ### Get the lm residuals: ###
    #############################
    
    mod_i = lm(formula = Resistomes_leave_one_out[,i] ~ as.matrix(confounders_leave_one_out))
    
    temp_df = data.frame(residual = mod_i$residuals)
    
    filename = paste0("lm_residuals_", colnames(Resistomes_leave_one_out)[i], ".xlsx")
    write_xlsx(temp_df, filename, col_names = T)
    
  }
  
}
