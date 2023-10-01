library(writexl)
library(readxl)
library(dplyr)
library(MASS)

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
Phages <- as.data.frame(read_xlsx("/path/to/RefSeq/map_results/phage_counts_norm.xlsx"))
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

#################################################################

for(i in 1:17)
{ 
  mod_i = lm(formula = Resistomes[,i] ~ as.matrix(confounders))
  
  temp_df = data.frame(residual = mod_i$residuals)
  
  setwd("path/to/save/Residuals/")
  filename = paste0("lm_residuals_", colnames(Resistomes)[i], ".xlsx")
  write_xlsx(temp_df, filename, col_names = T)
}