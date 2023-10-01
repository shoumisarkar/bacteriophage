library(readxl)
library(ade4)
library(dplyr)
library(adegraphics)

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

#####################################################################
################### Apply MBPLS #####################################
#####################################################################

#Blocks
blo = c(ncol(Phages),
        (ncol(Climate)-1),
        (ncol(Demographics)-1),
        (ncol(Landscape)-1)
) #cityCodes column is omitted, hence the -1.

tab.names = c("Phages", "Climate", "Demographics", "Landscape")

#store predictor blocks in ktab format
ktabX = ktab.data.frame(df = all_predictors,
                        blocks = blo,
                        tabnames = tab.names
)

dudiY = dudi.pca(df = Resistomes, center = T, scale = T, scannf = F, nf = 5)


mod4 <- mbpls(dudiY, ktabX, scale = TRUE, option = "uniform", scannf = FALSE, nf = 4)
summary(mod4)

# To bootstrap:
start_bs = Sys.time()
set.seed(123456)
boot.Y.199 <- randboot(mod4, nrepet = 199, optdim = 4)
end_bs = Sys.time()
elapsed_bs = end_bs - start_bs; elapsed_bs

boot.Y.199$vipc
boot.Y.199$bipc

plot5 <- plot(boot.Y.199$vipc, main = "vipc", plot = FALSE)
plot6 <- plot(boot.Y.199$bipc, main = "bipc", plot = FALSE)
ADEgS(list(plot5, plot6))


setwd("path/to/downstream_analysis/folder")

bipc_bootstrap = data.frame(bipc = boot.Y.199$bipc$obs, 
                            CI_lb = boot.Y.199$bipc$stats[,1],
                            CI_ub = boot.Y.199$bipc$stats[,2])
bipc_bootstrap$Block = rownames(bipc_bootstrap)
bipc_bootstrap = bipc_bootstrap[order(bipc_bootstrap$bipc, decreasing = T),c(4,1:3)]

setwd("path/to/downstream_analysis/Results/folder")
write_xlsx(bipc_bootstrap, "bipc_bootstrapped_nrepet_199.xlsx", col_names = T)

vipc_bootstrap = data.frame(vipc = boot.Y.199$vipc$obs, 
                            CI_lb = boot.Y.199$vipc$stats[,1],
                            CI_ub = boot.Y.199$vipc$stats[,2])
vipc_bootstrap$Variable = rownames(vipc_bootstrap)
vipc_bootstrap = vipc_bootstrap[order(vipc_bootstrap$vipc, decreasing = T),c(4,1:3)]

setwd("path/to/downstream_analysis/Results/folder")
write_xlsx(vipc_bootstrap, "vipc_bootstrapped_nrepet_199.xlsx", col_names = T)
