library(readxl)
library(ade4)
library(xtable)

setwd("path/to/MPLS_bootstrap_object_for_phage_families/")

## Use of the original (centered and scaled) datasets X and Y
load("MPLS_model_object_for_phage_families.Rdata") #this loads object `mod` (say). Then,
X = mod$tabX ; Y = mod$tabY

## Load bootstrapped XYcoef for each of the responses
load("boot.Y.phage.families.Rdata")

XYcoef = data.frame(aminoglycoside = boot.Y.phage.families$XYcoef$aminoglycoside$obs, 
                    beta.lactam = boot.Y.phage.families$XYcoef$beta.lactam$obs,
                    colistin = boot.Y.phage.families$XYcoef$colistin$obs,
                    fosfomycin = boot.Y.phage.families$XYcoef$fosfomycin$obs,
                    fusidicacid = boot.Y.phage.families$XYcoef$fusidicacid$obs,
                    glycopeptide = boot.Y.phage.families$XYcoef$glycopeptide$obs,
                    macrolide = boot.Y.phage.families$XYcoef$macrolide$obs,
                    misc = boot.Y.phage.families$XYcoef$misc$obs,
                    nitroimidazole = boot.Y.phage.families$XYcoef$nitroimidazole$obs,
                    oxazolidinone = boot.Y.phage.families$XYcoef$oxazolidinone$obs,
                    phenicol = boot.Y.phage.families$XYcoef$phenicol$obs,
                    pseudomonicacid = boot.Y.phage.families$XYcoef$pseudomonicacid$obs,
                    quinolone = boot.Y.phage.families$XYcoef$quinolone$obs,
                    rifampicin = boot.Y.phage.families$XYcoef$rifampicin$obs,
                    sulphonamide = boot.Y.phage.families$XYcoef$sulphonamide$obs,
                    tetracycline = boot.Y.phage.families$XYcoef$tetracycline$obs,
                    trimethoprim = boot.Y.phage.families$XYcoef$trimethoprim$obs)

Yhat = data.frame(matrix(ncol = ncol(Y), nrow = nrow(Y)))
rownames(Yhat) = rownames(Y)

for(i in 1:ncol(Y))
{
  Yhat[,i] = as.matrix(X) %*% (as.matrix(XYcoef[,i]))
}
colnames(Yhat) = colnames(Y)

###################################
######## To compute R^2:  #########
###################################

ybar = mean(as.matrix(Y))
Y_minus_ybar_squared <- (Y - ybar)^2
Yhat_minus_ybar_squared  <- (Yhat - ybar)^2
Y_minus_Yhat_squared <- (Y - Yhat)^2

SS_tot = sum(Y_minus_ybar_squared)
SS_exp = sum(Yhat_minus_ybar_squared)
SS_res = sum(Y_minus_Yhat_squared)
ESS = SS_tot - SS_exp
R_squared_overall =  data.frame(R_sq = sort(SS_exp/SS_tot * 100, decreasing = T)); R_squared_overall
