library(readxl)
library(tidyverse)

# The colors
BLUE <- "#076fa2"
RED <- "#E3120B"
BLACK <- "#202020"
GREY <- "grey50"
PINK <- "#E16fa8"

###################
#### MACROLIDE ####
###################

macrolide <- read_excel("path/to/downstream_analysis/Results/univariate_RF_importance_scores/importances_univariate_RF_macrolide.xlsx")[1:20,]

macrolide = macrolide[order(macrolide$Overall, decreasing = F),]
macrolide$specieName = factor(macrolide$specieName, levels = macrolide$specieName)

macrolide_plt <- ggplot(macrolide[1:20,]) +
  geom_col(aes(Overall, specieName), fill = BLUE, width = 0.6) + ggtitle("Macrolide") +
  labs(y = "Phages", x = "Mean decrease in accuracy") + theme_bw(); macrolide_plt

#####################
#### BETA.LACTAM ####
#####################

beta.lactam <- read_excel("path/to/downstream_analysis/Results/univariate_RF_importance_scores/importances_univariate_RF_beta.lactam.xlsx")[1:30,]
beta.lactam = (beta.lactam[-which(startsWith(beta.lactam$specieName, "Unknown")),])[1:20,]

beta.lactam = beta.lactam[order(beta.lactam$Overall, decreasing = F),]
beta.lactam$specieName = factor(beta.lactam$specieName, levels = beta.lactam$specieName)

beta.lactam_plt <- ggplot(beta.lactam) +
  geom_col(aes(Overall, specieName), fill = BLUE, width = 0.6) + ggtitle("Beta lactam") +
  labs(y = "Phages", x = "Mean decrease in accuracy") + theme_bw(); beta.lactam_plt

######################
#### TETRACYCLINE ####
######################

tetracycline <- read_excel("path/to/downstream_analysis/Results/univariate_RF_importance_scores/importances_univariate_RF_tetracycline.xlsx")[1:30,]
tetracycline = (tetracycline[-which(startsWith(tetracycline$specieName, "Unknown")),])[1:20,]

tetracycline = tetracycline[order(tetracycline$Overall, decreasing = F),]
tetracycline$specieName = factor(tetracycline$specieName, levels = tetracycline$specieName)

tetracycline_plt <- ggplot(tetracycline[1:20,]) +
  geom_col(aes(Overall, specieName), fill = BLUE, width = 0.6) + ggtitle("Tetracycline") +
  labs(y = "Phages", x = "Mean decrease in accuracy") + theme_bw(); tetracycline_plt

library(ggpubr)
RF_imp_three_classes = ggarrange(macrolide_plt, beta.lactam_plt, tetracycline_plt, nrow=3); RF_imp_three_classes

setwd("path/to/downstream_analysis/Results/Figures")
ggsave("RF_imp_three_classes.png", width = 10, height = 12)

