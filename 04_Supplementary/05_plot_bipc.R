library(ade4)
library(adegraphics)
library(ggplot2)
library(tidyverse)

imagePath <- "path/to/MPLS_bootstrapped_object"
dir()

load("path/to/boot.Y.phage.families.Rdata")
summary(boot.Y.phage.families)
K <- 21
P <- 451

#################################################################################
###############    to obtain and save the bipc interval plot:   ###################
#################################################################################

blk_imp = data.frame(boot.Y.phage.families$bipc$boot)

## block importance
blk_imp <- as.data.frame(boot.Y.phage.families$bipc$stats)
blk_imp <- cbind(blk_imp, boot.Y.phage.families$bipc$obs)
colnames(blk_imp) <- c("LCL", "UCL", "center")
blk_imp <- data.frame(id=rownames(blk_imp), blk_imp)
rownames(blk_imp) <- NULL

env_id = which(blk_imp$id %in% c("Climate", "Demographics", "Landscape"))
blk_imp = blk_imp[c(env_id, setdiff(1:nrow(blk_imp), env_id)),]

unclassif_id = which(blk_imp$id %in% c("Unclassified"))
blk_imp = blk_imp[c(setdiff(1:nrow(blk_imp), unclassif_id), unclassif_id),]

blk_imp$id <- factor(blk_imp$id, levels = blk_imp$id)

ggplot(blk_imp, aes(x = id, y = center)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = UCL, ymin = LCL)) +
  geom_hline(yintercept=1/K, linetype="dashed", size=0.9, color = "grey") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank() ,
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), #change legend text font size
        axis.text=element_text(size=14), # change axis text font size
        axis.title=element_text(size=18,face="bold"), # change axis title size
        strip.text.x = element_text(size = 14,face="bold"), #change facet_grid text size
        plot.title = element_text(hjust = 0.5, size=26,face="bold")) + # change size of plot title
  labs(x= "",
       y= "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(imagePath,"/phage_families_blk_imp.pdf", sep=""), width = 9.91, height = 6.96, dpi = 300,units="in")

#################################################################################

# variable importance
var_imp <- as.data.frame(boot.Y.phage.families$vipc$stats)
var_imp <- cbind(var_imp, boot.Y.phage.families$vipc$obs)
colnames(var_imp) <- c("LCL", "UCL", "center")
var_imp <- data.frame(id=rownames(var_imp), var_imp)
rownames(var_imp) <- NULL

var_imp <- var_imp[order(abs(var_imp$center), decreasing = T), ]
VI <- var_imp[sapply(strsplit(var_imp$id, "_"), "[[", 1, simplify=TRUE) != "NC", ]
# VI <- var_imp[1:50, ]
VI$id <- factor(VI$id, levels=VI$id)
ggplot(VI, aes(x = id, y = center)) +
  geom_point(size =2) +
  geom_errorbar(aes(ymax = UCL, ymin = LCL)) +
  geom_hline(yintercept=1/P, linetype="dashed", size=0.9, color = "grey") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank() ,
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), #change legend text font size
        axis.text.y=element_text(size=14), # change axis text font size
        axis.text.x=element_text(size=12), # change axis text font size
        axis.title=element_text(size=18,face="bold"), # change axis title size
        strip.text.x = element_text(size = 14,face="bold"), #change facet_grid text size
        plot.title = element_text(hjust = 0.5, size=26,face="bold")) + # change size of plot title
  labs(x= "", y= "",) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(imagePath,"/phage_families_var_imp.pdf", sep=""), width = 9.91, height = 6.96, dpi = 300,units="in")

# variable importance
var_imp2 <- var_imp[sapply(strsplit(var_imp$id, "_"), "[[", 1, simplify=TRUE)=="NC", ]
var_imp2 <- var_imp2[order(abs(var_imp2$center), decreasing=TRUE), ]
var_imp2 <- var_imp2[1:50, ]
var_imp2$id <- factor(var_imp2$id, levels=var_imp2$id)

ggplot(var_imp2, aes(x = id, y = center)) +
  geom_point(size =2) +
  geom_errorbar(aes(ymax = UCL, ymin = LCL)) +
  geom_hline(yintercept=1/P, linetype="dashed", size=0.9, color = "grey") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank() ,
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12), #change legend text font size
        axis.text.y=element_text(size=14), # change axis text font size
        axis.text.x=element_text(size=12), # change axis text font size
        axis.title=element_text(size=18,face="bold"), # change axis title size
        strip.text.x = element_text(size = 14,face="bold"), #change facet_grid text size
        plot.title = element_text(hjust = 0.5, size=26,face="bold")) + # change size of plot title
  labs(x= "", y= "",) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(imagePath,"/phage_families_var_imp2.pdf", sep=""), width = 9.91, height = 6.96, dpi = 300,units="in")

## variable importance for some resistomes
resistomes <- names(boot.Y.phage.families$XYcoef)

for(k in 1:length(resistomes)){
  t1 <- resistomes[k]
  t2 <- boot.Y.phage.families$XYcoef[[t1]]
  
  # variable importance
  vimp <- as.data.frame(t2$stats)
  vimp <- cbind(vimp, t2$obs)
  colnames(vimp) <- c("LCL", "UCL", "center")
  vimp <- data.frame(id=rownames(vimp), vimp)
  rownames(vimp) <- NULL
  
  # vimp <- vimp[sapply(strsplit(vimp$id, "_"), "[[", 1, simplify=TRUE)=="NC", ]
  vimp <- vimp[order(abs(vimp$center), decreasing=TRUE), ]
  vimp <- vimp[1:20, ]
  vimp$id <- factor(vimp$id, levels=vimp$id)
  
  ggplot(vimp, aes(x = id, y = center)) +
    geom_point(size =2) +
    geom_errorbar(aes(ymax = UCL, ymin = LCL)) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank() ,
          legend.key.size = unit(0.8, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=12), #change legend title font size
          legend.text = element_text(size=12), #change legend text font size
          axis.text.y=element_text(size=14), # change axis text font size
          axis.text.x=element_text(size=12), # change axis text font size
          axis.title=element_text(size=18,face="bold"), # change axis title size
          strip.text.x = element_text(size = 14,face="bold"), #change facet_grid text size
          plot.title = element_text(hjust = 0.5, size=26,face="bold")) + # change size of plot title
    labs(x= "", y= "",) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste0(imagePath,"/","phage_families_",t1,"_vimp",".pdf", sep=""), width = 9.91, height = 6.96, dpi = 300,units="in")
}
