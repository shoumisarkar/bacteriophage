library(tidyverse)

###### PHAGES SUMS BY SAMPLE #########

phage_df = data.frame(read_xlsx("/path/to/RefSeq/map_results/OTU_counts_norm.xlsx", col_names = T))
rownames(phage_df) = phage_df$samples
phage_df$samples = NULL

#convert character entries to numeric, converting phage_df to a matrix
phage_mat = (apply(as.matrix(phage_df), 2,          
                   function(x) as.numeric(as.character(x))))
rownames(phage_mat) = rownames(phage_df)

phage_sums = as.numeric(apply(phage_mat, 1, sum))
names(phage_sums) = rownames(phage_df)


###### AMR #########

amr_df = data.frame(read_xlsx("/path/to/ResFinder/map_results/AMR_counts_norm.xlsx", col_names = T))
rownames(amr_df) = amr_df$samples
amr_df$samples = NULL

#convert character entries to numeric, and keep amr_df as a matrix
amr_mat = (apply(as.matrix(amr_df), 2,          
                 function(x) as.numeric(as.character(x))))
rownames(amr_mat) = rownames(amr_df)

amr_df = data.frame(amr_df)

combined_wide = cbind(amr_df)
colnames(combined_wide) = c("Aminoglycoside", "Beta lactam", "Colistin",       
                            "Fosfomycin", "Fusidicacid", "Glycopeptide", "Macrolide", 
                            "Miscellanous", "Nitroimidazole", "Oxazolidinone", "Phenicol",       
                            "Pseudomonicacid", "Quinolone", "Rifampicin", "Sulphonamide",   
                            "Tetracycline", "Trimethoprim")

#Import mapping between cities and city codes
cityCodes <- read_excel("cityCodes.xlsx")

combined_wide$CityCode = as.factor(substr(rownames(combined_wide), 8,10))
combined_wide$Location = substr(rownames(combined_wide), 8,10)
combined_wide$Region = substr(rownames(combined_wide), 8,10)

for(i in 1:nrow(cityCodes))
{
  for(j in 1:nrow(combined_wide))
  {
    if(substr(rownames(combined_wide)[j], 8,10) == cityCodes$code[i])
    {
      combined_wide$Location[j] = cityCodes$city[i]
      combined_wide$Region[j] = cityCodes$region[i]
    }
  }
}

combined_wide$Location = as.factor(combined_wide$Location)
combined_wide$Region = as.factor(combined_wide$Region)

combined_long = gather(data = combined_wide, key = Resistomes, value=Abundances, Aminoglycoside, 
                       `Beta lactam`, Colistin, Fosfomycin, Fusidicacid,
                       Glycopeptide,  Macrolide, Miscellanous, Nitroimidazole,
                       Oxazolidinone, Phenicol, Pseudomonicacid, Quinolone, 
                       Rifampicin, Sulphonamide, Tetracycline, Trimethoprim)


combined_long$Abundances = as.numeric(combined_long$Abundances)

# Group by sum of multiple columns using across()
combined_long_2 <- combined_long %>% group_by(Region, Location, Resistomes) %>% 
  summarise(across(c(Abundances),sum),
            .groups = 'drop') %>%
  as.data.frame() ; combined_long_2

variable_names <- list(
  "Region" = "Region" ,
  "Location" = "Location",
  "Resistomes" = "Resistomes",
  "Abundances" = "Abundances"
)

region_names = levels(combined_long_2$Region)

variable_labeller2 <- function(variable,value){
  if (variable=='measure') {
    return(variable_names[value])
  } else {
    return(region_names)
  }
}

colorscheme1 = c('lightgoldenrod4', 'brown2', 'mediumblue', 'skyblue1',
                 'chocolate', 'springgreen2', 'violetred1', 'purple4', 'lightsteelblue1',
                 'turquoise4', 'grey20', 'lightgoldenrod1', 'hotpink3', 'royalblue2',
                 'darkolivegreen3', 'plum3', 'darkgreen')

colorscheme2 = c('goldenrod2', 'darkgreen', 'orangered1', 'mediumblue', 'skyblue1',
                 'chocolate', 'springgreen2', 'violetred1', 'purple4', 'lightsteelblue1',
                 'turquoise4', 'grey20', 'lightgoldenrod1', 'hotpink3', 'royalblue2',
                 'darkolivegreen3', 'plum3')

colorscheme3 = c('goldenrod2', 'darkgreen', 'orangered1', 'mediumblue', 'skyblue1',
                 'pink', 'gray58', 'chocolate', 'purple4', 'lightsteelblue1',
                 'turquoise4', 'grey20', 'darkolivegreen3', 'hotpink3', 'royalblue2',
                 'lightgoldenrod1', 'plum3')


## Plot:

p0 = ggplot(combined_long_2, aes(x = Location, y = Abundances, fill = Resistomes, label = Abundances)) +
  geom_bar(position="fill", stat = "identity", width=0.7) +
  labs( #title="Relative abundance of antibiotic classes and phages for each sample",
    x ="Locations", y = "Relative Abundance") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5,
                               hjust=1, size=15),
    legend.direction = 'horizontal', legend.position = 'bottom',
    strip.text.x = element_text(size = 12)) +
  facet_grid(.~Region, scales="free", space="free_x",  labeller= variable_labeller2) +
  scale_fill_manual(values=colorscheme3) ; p0


p = ggplot(combined_long_2, aes(x = Location, y = Abundances, fill = Resistomes, label = Abundances)) +
  geom_bar(position="fill", stat = "identity", width=0.7) +
  labs( #title="Relative abundance of antibiotic classes and phages for each sample",
    x ="Locations", y = "Relative Abundance") +
  facet_grid(.~Region, scales="free", space="free_x",  labeller= variable_labeller2) +
  scale_fill_manual(values=colorscheme3) + 
  theme_bw(base_size = 15); p

###############################################
###############################################
###### Relative abundance of resistomes #######
###############################################
###############################################

resistome <- read_excel("/path/to/ResFinder/map_results/AMR_counts_norm.xlsx")
resistome[,-c(1)] = apply(resistome[,-c(1)], 2, as.numeric)

resistome_abundance = data.frame(abundance = colSums(resistome[,-c(1)])[order(colSums(resistome[,-c(1)]), decreasing = T)])
resistome_abundance = resistome_abundance/sum(resistome_abundance) #relative
resistome_abundance$class = factor(rownames(resistome_abundance), levels = rownames(resistome_abundance))
resistome_abundance = resistome_abundance[order(resistome_abundance$class, decreasing = T),]
resistome_abundance$id = 1:17


resistome_plt <- ggplot(data=resistome_abundance) +
  aes(x=class,y=abundance, fill=class) +
  geom_bar(stat="identity", position=position_dodge(), fill = PINK, width = 0.6) +
  labs(x = "Antibiotic classes", y = "Relative abundance") + theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.2), angle = 90, vjust = 0.5, hjust=1)); resistome_plt

setwd("path/to/downstream_analysis/Results/Figures")
ggsave("relative_abundance_resistomes.png", width = 10, height = 6)
