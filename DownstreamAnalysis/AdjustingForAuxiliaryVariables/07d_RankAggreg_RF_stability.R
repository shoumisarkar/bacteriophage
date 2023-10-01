library(readxl)
library(writexl)

setwd("path/to/rankAggreg_final_ordering")
dir()
rankAggreg_all_samples = read_xlsx("Rank_aggreg_RF_20.xlsx", col_names = T)

count_mat = matrix(numeric(20), nrow=1, ncol=20)
count_df = as.matrix(count_mat)
colnames(count_df) = rankAggreg_all_samples$accNo

for(i in 1:nrow(rankAggreg_all_samples))
{
  search_term = rankAggreg_all_samples$accNo[i] 
  
  count = 0
  
  for(leave_out_row in 1:366)
  {
    path = "/path/to/Results/leave_one_out/rankAggreg_final_ordering"   
    setwd(path)
    dir()
    filename = paste0("RankAggreg_RF_20_leave_out_row_", leave_out_row, ".xlsx")
    rankAggreg_excluding_i = read_xlsx(filename, col_names = T)
    
    
    if(search_term %in% rankAggreg_excluding_i$accNo){
      count=count+1}
    
    print(leave_out_row)
  }
  
  count_df[1,i] = count
  print(i)
}

prop_df = count_df/366

setwd("path/to/downstream_analysis/Results/")
write_xlsx(data.frame(count_df), "rankAggreg_agreement_list_size_20_counts.xlsx", col_names = T)
write_xlsx(data.frame(prop_df), "rankAggreg_agreement_list_size_20_proportions.xlsx", col_names = T)
