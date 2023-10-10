###############################################################
## Truncate/format file names, make them easier to read 
###############################################################

filenames <- read_table2("path/to/all/all_filenames.txt", 
                              col_names = FALSE) #all_filenames.txt contain the names of all .fastq.gz files obtained from the CAMDA server

#Create a vector of all names:
names_vec = c()
for(i in ncol(filenames))
{
  names_vec = c(names_vec, filenames[,i])
}

names_vec = na.omit(c(filenames$X1,
                      filenames$X2,
                      filenames$X3))

startpos = numeric(length(names_vec)) + nchar("camda2023/") + 1
endpos = nchar(names_vec) - nchar('_1.fastq.gz')

truncated_names = unique(substr(names_vec, startpos, endpos))

write.table(truncated_names, 
            file = "/path/to/store/truncated_filenames.txt",
            sep = "\t",
            row.names = F, col.names = F, quote = F)
