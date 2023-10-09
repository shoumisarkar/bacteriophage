#!/bin/sh
#SBATCH --job-name=ungzip		# Job name
#SBATCH --account=group_name			# your own sponser or account
#SBATCH --qos=group_name                # your own sponser or account
#SBATCH --mail-type=ALL               	# Mail events
#SBATCH --mail-user=name@email.com      # Where to send email	
#SBATCH --ntasks=1                    	# Run on a single machine (node)
#SBATCH --cpus-per-task 1		  		# Run on a several CPU
#SBATCH --mem-per-cpu=1gb             # Memory limit
#SBATCH --time=01:00:00               	# Time: hrs:min:sec
#SBATCH --output=output_%j.out  	 	# Output and error log 
pwd; hostname; date 

dos2unix /path/to/truncated_filenames.txt

for file in $(cat /path/to/truncated_filenames.txt)
do 

   ####################################################	
   ####  download the file from the CAMDA server:  #### 
   ####################################################
  
   sep1="_1.fastq.gz"
   sep2="_2.fastq.gz"
   file1="${file}${sep1}"
   file2="${file}${sep2}"
    
   ####################################################
   ##  Begins the loop for processing with KneadData ##	
   ####################################################
  
   cd /path/to/fastq_files

   # If files are gzipped, you unzip the files and replace the original gzipped file
   gzip -d -f ${file}_1.fastq.gz > ${file}_1.fastq
   gzip -d -f ${file}_2.fastq.gz > ${file}_2.fastq

done

date
