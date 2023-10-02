#!/bin/sh
#SBATCH --job-name=ungzip		# Job name
#SBATCH --account=somnath.datta			# your own sponser or account
#SBATCH --qos=somnath.datta-b                # your own sponser or account
#SBATCH --mail-type=ALL               	# Mail events
#SBATCH --mail-user=shoumisarkar@ufl.edu      # Where to send email	
#SBATCH --ntasks=1                    	# Run on a single machine (node)
#SBATCH --cpus-per-task 20		  		# Run on a several CPU
#SBATCH --mem-per-cpu=5gb	            # Memory limit
#SBATCH --time=24:00:00               	# Time: hrs:min:sec
#SBATCH --output=ungzip_missing_%j.out   		# Output and error log 

pwd; hostname; date 

max_cpus=20

module load dsrc/2.0 
module load gcc/5.2.0
module load kneaddata/0.7.4
module load bowtie2/2.3.5.1
module load trimmomatic/0.32
module load bbmap/38.90

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
   #Commenting this as gzip downloads are already done.
   #echo "Downloading ${file} files..."
   #lftp -u camda,camdac4md4 sftp://metasub.ethz.ch -e "cd camda2023; get ${file1} -o #/orange/somnath.datta/CAMDA2023/CAMDA_Phages_Forensic/fastq_files_1; get ${file2} -o #/orange/somnath.datta/CAMDA2023/CAMDA_Phages_Forensic/fastq_files_1; get ${file}_2.fastq.gz -o #/path/to/store/fastq_files; quit" 
 
   ####################################################
   ##  Begins the loop for processing with KneadData ##	
   ####################################################
  
   cd /path/to/fastq_files

   # If files are gzipped, you unzip the files and replace the original gzipped file
   gzip -d -f ${file}_1.fastq.gz > ${file}_1.fastq
   gzip -d -f ${file}_2.fastq.gz > ${file}_2.fastq

done

date
