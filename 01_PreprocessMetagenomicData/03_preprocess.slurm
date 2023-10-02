#!/bin/sh
#SBATCH --job-name=pre1		# Job name
#SBATCH --job-name=phage				# Job name
#SBATCH --account=group_name			# your own sponser or account
#SBATCH --qos=group_name                # your own sponser or account
#SBATCH --mail-type=ALL               	# Mail events
#SBATCH --mail-user=name@email.com      # Where to send email	
#SBATCH --ntasks=1                    	# Run on a single machine (node)
#SBATCH --cpus-per-task 30		  		# Run on a several CPU
#SBATCH --mem-per-cpu=500mb             # Memory limit
#SBATCH --time=01:00:00               	# Time: hrs:min:sec
#SBATCH --output=output_%j.out  	 	# Output and error log 

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

   echo "Downloading ${file} files..."
   lftp -u camda,camdac4md4 sftp://metasub.ethz.ch -e "cd camda2023; get ${file1} -o /path/to/store/fastq_files; get ${file2} -o /path/to/store/fastq_files; get ${file}_2.fastq.gz -o /path/to/store/fastq_files; quit" 
 
   ####################################################
   ##  Begins the loop for processing with KneadData ##	
   ####################################################
  
   cd /path/to/store/fastq_files

   # If files are gzipped, you unzip the files and replace the original gzipped file
   gzip -d -f ${file}_1.fastq.gz > ${file}_1.fastq
   gzip -d -f ${file}_2.fastq.gz > ${file}_2.fastq

   for f in $(ls *.fastq | sed -e 's/_1.fastq//' -e 's/_2.fastq//' | sort -u)
   do
  
   # removal of contaminated reads and trimming:

   echo "Preprocessing ${f}..."

   kneaddata -i ${f}_1.fastq -i ${f}_2.fastq -o /path/to/store/QC_results \
   -db /blue/somnath.datta/share/Homo_Sapiens --trimmomatic /path/to/Trimmomatic-0.36 \
   -t ${max_cpus} --trimmomatic-options "ILLUMINACLIP:/path/to/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10:8:keepBothReads SLIDINGWINDOW:4:30 MINLEN:60" \
   --bowtie2-options "--very-sensitive --dovetail"  

   echo "Finished preprocessing ${f}"
   done
   # End of the loop for the KneadData QCing

   echo "KneadData QCing complete"

   #Remove files that are no longer needed:
   echo "Deleting original .fastq files..."
   cd /path/to/fastq_files
   rm *
   echo ".fastq files removed"

   echo "Removing unneccessary kneaddata results..."
   cd /path/to/QC_results/
   rm *_kneaddata.trimmed.*
   rm *bowtie2*
   rm *unmatched*
   rm *kneaddata.log		    
   echo "Unnecessary kneaddata results removed"

   ##############################################
   ##### Repair disordered paired-end files #####
   ##############################################

   echo "Re-pair files using bbmap..."

   cd /path/to/QC_results

   for f in $(ls *.fastq | sed -e 's/_1.fastq//' -e 's/_2.fastq//' | sort -u)
   do
	repair.sh in1=${f}_1.fastq in2=${f}_2.fastq \
	out1=/path/to/store/bbmap_repaired_files/${f}_ORDERED_1.fastq \
	out2=/path/to/store/bbmap_repaired_files/${city}/${f}_ORDERED_2.fastq \
	outs=/path/to/store/bbmap_repaired_files/${city}/${f}_ORDERED_singleton.fastq repair
        echo "Repaired ${f}"
   done

   cd /path/to/bbmap_repaired_files
   rm *ORDERED_singleton.fastq
   echo "Removed *ORDERED_singleton.fastq files"

   ## Add comment to remove the *kneaddata_paired* files which are now unneeded.
   cd /path/to/QC_results
   rm *
   echo "KneadData outputs deleted."

done

### Keep the following lines if you would like the script to delete the unnecessary folders:
#echo "Deleting folders from the intermediate steps…"
# cd /path/to/parent/directory/
#rm -r fastq_files
#rm -r QC_results

date
