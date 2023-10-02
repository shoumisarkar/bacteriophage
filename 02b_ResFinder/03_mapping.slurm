#!/bin/bash
#SBATCH --job-name=KMAmap    # Job name
#SBATCH --account=group_name			# your own sponser or account
#SBATCH --qos=group_name                # your own sponser or account
#SBATCH --mail-type=ALL               	# Mail events
#SBATCH --mail-user=name@email.com      # Where to send email	
#SBATCH --ntasks=1                    	# Run on a single machine (node)
#SBATCH --cpus-per-task 30		  		# Run on a several CPU
#SBATCH --mem-per-cpu=2gb             	# Memory limit
#SBATCH --time=01:00:00               	# Time: hrs:min:sec
#SBATCH --output=output_%j.out   		# Output and error log 

pwd; hostname; date 

# loads dependencies
module purge
module load kma/20230106

echo "Running Kaiju script for multiple fastq files on several CPU cores" 

# General options
max_cpus=30
readPath="/path/to/bbmap_output_files"
ResFinderDBindex="/path/to/ResFinder/kma_indexing/output"
output_dir="/path/to/ResFinder/map_results"
separator="_"

# For each class of AMR gene, we have a different database. We use a nested for loop
# to run each pair-end sample against each database index.

cd ${readPath}
for seqID in $(ls *.fastq | sed -e 's/_1.fastq//' -e 's/_2.fastq//' | sort -u)
do
	printf "KMA mapping for ${seqID} \n"
	cd ${ResFinderDBindex}
	
	for dbID in $(ls *.length.b | sed -e 's/.length.b//' | sort -u)
	do
		out_name="$seqID$separator$dbID" # name for storing the output
		
		kma -ipe ${readPath}/${seqID}_1.fastq ${readPath}/${seqID}_2.fastq \
			-t_db ${ResFinderDBindex}/${dbID} -t ${max_cpus} \
			-o ${output_dir}/$out_name -1t1 -ef -verbose
		
		# remove redundant files
		rm ${output_dir}/${out_name}.aln
		rm ${output_dir}/${out_name}.fsa
		rm ${output_dir}/${out_name}.frag.gz
		rm ${output_dir}/${out_name}.res

	done
done

# the 7th row of the .mapstat files has the column names begin with '# '. 
# R will not recognize this row since it begins with '#'
# remove the first two characters from the row
cd ${output_dir}
sed -i '7s/^..//' *.mapstat

date

