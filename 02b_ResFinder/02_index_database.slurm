#!/bin/bash
#SBATCH --job-name=KMAindex    # Job name
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
ResFinderDB="/path/to/ResFinder/database"
output_dir="/path/to/ResFinder/kma_indexing/output"

printf "Index the database \n"

cd ${ResFinderDB}

for seqID in $(ls *.fsa | sed -e 's/.fsa//' | sort -u)
do
	echo "KMA indexing for ${seqID}"
	kma index -i ${ResFinderDB}/${seqID}.fsa -o ${output_dir}/${seqID}
done

date
