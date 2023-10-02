#!/bin/bash
#SBATCH --job-name=KMAindex    # Job name
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

# loads dependencies
module purge
module load kma/20230106

# General options
phageDB="/path/to/RefSeq/phage/database"
output_dir="/path/to/RefSeq/kma_indexing/output"

printf "Index the database for the combined sequences of bacteriophages \n"

cd ${phageDB}

for seqID in $(ls *_comb_refseq.fa | sed -e 's/_comb_refseq.fa//' | sort -u)
do
	echo "KMA indexing for ${seqID}"
	kma index -i ${phageDB}/${seqID}_comb_refseq.fa -o ${output_dir}/${seqID}_comb
done

date
