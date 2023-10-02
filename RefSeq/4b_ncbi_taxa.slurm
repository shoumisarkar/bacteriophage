#!/bin/bash
#SBATCH --job-name=acc2taxIDs    # Job name
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
module load edirect/12.2

# General options
WD="/path/to/RefSeq/phage/taxon_names"

printf "Read the .txt file of accession numbers and get corresponding NCBI taxa ID \n"

cd ${WD}


# checks if output file exists and remove. This avoids replicates
if test -f acc_NCBI_taxa.txt; then
	rm acc_NCBI_taxa.txt
fi	

separator="\t"

input_file="uniq_accs.txt"
for acc in $(cat "$input_file")
do

	# get the NCBI taxa ID using esearch, stores the result in res
	res=$(esearch -db nucleotide -query $acc | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId)
	
	# append the acc and corresponding NCBI taxon id to file
	acc_res="$res$separator$acc"
	echo -e "$acc_res" >> acc_NCBI_taxa.txt
done

date
