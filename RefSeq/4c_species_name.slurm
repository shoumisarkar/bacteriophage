#!/bin/bash
#SBATCH --job-name=taxID2species    # Job name
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
module load taxonkit/0.13.0

# General options
max_cpus=4
WD="/path/to/RefSeq/phage/taxon_names"

# checks if output file exists and remove. This avoids replicates
if test -f NCBI_taxa_species.txt; then
	rm NCBI_taxa_species.txt
fi	

printf "Retrieve species name corresponding NCBI taxon ID \n"

cd ${WD}
input_file="acc_NCBI_taxa.txt"

# selects the first column of the input_file
cut -f 1 -d $'\t' $input_file >> NCBI_taxa.txt

# retrieve the species
taxonkit lineage --data-dir ${WD} NCBI_taxa.txt -j ${max_cpus} >> NCBI_taxa_species.txt

# rm temp file
rm NCBI_taxa.txt

date
