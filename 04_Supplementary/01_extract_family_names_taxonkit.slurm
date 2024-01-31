#!/bin/bash
#SBATCH --job-name=taxID2family    # Job name
#SBATCH --account=group_name			# your own sponser or account
#SBATCH --qos=group_name                # your own sponser or account
#SBATCH --mail-type=ALL               	# Mail events
#SBATCH --mail-user=name@email.com      # Where to send email	
#SBATCH --ntasks=1                    	# Run on a single machine (node)
#SBATCH --cpus-per-task 1		  		# Run on a several CPU
#SBATCH --mem-per-cpu=1gb             # Memory limit
#SBATCH --time=01:00:00               	# Time: hrs:min:sec
#SBATCH --output=taxID2family_output_%j.out  	 	# Output and error log 

pwd; hostname; date 

# Load dependencies
module purge
module load taxonkit/0.13.0

# General options
max_cpus=4
WD="/path/to/RefSeq/phage/taxon_name_information"

# Navigate to the working directory
cd ${WD}

# Define input and output file names
input_file="acc_NCBI_taxa.txt" #file containing phage accession numbers and corresponding NCBI taxa IDs
output_file="NCBI_taxa_families.txt"

# Check if the output file exists and remove it to avoid duplicates
if test -f ${output_file}; then
    rm ${output_file}
fi

printf "Extracting family from NCBI taxon IDs along with accession numbers\n"

# Extract TaxIds and accession numbers, then use taxonkit to get family
awk -F '\t' '{print $1}' $input_file | \
taxonkit reformat -I 1 -f "{f}" -j ${max_cpus} | \
paste <(awk -F '\t' '{print $2}' $input_file) - > ${output_file}

# Finish and print date
date
