#!/bin/bash
#SBATCH --job-name=bp_taxID2family    # Job name
#SBATCH --account=group_name			# your own sponser or account
#SBATCH --qos=group_name                # your own sponser or account
#SBATCH --mail-type=ALL               	# Mail events
#SBATCH --mail-user=name@email.com      # Where to send email	
#SBATCH --ntasks=1                    	# Run on a single machine (node)
#SBATCH --cpus-per-task 1		  		# Run on a several CPU
#SBATCH --mem-per-cpu=1gb             # Memory limit
#SBATCH --time=01:00:00               	# Time: hrs:min:sec
#SBATCH --output=log_bp_taxID2family.out  # Output and error log 

pwd; hostname; date 

# Load Python module
module load python

# Run the script
python fetch_phage_families.py

# Finish and print date
date
