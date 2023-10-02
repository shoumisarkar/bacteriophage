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

# General options
max_cpus=40
readPath="/path/to/bbmap_output_files"
phageDBindex="/path/to/kma_indexing/output"
output_dir="/path/to/RefSeq/phage/map_results"

separator="_"

# we have an indexed combined database of bacteriophages (bacteria and archaea viruses)
# goal: map pre-processed metagenomic reads to this database

cd ${readPath}
for seqID in $(ls *.fastq | sed -e 's/_1.fastq//' -e 's/_2.fastq//' | sort -u)
do
	printf "KMA mapping for ${seqID} \n"
	cd ${phageDBindex}
	
	for dbID in $(ls *.seq.b | sed -e 's/.seq.b//' | sort -u)
	do
		out_name="$seqID$separator$dbID" # name for storing the output
		
		kma -ipe ${readPath}/${seqID}_1.fastq ${readPath}/${seqID}_2.fastq \
			-t_db ${phageDBindex}/${dbID} -t ${max_cpus} \
			-o ${output_dir}/$out_name -1t1 -ef -verbose
	done

done


# the 7th row of the .mapstat files has the column names begin with '# '. 
# R will not recognize this row since it begins with '#'
# remove the first two characters from the row
cd ${output_dir}
sed -i '7s/^..//' *.mapstat

date
