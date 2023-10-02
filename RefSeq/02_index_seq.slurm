#!/bin/bash
#SBATCH --job-name=KMAindex    # Job name
#SBATCH --account=biostat-dept            # your own sponser or account from this class
#SBATCH --qos=biostat-dept                 # your own sponser or account from this class
#SBATCH --mail-type=ALL               # Mail events
#SBATCH --mail-user=sanyasosamuel@ufl.edu          # Where to send email	
#SBATCH --ntasks=1                    # Run on a single machine (node)
#SBATCH --cpus-per-task 1		  # Run on a several CPU
#SBATCH --mem-per-cpu=50gb                    # Memory limit
#SBATCH --time=2:00:00               # Time: hrs:min:sec
#SBATCH --output=log_KMA_index.out   # Output and error log 

pwd; hostname; date 

# loads dependencies
module purge
module load kma/20230106

# General options
phageDB="/blue/somnath.datta/share/RefSeq/phage"
output_dir="/blue/somnath.datta/share/RefSeq/phage/kma_indexing"

printf "Index the database for the combined sequences of bacteriophages \n"

cd ${phageDB}

for seqID in $(ls *_comb_refseq.fa | sed -e 's/_comb_refseq.fa//' | sort -u)
do
	echo "KMA indexing for ${seqID}"
	kma index -i ${phageDB}/${seqID}_comb_refseq.fa -o ${output_dir}/${seqID}_comb
done

date
