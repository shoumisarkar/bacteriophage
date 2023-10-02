#!/bin/sh
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

# loads dependencies
module purge
module load edirect/12.2

echo "Running Kaiju script for multiple fastq files on several CPU cores" 

# General options
max_cpus=1
phageDB="/path/to/download/RefSeq/phage/database/"
cd ${phageDB}

# Download reference sequences of bacteriophage from NCBI DB of viral genomes. Filter hosts (by bacteria and archaea)
# https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239
# get the accession numbers and download the reference sequence corresponding to each accession number.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
printf "(1) Search for all bacteriophage (filtered by archaea) reference genomes in NCBI's nucleotide database \n"
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# checks if a the accession list exists and remove. This avoids replicates
if test -f bacteriophage_archaea_accessions.txt; then
	rm bacteriophage_archaea_accessions.txt
fi	

# get the accession numbers
esearch -db nucleotide -query '((Viruses[Organism]) AND (srcdb_refseq[PROP]) NOT wgs[PROP] NOT cellular organisms[ORGN] NOT AC_000001:AC_999999[PACC] AND ("vhost archaea"[Filter]))' | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion >> bacteriophage_archaea_accessions.txt

# Download the reference genome sequence for each bacteriophage
touch bacteriophage_archaea_refseq.fa # create empty/null file

input_file="bacteriophage_archaea_accessions.txt"
for acc in $(cat "$input_file")
do
  # Check if the line starts with "<"
  if [[ ${acc:0:1} != "<" ]]
  then
	# Process the line here
    echo "Processing line: $acc"
	
	# Download the reference genome sequence
    esearch -db nucleotide -query "$acc" | efetch -format fasta_cds_na >> ${acc}_archaea.fa
  fi
done

# concatenate all reference sequence into a single file
cat *archaea.fa* > bacteriophage_archaea_refseq.fa
# remove the individual reference sequences
rm *archaea.fa* 

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
printf "(2) Search for all bacteriophage (filtered by bacteria) reference genomes in NCBI's nucleotide database \n"
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# checks if a the accession list exists and remove. This avoids replicates
if test -f bacteriophage_bacteria_accessions.txt; then
	rm bacteriophage_bacteria_accessions.txt
fi	

# get the accession numbers
esearch -db nucleotide -query '((Viruses[Organism]) AND (srcdb_refseq[PROP]) NOT wgs[PROP] NOT cellular organisms[ORGN] NOT AC_000001:AC_999999[PACC] AND ("vhost bacteria"[Filter]))' | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion >> bacteriophage_bacteria_accessions.txt

# Download the reference genome sequence for each bacteriophage
touch bacteriophage_bacteria_refseq.fa # create empty/null file

input_file="bacteriophage_bacteria_accessions.txt"
for acc in $(cat "$input_file")
do
  # Check if the line starts with "<"
  if [[ ${acc:0:1} != "<" ]]
  then
	# Process the line here
    echo "Processing line: $acc"
	
	# Download the reference genome sequence
    esearch -db nucleotide -query "$acc" | efetch -format fasta_cds_na >> ${acc}_bacteria.fa
  fi
done

# concatenate all reference sequence into a single file
cat *bacteria.fa* > bacteriophage_bacteria_refseq.fa
# remove the individual reference sequences
rm *bacteria.fa* 

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
printf "(3) Concantenate the both files into a single file of bacteriophage reference genomes \n"
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cat bacteriophage_archaea_refseq.fa bacteriophage_bacteria_refseq.fa > bacteriophage_comb_refseq.fa

date


