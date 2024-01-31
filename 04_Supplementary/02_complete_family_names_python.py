import time
from Bio import Entrez
import pandas as pd

# Update this email to your email address
Entrez.email = "name@email.com"  

def fetch_family(accession_number):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_number, retmode="xml")
        record = Entrez.read(handle)
        taxonomy = record[0]['GBSeq_taxonomy'].split('; ')
        # Search for the family name in the taxonomy list
        family = next((taxon for taxon in taxonomy if taxon.endswith('viridae')), "Unknown")
        return family
    except Exception as e:
        print(f"Error fetching data for {accession_number}: {e}")
        return "Error"

def main():
    # Replace 'NCBI_taxa_families.xlsx' with the path to your Excel file
    phage_data = pd.read_excel('NCBI_taxa_families.xlsx', header=None)
    phage_data.columns = ['Accession', 'TaxonID', 'Family']

    total = len(phage_data)
    print(f"Starting to fetch family information for {total} phages.")

    for index, row in phage_data.iterrows():
        if pd.isna(row['Family']):
            print(f"Fetching family for accession number: {row['Accession']} ({index+1}/{total})")
            family = fetch_family(row['Accession'])
            phage_data.at[index, 'Family'] = family
            print(f"Updated {row['Accession']} with family: {family}")
            time.sleep(0.5)  # Sleep to avoid overloading NCBI servers
        else:
            print(f"Family already present for {row['Accession']} ({index+1}/{total})")

    # Replace 'updated_NCBI_taxa_families.xlsx' with your desired output file name
    phage_data.to_excel('updated_phage_data_families.xlsx', index=False)
    print("Family information update complete.")

if __name__ == "__main__":
    main()
