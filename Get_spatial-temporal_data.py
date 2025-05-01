from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

# Set Entrez email
Entrez.email = "3331606993@qq.com"

# Define search criteria
tax_id = "39054"
start_date = "1000/01/01"
end_date = "2023/12/31"
search_term = f"txid{tax_id}[Organism:exp] AND {start_date}[PDAT] : {end_date}[PDAT]"

# Search NCBI database
handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1, usehistory="y")
record = Entrez.read(handle)
handle.close()

# Get total record count, WebEnv, and query_key
total_count = int(record["Count"])
webenv = record["WebEnv"]
query_key = record["QueryKey"]

# Print total record count
print(f"Total number of sequences found: {total_count}")

if total_count == 0:
    print("No sequences found. Please check your search term.")
else:
    # Define the number of records to download per batch
    batch_size = 10000

    # Initialize sequence list
    sequences = []

    # Download sequences in batches and save to a FASTA file
    with open("ev7111111_sequences.fasta", "w") as out_handle:
        for start in tqdm(range(0, total_count, batch_size), desc="Downloading batches"):
            end = min(start + batch_size, total_count)
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
            out_handle.write(handle.read())
            handle.close()
            print(f"Downloaded and saved sequences from {start} to {end}")

    print(f"\nData saved to ev71_sequences.fasta")

    # Download sequences again in batches to extract month information
    for start in tqdm(range(0, total_count, batch_size), desc="Processing sequences for month extraction"):
        end = min(start + batch_size, total_count)
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
        records = SeqIO.parse(handle, "genbank")
        for record in records:
            # Extract month information
            collection_date = None
            country = None
            for feature in record.features:
                if feature.type == "source":
                    if "collection_date" in feature.qualifiers:
                        date_str = feature.qualifiers["collection_date"][0]
                        try:
                            collection_date = pd.to_datetime(date_str)
                        except ValueError:
                            collection_date = None
                    if "geo_loc_name" in feature.qualifiers:
                        country = feature.qualifiers["geo_loc_name"][0]
            
            if collection_date is not None:
                # Check if the date contains month information and is before December 31, 2023
                if collection_date <= pd.to_datetime('2023-12-31'):
                    if collection_date.month != 1 or collection_date.day != 1:
                        year_month = collection_date.strftime('%Y-%m')                    
                        sequences.append({
                                "Accession": record.id,
                                "Sequence": str(record.seq),
                                "Month": year_month,
                                "Nation": country
                            })
                        print(f"Added sequence: {record.id}, Year-Month: {year_month}, Nation: {country}")  # Debug information
        handle.close()

# Save results to a DataFrame
df = pd.DataFrame(sequences)

# Print results
print("\nProcessed DataFrame:")
print(df)

# Save DataFrame to an Excel file
if not df.empty:
    df.to_excel("ev71_sequences.xlsx", index=False)
    print("Data saved to ev71_sequences.xlsx")
else:
    print("No data to save to Excel.")

    # Plot month distribution
    if not df.empty:
        plt.figure(figsize=(12, 6))
        df['Month'].value_counts().sort_index().plot(kind='bar')
        plt.title('Distribution of CVA6 Sequences by Month (China)')
        plt.xlabel('Year-Month')
        plt.ylabel('Number of Sequences')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.show()
    else:
        print("No data to plot.")
