from Bio import SeqIO

def filter_protein_fasta(input_file, output_file):
    valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy')
    with open(output_file, 'w') as outfile:
        for record in SeqIO.parse(input_file, 'fasta'):
            if all(aa in valid_amino_acids for aa in record.seq):
                SeqIO.write(record, outfile, 'fasta')
            else:
                # Identify invalid amino acids
                invalid_amino_acids = set(aa for aa in record.seq if aa not in valid_amino_acids)
                print(f"Deleted sequence ID: {record.id}")
                print(f"Invalid amino acids: {invalid_amino_acids}")
                print(f"Sequence: {record.seq}\n")

# Example usage
input_file = 'sequences-CVA6-VP1aa.fas'
output_file = 'sequences-CVA6-VP1aa_modified.fas'
filter_protein_fasta(input_file, output_file)