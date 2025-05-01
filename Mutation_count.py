from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import re
import time
import os
from operator import itemgetter

# Input file locations (nucleotide as filename3 and amino acid as filename)
# ——————————————————————————————————————————————————————————————————————————————————#
reference = "AF081297.1"  # Reference strain
filename = "sequences-CVA6-VP1aa_modified.fas"  # Amino acid
filename3 = "sequences-CVA6-VP1.fasta"  # Nucleotide
# ——————————————————————————————————————————————————————————————————————————————————#

def filter_high_mutation_sequences1(filename, reference_id=reference, threshold=0.5):  # Remove amino acid sequences with mutation rate above threshold
    records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    if reference_id not in records:
        raise ValueError(f"Reference sequence {reference_id} not found in file.")
    reference = records[reference_id].seq
    filtered_records = []
    for record_id, record in records.items():
        if record_id != reference_id:
            seq = record.seq
            mutations = sum(1 for a, b in zip(seq, reference) if a != b)
            if mutations / len(seq) <= threshold:
                filtered_records.append(record)
        else:
            filtered_records.append(record)
    SeqIO.write(filtered_records, "filtered1.fasta", "fasta")

def filter_high_mutation_sequences2(filename, reference_id=reference, threshold=0.5):  # Remove nucleotide sequences with mutation rate above threshold
    records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    if reference_id not in records:
        raise ValueError(f"Reference sequence {reference_id} not found in file.")
    reference = records[reference_id].seq
    filtered_records = []
    for record_id, record in records.items():
        if record_id != reference_id:
            seq = record.seq
            mutations = sum(1 for a, b in zip(seq, reference) if a != b)
            if mutations / len(seq) <= threshold:
                filtered_records.append(record)
        else:
            filtered_records.append(record)
    SeqIO.write(filtered_records, "filtered2.fasta", "fasta")

def uppercase(input_fasta, output_fasta):  # Convert lowercase sequences to uppercase
    records = SeqIO.parse(input_fasta, "fasta")
    with open(output_fasta, "w") as output_handle:
        for record in records:
            new_seq = str(record.seq).translate(str.maketrans('abcdefghijklmnopqrstuvwxyz', 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'))
            record.seq = Seq(new_seq)
            SeqIO.write(record, output_handle, "fasta")

def cleaner(filename):  # Remove sequences containing non-ACGT bases
    input_fasta = filename
    output_fasta = "output.fasta"
    records = list(SeqIO.parse(input_fasta, "fasta"))
    valid_records = [
        record for record in records
        if all(base in "ACGT" for base in str(record.seq))
    ]
    SeqIO.write(valid_records, output_fasta, "fasta")

def print_sequence_lengths(filename):  # Get sequence lengths
    records = SeqIO.parse(filename, "fasta")
    for record in records:
        length = len(record)
    return length

def count_base_frequencies(filename):  # Count amino acid frequencies at each position
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y']
    base_counts = {aa: [0] * length_count for aa in amino_acids}
    records = SeqIO.parse(filename, "fasta")
    num_sequences = 0  
    for record in records:
        seq = str(record.seq)      
        if len(seq) == length_count:
            for i, base in enumerate(seq):
                if base in base_counts:
                    base_counts[base][i] += 1
            num_sequences += 1
    base_frequencies = {base: [count / num_sequences for count in counts] for base, counts in base_counts.items()}
    return base_frequencies

def count_base_frequencies2(filename):  # Count nucleotide frequencies at each position
    amino_acids = ['A', 'C', 'G', 'T']
    base_counts = {aa: [0] * length_count2 for aa in amino_acids}
    records = SeqIO.parse(filename, "fasta")
    num_sequences = 0 
    for record in records:
        seq = str(record.seq)
        if len(seq) == length_count2:
            for i, base in enumerate(seq):
                if base in base_counts:
                    base_counts[base][i] += 1
            num_sequences += 1
    base_frequencies2 = {base: [count / num_sequences for count in counts] for base, counts in base_counts.items()}
    return base_frequencies2

def calculate_amino_acid_variants(filename, reference_id=reference):  # Calculate amino acid mutation rates
    amino_acid_variants = [{amino_acid: 0 for amino_acid in "ACDEFGHIKLMNPQRSTVWXY"} for _ in range(length_count)]
    records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    if reference_id not in records:
        raise ValueError(f"Reference sequence {reference_id} not found in file.")
    reference = records[reference_id].seq
    for record_id, record in records.items():
        if record_id != reference_id:
            seq = record.seq
            for i in range(len(seq)):
                if seq[i] != reference[i]:  
                    amino_acid_variants[i][seq[i]] += 1
    variant_counts = [sum(1 for aa in amino_acids if amino_acids[aa] > 0) for amino_acids in amino_acid_variants]
    df = pd.DataFrame({'Variant Count': variant_counts})
    df.index = range(1, len(variant_counts) + 1)
    return df

def calculate_amino_acid_variants2(filename, reference_id=reference):  # Calculate nucleotide mutation rates
    amino_acid_variants = [{amino_acid: 0 for amino_acid in "ACGTYRSWHNMK"} for _ in range(length_count2)]
    records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    if reference_id not in records:  # Ensure reference sequence exists
        raise ValueError(f"Reference sequence {reference_id} not found in file.")
    reference = records[reference_id].seq
    for record_id, record in records.items():
        if record_id != reference_id:
            seq = record.seq
            for i in range(len(seq)):
                if seq[i] != reference[i]:  
                    amino_acid_variants[i][seq[i]] += 1
    variant_counts = [sum(1 for aa in amino_acids if amino_acids[aa] > 0) for amino_acids in amino_acid_variants]
    df = pd.DataFrame({'Variant Count': variant_counts})
    df.index = range(1, len(variant_counts) + 1)
    return df

def normalize_date(date_tuple):  # Normalize date
    year, month, day = date_tuple
    return (year, month or 1, day or 1)

def time_order(file_path):  # Order by time and accumulate mutation sites
    records = list(SeqIO.parse(file_path, "fasta"))
    date_record_pairs = []
    for record in records:
        match = re.search(r'\|(\d{4})(?:-(\d{2}))?(?:-(\d{2}))?', record.description)
        if match:
            year = int(match.group(1))
            if year >= 1970:  # Filter records before 1970
                month = int(match.group(2)) if match.group(2) else None
                day = int(match.group(3)) if match.group(3) else None
                date_info = (year, month, day)  # Convert date information to tuple for sorting
                date_record_pairs.append((date_info, record))
                date_info = (year, month, day)
                date_record_pairs.append((date_info, record))
    normalized_date_record_pairs = [(normalize_date(date_info), record) for date_info, record in date_record_pairs]
    sorted_pairs = sorted(normalized_date_record_pairs, key=lambda x: x[0])
    SeqIO.write([record for _, record in sorted_pairs], "ordered.fasta", "fasta")

    time.sleep(1)

    records = list(SeqIO.parse("ordered.fasta", "fasta"))
    reference = records[0]
    year_variants = {}
    date_pattern = r'\|(\d{4})[-/]?(\d{2})?[-/]?(\d{2})?'  # Regular expression to match date format
    for record in records:
        match = re.search(date_pattern, record.description)
        if match:
            year = int(match.group(1))
            month = match.group(2)
            day = match.group(3)
            if not month and not day:  # If no month and day, use only year
                year_key = year
            elif month and not day:  # If month is present, use year-month as key
                year_key = f"{year}-{month}"
            else:  # If complete date is present, use year-month-day as key
                year_key = f"{year}-{month}-{day}"
            variants = [i for i, (a, b) in enumerate(zip(reference.seq, record.seq)) if a != b]  # Compare sequences to find variant sites
            year_key = str(year)  # Ensure key is always a string type
            if year_key in year_variants:
                year_variants[year_key].update(set(variants))
            else:
                year_variants[year_key] = set(variants)
    cumulative_variants = {}  # Calculate cumulative variants
    total_variants = set()
    for year_key in sorted(year_variants.keys()):  # Now all keys are strings, can sort normally
        total_variants.update(year_variants[year_key])
        cumulative_variants[year_key] = len(total_variants)

    df = pd.DataFrame(list(cumulative_variants.items()), columns=['Date', 'NumVariants'])
    df['Variability'] = df['NumVariants'] / print_sequence_lengths(file_path)

    return df

uppercase(filename, 'uppered.fasta')  # Convert lowercase sequences to uppercase
uppercase(filename3, 'uppered1.fasta')
filename = os.path.abspath("uppered.fasta")

cleaner("uppered1.fasta")  # Remove sequences containing non-ACGT bases
filename2 = os.path.abspath("output.fasta")

filter_high_mutation_sequences1(filename)  # Remove amino acid sequences with mutation rate above threshold
filter_high_mutation_sequences2(filename2)
filename = os.path.abspath("filtered1.fasta")
filename2 = os.path.abspath("filtered2.fasta")

length_count = print_sequence_lengths(filename)  # Get sequence lengths
length_count2 = print_sequence_lengths(filename2)

base_frequencies = count_base_frequencies(filename)  # Count amino acid frequencies at each position
base_frequencies2 = count_base_frequencies2(filename2)
df = pd.DataFrame.from_dict(base_frequencies)  # Generate list file
df2 = pd.DataFrame.from_dict(base_frequencies2)
df.columns = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y']
df.index = range(1, df.shape[0] + 1)
df2.columns = ['A', 'C', 'G', 'T']
df2.index = range(1, df2.shape[0] + 1)

mutation_rates_df = calculate_amino_acid_variants(filename)  # Calculate amino acid mutation rates
mutation_rates_df2 = calculate_amino_acid_variants2(filename2)

gf = time_order(filename)  # Order by time and accumulate mutation sites
gf2 = time_order(filename2)

# Write DataFrame to Excel file
with pd.ExcelWriter('amino acid mutation.xlsx', engine='openpyxl') as writer:
    df.to_excel(writer, sheet_name='frequency')
    mutation_rates_df.to_excel(writer, sheet_name='mutation ratio')
    gf.to_excel(writer, sheet_name='accumulate mutation rate')

with pd.ExcelWriter('nucleic acid mutation.xlsx', engine='openpyxl') as writer:
    df2.to_excel(writer, sheet_name='frequency')
    mutation_rates_df2.to_excel(writer, sheet_name='mutation ratio')
    gf2.to_excel(writer, sheet_name='accumulate mutation rate')

print('Output success')