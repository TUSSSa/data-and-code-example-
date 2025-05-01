from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

# 设置Entrez邮箱
Entrez.email = #"yourmail.com"

# 定义搜索条件
tax_id = "39054"
start_date = "1000/01/01"
end_date = "2023/12/31"
search_term = f"txid{tax_id}[Organism:exp] AND {start_date}[PDAT] : {end_date}[PDAT]"

# 搜索NCBI数据库
handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1, usehistory="y")
record = Entrez.read(handle)
handle.close()

# 获取总记录数、WebEnv 和 query_key
total_count = int(record["Count"])
webenv = record["WebEnv"]
query_key = record["QueryKey"]

# 打印总记录数
print(f"Total number of sequences found: {total_count}")

if total_count == 0:
    print("No sequences found. Please check your search term.")
else:
    # 定义每批次下载的数量
    batch_size = 10000

    # 初始化序列列表
    sequences = []

    # 分批下载序列并保存到FASTA文件
    with open("ev7111111_sequences.fasta", "w") as out_handle:
        for start in tqdm(range(0, total_count, batch_size), desc="Downloading batches"):
            end = min(start + batch_size, total_count)
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
            out_handle.write(handle.read())
            handle.close()
            print(f"Downloaded and saved sequences from {start} to {end}")

    print(f"\nData saved to ev71_sequences.fasta")

    # 再次分批下载序列以提取月份信息
    for start in tqdm(range(0, total_count, batch_size), desc="Processing sequences for month extraction"):
        end = min(start + batch_size, total_count)
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
        records = SeqIO.parse(handle, "genbank")
        for record in records:
            # 提取月份信息
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
                # 检查日期是否包含月份信息并且在2023年12月31日之前
                if collection_date <= pd.to_datetime('2023-12-31'):
                    if collection_date.month != 1 or collection_date.day != 1:
                        year_month = collection_date.strftime('%Y-%m')                    
                        sequences.append({
                                "Accession": record.id,
                                "Sequence": str(record.seq),
                                "Month": year_month,
                                "Nation": country
                            })
                        print(f"Added sequence: {record.id}, Year-Month: {year_month}, Nation: {country}")  # 调试信息
        handle.close()

# 将结果保存到DataFrame
df = pd.DataFrame(sequences)

# 打印结果
print("\nProcessed DataFrame:")
print(df)

# 将DataFrame保存到Excel文件
if not df.empty:
    df.to_excel("ev71_sequences.xlsx", index=False)
    print("Data saved to ev71_sequences.xlsx")
else:
    print("No data to save to Excel.")


    # 绘制月份分布图
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
