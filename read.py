from Bio import SeqIO
import pandas as pd
import gzip

records = []
with gzip.open(r"./data/uniprot_sprot.fasta.gz", "rt", encoding="utf-8") as file:
    for record in SeqIO.parse(file, "fasta"):
        records.append({
            'id': record.id,
            'description': record.description,
            'sequence': str(record.seq)
        })

df = pd.DataFrame(records)
print(f"Loaded {len(df)} records")
print(df.head())
df.to_csv(r"./data/uniprot_data.csv", index=False)
print("Data saved as .\\data\\uniprot_data.csv")
