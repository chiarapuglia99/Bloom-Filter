from Bio import SeqIO
import pandas as pd

# Percorso del file FASTQ caricato
input_file = "/Progetti/Bloom-Filter/balaenoptera.fastq"
output_file = "/Progetti/Bloom-Filter/balaenoptera.csv"

# Parsing del FASTQ
records = []
for record in SeqIO.parse(input_file, "fastq"):
    records.append({
        "id": record.id,
        "sequence": str(record.seq),
        "quality": record.letter_annotations["phred_quality"]
    })

# Creazione DataFrame
df = pd.DataFrame(records)

# Salvataggio in CSV
df.to_csv(output_file, index=False)

print(f"File CSV salvato in: {output_file}")
