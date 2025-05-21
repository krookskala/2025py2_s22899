from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt

def main():
    Entrez.email = input("Email: ")
    Entrez.api_key = input("API key: ")
    taxid = input("TaxID: ")
    min_len = int(input("Min length: "))
    max_len = int(input("Max length: "))

    try:
        h = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism]", usehistory="y")
        r = Entrez.read(h)
        count, webenv, query_key = int(r["Count"]), r["WebEnv"], r["QueryKey"]
    except Exception as e:
        print("Search error:", e)
        return

    records, step, limit = [], 100, 1000
    for start in range(0, min(count, limit), step):
        try:
            h = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                              retstart=start, retmax=step, webenv=webenv, query_key=query_key)
            for rec in SeqIO.parse(h, "gb"):
                l = len(rec.seq)
                if min_len <= l <= max_len:
                    records.append({"Accession": rec.id, "Length": l, "Description": rec.description})
            h.close()
        except:
            continue

    if not records:
        print("No records in specified range.")
        return

    df = pd.DataFrame(records).sort_values("Length", ascending=False)
    df.to_csv("records.csv", index=False)

    plt.figure(figsize=(12, 6))
    plt.plot(df["Accession"], df["Length"], marker='o')
    plt.xticks(range(0, len(df), 10), df["Accession"].iloc[::10], rotation=45)
    plt.tight_layout()
    plt.savefig("sequence.png")

    print("Saved: records.csv and sequence.png")

if __name__ == "__main__":
    main()
