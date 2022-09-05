from lib2to3.pgen2.pgen import generate_grammar
from tokenize import Name
from Bio import SeqIO
import pandas as pd

print("hi")
if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        records = list(SeqIO.parse(snakemake.input["fasta"], "fasta"))
        seq_df = pd.read_csv(snakemake.input["csv"], index_col=0)
        for each in records:
            gid = each.id.split("::")[0]
            each.id = seq_df.loc[gid, "name"]
            each.name = seq_df.loc[gid, "name"]
        SeqIO.write(records, snakemake.output["fasta"], "fasta")
