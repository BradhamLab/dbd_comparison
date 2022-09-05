import re
import pandas as pd
from Bio import SeqIO


def merge_domains(data):
    # Sort the array on the basis of start values of intervals.
    data.sort_values("seq_from", inplace=True)
    intervals = data[["seq_from", "seq_to"]].values
    stack = []
    info = []
    # insert first interval into stack
    stack.append(intervals[0])

    def get_info(i_idx):
        return data.iloc[i_idx, :][["query_name", "evalue"]]

    info.append(get_info(0))

    for i, start_stop in enumerate(intervals[1:]):
        # Check for overlapping interval,
        # if interval overlap
        new_info = get_info(i + 1)
        if stack[-1][0] <= start_stop[0] <= stack[-1][-1]:
            stack[-1][-1] = max(stack[-1][-1], start_stop[-1])
            if new_info["evalue"] < info[-1]["evalue"]:
                info[-1] = new_info
        else:
            stack.append(start_stop)
            info.append(new_info)
    merged_df = pd.DataFrame(info)
    merged_df[["seq_from", "seq_to"]] = stack
    return merged_df.rename(columns={"query_name": "domain"})


def extract_domain_sequence(sequences, data):
    entry = sequences[data["gene"]][data["seq_from"] : data["seq_to"]]
    entry.id += f':{data["seq_from"]}-{data["seq_to"]}_{data["domain"]}'
    entry.name = entry.id
    entry.description = ""
    return entry


if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        df = pd.read_csv(snakemake.input["sequence_domains"]).query(
            f"evalue < {snakemake.params['e_val']}"
        ).copy()
        allowed_domains = []
        for each in snakemake.input["domain_list_files"]:
            with open(each, "r") as f:
                allowed_domains += [x.strip() for x in f.readlines()]
        regex = [re.compile(x) for x in snakemake.params["domain_regex"]]

        def domain_match(x):
            return any([pattern.search(x) is not None for pattern in regex])

        df["allowed"] = df.query_name.isin(allowed_domains) | df.query_name.apply(
            lambda x: domain_match(x)
        )
        groups = dict(list(df.query("allowed").groupby("target_name")))
        summarized = []
        for key, data in groups.items():
            merged = merge_domains(data)
            merged["gene"] = key
            summarized.append(merged)
        retained = pd.concat(summarized)[["gene", "domain", "seq_from", "seq_to"]]
        retained.to_csv(snakemake.output[0], index=False)

        df[~df.target_name.isin(retained.gene)].query("not allowed").sort_values(
            ["target_name", "seq_from"]
        ).to_csv(snakemake.output["filtered"], index=False)

        # extract sequences
        records = {x.name: x for x in SeqIO.parse("output/tf_proteins.fasta", "fasta")}
        out_fasta = []
        for idx in retained.index:
            domain_seq = extract_domain_sequence(records, retained.loc[idx, :])
            out_fasta.append(domain_seq)
        SeqIO.write(out_fasta, snakemake.output["fasta"], "fasta")
