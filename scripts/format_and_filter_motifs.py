import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import logging
from pathlib import Path
from Bio import motifs


def scale_counts(pwm, scalar=100):
    for nuc in "ACGT":
        pwm.counts[nuc] = [x * scalar for x in pwm.counts[nuc]]


def pident_plot(df, threshold=0.7):
    with sns.plotting_context("poster"):
        fig = plt.figure(figsize=(6, 6))
        plot = sns.displot(
            data=dbd_df.sort_values("TF"),
            x="pident",
            hue="TF",
            col="TF",
            col_wrap=5,
            facet_kws={"sharey": False, "sharex": True},
            legend=False,
            bins=12,
        )
        for ax in plot.axes:
            ax.axvline(x=threshold, ls="--", c="black", lw=2)
            ax.yaxis.set_major_formatter(FormatStrFormatter("%d"))


HEADER = [
    "#motif_id",
    "motif_name",
    "motif_description",
    "source_name",
    "source_version",
    "gene_name",
    "motif_similarity_qvalue",
    "similar_motif_id",
    "similar_motif_description",
    "orthologous_identity",
    "orthologous_gene_name",
    "orthologous_species",
    "description",
]

ENTRY = [
    "{}",
    "{}",
    "{}",
    "jaspar",
    "1.1",
    "{}",
    "0.000000",
    "None",
    "None",
    "1.000000",
    "None",
    "None",
    "orthology",
]


def get_line(motif_id, gene_name, consensus):
    return (
        "\t".join(ENTRY).format(
            motif_id, "_".join(motif.split("_")[2:]), consensus, gene_name
        )
        + "\n"
    )


if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        logging.basicConfig(filename=snakemake.log[0], level=logging.WARNING)
        lookup = pd.read_csv(snakemake.input["id2name"])
        dbd_df = pd.read_csv(snakemake.input["csv"], index_col=0).assign(
            TF=lambda df_: df_.query_id.map(lookup.set_index("gene_id").name).apply(
                lambda x: x.replace("Sp-", "Lv-")
            )
        )
        # plot distribution pident of DBDs for each gene
        pident_plot(dbd_df)
        plt.savefig(snakemake.output["png"])
        plt.close()

        input_dir = Path(snakemake.input["motif_dir"])
        output_dir = Path(snakemake.output["motif_dir"])
        if not output_dir.exists():
            output_dir.mkdir()
        pident_threshold = snakemake.params["pident"]
        motif_table = Path(snakemake.output["motif_table"])
        # set up motif table file
        if motif_table.exists():
            motif_table.unlink()
        selected = dbd_df.query("pident >= @pident_threshold")
        missed_TFs = set(dbd_df.TF) - set(selected.TF)

        for idx in selected.index:
            motif = selected.loc[idx, "motif_id"]
            gene = selected.loc[idx, "TF"].replace("Lv-", "Sp-")
            with open(input_dir.joinpath(motif + ".txt"), "r") as handle:
                try:
                    pwm = motifs.read(handle, "pfm-four-columns")
                except:
                    logging.warning(
                        "Unable to open file: %s", input_dir.joinpath(motif + ".txt")
                    )
                pwm.id = pwm.consensus
                pwm.name = pwm.consensus
            # have to scale counts bc otherwise all 0 in clusterbuster format
            scale_counts(pwm, 100)
            with open(output_dir.joinpath(motif + ".cb"), "w") as handle:
                handle.write(pwm.format("clusterbuster"))
            with open(motif_table, "a") as handle:
                handle.write(
                    ",".join([motif, "_".join(motif.split("_")[2:]), gene]) + "\n"
                )
                # handle.write(get_line(motif, gene, pwm.consensus))
