import logging
from pathlib import Path
import itertools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import motifs
from matplotlib.ticker import FormatStrFormatter


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


def score_motif(motif, background={"G": 0.35, "A": 0.15, "T": 0.15, "C": 0.35}):
    motif.background = background
    motif.pseudocounts = {k: 0.25 * v for k, v in background.items()}
    info = 0
    for nt, pos in itertools.product("GATC", range(motif.length)):
        info += (
            -1 * motif.pwm[nt][pos] * np.log(motif.pwm[nt][pos] / motif.background[nt])
        )
    return info


def plot_motif_scores(*args, **kwargs):
    data = kwargs.pop("data")
    plot_df = data.assign(
        Rank=lambda df_: df_.score.rank(ascending=True), is_min=False
    ).rename(columns={"score": "I.C."})
    plot_df.loc[plot_df["I.C."].idxmin(), "is_min"] = True
    ax = plt.gca()
    ax.set(aspect=1)
    ax.axis("off")
    scatter_ax = ax.inset_axes([0, 0, 0.8, 1])
    hist_ax = ax.inset_axes([0.85, 0, 0.2, 1])
    for i in range(2):
        sns.scatterplot(
            data=plot_df[plot_df.is_min == i],
            x="Rank",
            y="I.C.",
            hue="is_min",
            palette={True: "#d62728", False: "#7f7f7f"},
            hue_order=[True, False],
            legend=False,
            ax=scatter_ax,
            zorder=i,
        )
    bins = "auto"
    if plot_df["I.C."].nunique() < 5:
        bins = 5
    sns.histplot(data=plot_df, y="I.C.", ax=hist_ax, bins=bins, color="#7f7f7f")
    for i, p in enumerate(hist_ax.patches):
        if i == 0:
            p.set_facecolor("#d62728")
            p.set_zorder(1)
        else:
            p.set_zorder(0)
    hist_ax.set_ylabel("")
    hist_ax.set_yticks([])

    for spine in ["top", "right"]:
        scatter_ax.spines[spine].set_visible(False)
        hist_ax.spines[spine].set_visible(False)
    # hist_ax.spines["bottom"].set_visible(False)
    for each in [scatter_ax, hist_ax]:
        each.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))


def score_plot(df):
    with sns.plotting_context(
        "talk",
        rc={
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "axes.labelsize": 14,
            "xtick.majorsize": 1,
            "ytick.majorsize": 1,
        },
    ):
        fig = plt.figure(figsize=(6, 8))
        g = sns.FacetGrid(df.sort_values("TF"), col="TF", col_wrap=4)
        g.map_dataframe(plot_motif_scores)
        fig.subplots_adjust(hspace=0.05, wspace=0.05)


def close_plot():
    plt.cla()
    plt.clf()
    plt.close()


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
        plt.savefig(snakemake.output["pident_plot"])
        close_plot()

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
        motif_scores = {"motif": [], "score": [], "TF": []}
        best_motifs = {}
        for idx in selected.index:
            motif = selected.loc[idx, "motif_id"]
            gene = selected.loc[idx, "TF"]
            with open(input_dir.joinpath(motif + ".txt"), "r") as handle:
                try:
                    pwm = motifs.read(handle, "pfm-four-columns")
                except:
                    logging.warning(
                        "Unable to open file: %s", input_dir.joinpath(motif + ".txt")
                    )
                    continue
                pwm.id = pwm.consensus
                pwm.name = pwm.consensus
                pwm_info = score_motif(pwm)
                motif_scores["motif"].append(motif)
                motif_scores["score"].append(pwm_info)
                motif_scores["TF"].append(gene)
                if gene not in motif_scores:
                    best_motifs[gene] = {"motif": pwm, "info": pwm_info, "name": motif}
                if pwm_info < best_motifs[gene]["info"]:
                    best_motifs[gene] = {"motif": pwm, "info": pwm_info, "name": motif}

        # plot motif scores
        score_df = pd.DataFrame(motif_scores)
        score_plot(score_df)
        plt.savefig(snakemake.output["score_plot"])
        close_plot()

        for tf, each in best_motifs.items():
            # have to scale counts bc otherwise all 0 in clusterbuster format
            scale_counts(each["motif"], 100)
            with open(output_dir.joinpath(motif + ".cb"), "w") as handle:
                handle.write(each["motif"].format("clusterbuster"))
            with open(motif_table, "a") as handle:
                handle.write(
                    ",".join(
                        [
                            each["name"],
                            "_".join(each["name"].split("_")[2:]),
                            tf.replace("Lv-", "Sp-"),
                        ]
                    )
                    + "\n"
                )
                # handle.write(get_line(motif, gene, pwm.consensus))
