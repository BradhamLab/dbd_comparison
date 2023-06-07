import logging
from pathlib import Path
import itertools
import shutil

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
    with sns.plotting_context(
        "paper",
        rc={
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.labelsize": 12,
            "xtick.majorsize": 1,
            "ytick.majorsize": 1,
            "font.sans-serif": "inter",
            "savefig.dpi": 300,
        },
    ):
        fig = plt.figure()
        # dbd_df.TF = dbd_df.TF.str.replace("Lv-", "")
        plot = sns.displot(
            data=dbd_df.sort_values("TF"),
            x="pident",
            hue="TF",
            col="TF",
            col_wrap=5,
            facet_kws={"sharey": False, "sharex": True},
            legend=False,
            bins=12,
            height=1.5,
        )
        for ax in plot.axes:
            ax.axvline(x=threshold, ls="--", c="black", lw=1)
            ax.yaxis.set_major_formatter(FormatStrFormatter("%d"))
            ax.set_title(ax.get_title().replace("TF = ", ""), fontsize=10)
        fig.subplots_adjust(hspace=0.02, wspace=0.02)


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


def plot_motif_scores(*args, **kwargs):
    data = kwargs.pop("data")
    plot_df = data.assign(
        Rank=lambda df_: df_.score.rank(ascending=True), is_min=False
    ).rename(columns={"score": "I.C."})
    # plot_df.loc[plot_df["I.C."].idxmin(), "is_min"] = True
    ax = plt.gca()
    ax.set(aspect=1)
    ax.axis("off")
    scatter_ax = ax.inset_axes([0, 0, 0.8, 1])
    hist_ax = ax.inset_axes([0.85, 0, 0.2, 1])
    for i in range(2):
        sns.scatterplot(
            data=plot_df[plot_df.selected == i],
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
    sns.histplot(
        data=plot_df.reset_index(), y="I.C.", ax=hist_ax, bins=bins, color="#7f7f7f"
    )
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
            # "xtick.labelsize": 12,
            # "ytick.labelsize": 12,
            # "axes.labelsize": 14,
            # "xtick.majorsize": 1,
            # "ytick.majorsize": 1,
            "font.sans-serif": "inter",
            "savefig.dpi": 300,
        },
    ):
        fig = plt.figure()
        g = sns.FacetGrid(df.sort_values("TF"), col="TF", col_wrap=4)
        g.map_dataframe(plot_motif_scores)
        fig.subplots_adjust(hspace=0.05, wspace=0.05)


def close_plot():
    plt.cla()
    plt.clf()
    plt.close()


def read_pwm(pwm_in):
    """Read PWM file."""
    with open(pwm_in, "r") as pwm_handle:
        try:
            return motifs.read(pwm_handle, "pfm-four-columns")
        except:
            logging.warning("Unable to open file: %s", pwm_in)
            return None


def score_motif(motif_fn, background={"G": 0.35, "A": 0.15, "T": 0.15, "C": 0.35}):
    motif_pwm = read_pwm(motif_fn)
    if motif_pwm is None:
        return np.nan
    motif_pwm.background = background
    motif_pwm.pseudocounts = {k: 0.25 * v for k, v in background.items()}
    info = 0
    for nt, pos in itertools.product("GATC", range(motif_pwm.length)):
        info += (
            -1
            * motif_pwm.pwm[nt][pos]
            * np.log(motif_pwm.pwm[nt][pos] / motif_pwm.background[nt])
        )
    return info


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

        with open(motif_table, "w") as handle:
            handle.write("matrix_id,base_id,gene" + "\n")
        dbd_df.set_index("motif_id", inplace=True)
        dbd_df["score"] = dbd_df.apply(
            lambda x: score_motif(str(input_dir.joinpath(x.name + ".txt"))), axis=1
        )
        dbd_df = dbd_df[~dbd_df.score.isna()]

        purps = dbd_df.assign(
            is_purp=lambda df_: df_.index.str.match("Strongylocentrotus_purpuratus")
        ).query("is_purp")

        selected = (
            (
                dbd_df.query("pident >= @pident_threshold").sort_values(
                    ["pident", "score"], ascending=[True, False]
                )
                # .groupby("query_id")
                # .tail(10)
                # .sort_values("score", ascending=False)
                # .groupby("query_id")
                # .tail(1)
            )
            .reset_index()
            .set_index("motif_id")
        )
        selected.to_csv(snakemake.output["score_df"])
        # denote selected motifs, plot
        dbd_df["selected"] = dbd_df.index.isin(selected.index)
        score_plot(dbd_df.query("pident >= @pident_threshold"))
        plt.savefig(snakemake.output["score_plot"])
        close_plot()

        retained_motifs = list(set(purps.index.values).union(selected.index.values))
        missed_TFs = set(dbd_df.TF) - set(dbd_df.loc[retained_motifs, "TF"])
        logging.warning("No motif selected for %s", ", ".join(missed_TFs))

        # convert normalized PWMs to clusterbuster format for retained motifs
        # log ids + gene names
        rows = []
        for motif in retained_motifs:
            gene = dbd_df.loc[motif, "TF"]
            if not isinstance(gene, str):
                gene = gene.values[0]
            pwm = read_pwm(input_dir.joinpath(motif + ".txt"))
            if pwm is None:
                continue
            pwm.id = pwm.consensus
            pwm.name = pwm.consensus
            scale_counts(pwm, 100)
            with open(output_dir.joinpath(motif + ".cb"), "w") as pwm_out:
                pwm_out.write(pwm.format("clusterbuster"))
            rows.append(
                [motif, "_".join(motif.split("_")[2:]), gene.replace("Lv-", "Sp-")]
            )

        # join retained motifs with extra motifs
        pd.concat(
            [
                pd.DataFrame(rows, columns=["matrix_id", "base_id", "gene"]),
                pd.read_csv(snakemake.input["extra_motifs_csv"]),
            ]
        ).to_csv(motif_table, index=None)

        for each in Path(snakemake.input["extra_motifs"]).glob("*.cb"):
            shutil.copy(each, output_dir.joinpath(each.stem + ".cb"))
