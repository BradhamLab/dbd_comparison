import pandas as pd

if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        results = pd.read_csv(snakemake.input["blast"])
        anno_sep = ","
        if ".txt" in snakemake.input["anno"]:
            anno_sep = "\t"
        annos = pd.read_csv(snakemake.input["anno"], sep=anno_sep, index_col=0)
        selected = results.query("evalue < 0.05 and pident > 50").copy()
        selected["gene_id"] = selected["qseqid"].apply(lambda x: x.split(":")[0])
        selected["matches_annotation"] = (
            annos.loc[selected.gene_id, "uniprot.hit"].values
            == selected.sseqid.apply(lambda x: x.split(".")[0]).values
        )
        selected["passes_pident"] = selected.pident > 70
        keep = selected.query("matches_annotation and passes_pident")
        missed = selected[~selected.gene_id.isin(keep.gene_id)].gene_id.unique()
        new_rows = []
        for gene in missed:
            matches = selected.query("gene_id == @gene").sort_values(
                ["pident", "matches_annotation"], ascending=False
            )
            # alignment to annotation has highest pident, keep
            if matches.iloc[0, -2]:
                new_rows.append(matches.iloc[0, :])
                continue
            # keep highest pident entry
            new_rows.append(matches.iloc[0, :])
            # has *an* alignment that matches annotation, keep highest pident
            if matches.matches_annotation.sum() > 0:
                new_rows.append(matches.query("matches_annotation").iloc[0, :])
        combined = pd.DataFrame(new_rows)
        pd.concat(
            [
                x.sort_values(
                    ["gene_id", "passes_pident", "matches_annotation"],
                    ascending=[True, False, False],
                )
                for x in [keep, combined]
            ]
        ).reset_index(drop=True)[
            [
                "gene_id",
                "qseqid",
                "passes_pident",
                "matches_annotation",
                "pident",
                "length",
                "sseqid"
            ]
        ].replace(
            {'gene_id': annos.loc[selected.gene_id, 'SPU.common_name']}
        ).rename(
            columns={'gene_id': 'gene_name', 'qseqid': 'domain', 'sseqid': "anno"}
        ).to_csv(
            snakemake.output['csv'],
            index=False
        )
