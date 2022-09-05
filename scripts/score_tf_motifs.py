from multiprocessing.sharedctypes import Value
import pandas as pd
import pathlib
from Bio import motifs, SeqIO
import itertools

def calc_pseudocount(motif, p = 0.001):
    top_count = max([max(c) for c in motif.counts.values()])
    return top_count * p

if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:

        domain_df = pd.read_csv(snakemake.input['domains']).assign(
            nt_start=lambda x: x.seq_from * 3,
            nt_end=lambda x: x.nt_start + (x.seq_to - x.seq_from) * 3,
        )
        domain_df['index'] = domain_df.apply(
            lambda x: f"{x.gene}:{x.nt_start}-{x.nt_end}",
            axis=1
        )
        domain_df.set_index('index', inplace=True)
        motif_dir = pathlib.Path(snakemake.input['motif_dir'])

        tf_motifs = {}
        ext_to_format = {'.cb': 'clusterbuster', '.jaspar': 'jaspar'}
        for each in motif_dir.glob("*"):
            with open(each) as handle:
                motif = motifs.read(handle, ext_to_format[each.suffix])
                pseudo = calc_pseudocount(motif, 10**(-3))
                motif.pseudocounts = {x: pseudo for x in 'GATC'}
                tf_motifs[each.stem] = {
                    "motif": motif,
                    "threshold": motif.pssm.distribution().threshold_fpr(10**(-6)),
                }

        transcripts = {x.id: x for x in SeqIO.parse(snakemake.input['transcripts'], "fasta")}
        no_binding = set(transcripts.keys()) - set(domain_df.gene)
        motif_df = pd.DataFrame(
            index=pd.MultiIndex.from_tuples(
                itertools.product(domain_df.index.values, sorted(tf_motifs.keys()))
            ),
            columns=["max_score", "passed"],
        )
        motif_df.index = motif_df.index.rename(['gene_loc', 'motif'])
        # motif_df['gene'] = motif_df.apply(lambda x: x.split(':')[0], axis=1)
        motif_df.max_score = pd.to_numeric(motif_df.max_score)
        motif_df.passed = ~motif_df.passed.astype(bool)

        for (idx, gene) in zip(domain_df.index, domain_df.gene.values):
            dbd = transcripts[gene][slice(*domain_df.loc[idx, ["nt_start", "nt_end"]])]
            for name, mdict in tf_motifs.items():
                score = max(
                    mdict["motif"].pssm.calculate(dbd.seq).max(),
                    mdict["motif"].pssm.reverse_complement().calculate(dbd.seq).max(),
                )
                motif_df.loc[(idx, name), ["max_score", "passed"]] = [
                    score,
                    score > mdict["threshold"],
                ]
        motif_df.reset_index().to_feather(snakemake.output["feather"])
