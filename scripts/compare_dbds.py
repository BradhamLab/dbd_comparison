import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices
from tqdm import tqdm
import logging


blosum62 = substitution_matrices.load("BLOSUM62")


def get_gene_name(seq_id, lookup=None):
    if lookup is None:
        return seq_id.split(":")[0]
    return lookup[seq_id.split(":")[0]]


def p_identity(alignemnt):
    n = len(alignemnt.seqA)
    return sum([alignemnt.seqA[i] == alignemnt.seqB[i] for i in range(n)]) / n


if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        motif_df = pd.read_csv(snakemake.input["motif_table"])
        reference_dbds = SeqIO.parse(snakemake.input["ensembl_fa"], "fasta")
        records = SeqIO.parse(snakemake.input["tf_fa"], "fasta")
        lookup = (
            pd.read_csv(snakemake.input["id2name"])
            .set_index("gene_id")
            .to_dict()["name"]
        )
        logging.basicConfig(filename=snakemake.log[0], level=logging.WARNING)
        query_dbds = {}
        for each in records:
            gene = get_gene_name(each.id, lookup)
            if gene not in query_dbds:
                query_dbds[gene] = [each]
            else:
                query_dbds[gene].append(each)
        reference_ids = []
        motif_ids = []
        query_ids = []
        percent_identities = []
        missed_genes = []
        for dbd in tqdm(reference_dbds):
            ensembl_id = get_gene_name(dbd.id)
            possible_motifs = motif_df.query("gene_id == @ensembl_id")
            related_motifs = possible_motifs.motif.values.tolist()
            unique_genes = possible_motifs.gene.unique()
            for ref_gene in unique_genes:
                try:
                    to_query = query_dbds[ref_gene]
                except KeyError:
                    to_query = []
                    missed_genes.append(ref_gene)

                for query_seq in to_query:
                    # default blast parameters
                    try:
                        alignment = pairwise2.align.globalds(
                            dbd.seq, query_seq.seq, blosum62, -11, -1
                        )
                    except SystemError:
                        continue
                    if len(alignment) > 0:
                        reference_ids += [dbd.id] * len(related_motifs)
                        motif_ids += related_motifs
                        query_ids += [query_seq.id] * len(related_motifs)
                        percent_identities += [p_identity(alignment[0])] * len(
                            related_motifs
                        )

        pd.DataFrame(
            {
                "ensembl_dbd": reference_ids,
                "motif_id": motif_ids,
                "query_dbd": query_ids,
                "pident": percent_identities,
            }
        ).assign(
            gene_id=lambda df_: [x.split(":")[0] for x in df_.ensembl_dbd],
            query_id=lambda df_: [x.split(":")[0] for x in df_.query_dbd],
        ).to_csv(
            snakemake.output["csv"]
        )
        logging.warning("Missed Sequences: %s", ", ".join(set(missed_genes)))
