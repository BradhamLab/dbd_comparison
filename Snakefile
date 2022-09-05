configfile: "config.yaml"


shell.prefix(
    "module load hmmer; module load muscle; module load seqtk; module load blast+; "
)


def yaml_to_config_dict(yaml_file):
    import yaml

    with open(yaml_file, "r") as f:
        return yaml.safe_load(f)


rule all:
    input:
        selected_motifs=directory("output/retained_motifs"),
        matched="output/matched_sequences.csv",
        missed="output/missed_sequences.csv",
        # blast="output/retained_blast_out.csv"


rule extract_protein_fasta:
    input:
        ids=config["tf_ids"],
        fasta=config["protein_fasta"],
    output:
        "output/tf_proteins.fasta",
    shell:
        "seqtk subseq {input.fasta} {input.ids} > {output}"


rule find_uniprot_orthologues:
    input:
        ids="data/grn_genes.txt",
        fasta=config["protein_fasta"],
    output:
        fasta=temp("output/genes.fasta"),
        ortho="output/orthologue_blast.csv",
    params:
        db=config["blast_db"],
        columns=",".join(
            [
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
            ]
        ),
    shell:
        "seqtk subseq {input.fasta} {input.ids} > {output.fasta}; "
        "echo '{params.columns}' > {output.ortho}; "
        "blastp -query {output.fasta} -db {params.db} -outfmt 6 | "
        "sed 's/\\t/,/g' >> {output.ortho}"


rule run_hmmsearch:
    input:
        hmm=config["hmm"],
        fasta="output/tf_proteins.fasta",
    output:
        "output/tf_domains.domtblout",
    shell:
        "hmmsearch --domtblout {output} {input.hmm} {input.fasta}"


rule format_domtblout_pfam:
    input:
        "output/tf_domains.domtblout",
    output:
        "output/pfam_domains.csv",
    params:
        columns=",".join(
            [
                "target_name",
                "tlen",
                "query_name",
                "qlen",
                "evalue",
                "hmm_from",
                "hmm_to",
                "seq_from",
                "seq_to",
            ]
        ),
    shell:
        "echo '{params.columns}' > {output}; "
        "sh scripts/parse_hmmsearch.sh {input} >> {output}"


rule matched_sequences:
    input:
        "output/pfam_domains.csv",
    output:
        "output/matched_sequences.csv",
    shell:
        "cat {input} | cut -d',' -f1 | tail -n +2  | sort | uniq > {output}"


rule missed_sequences:
    input:
        f1=config["tf_ids"],
        f2="output/matched_sequences.csv",
    output:
        "output/missed_sequences.csv",
    shell:
        "sort {input.f1} | comm -23 - {input.f2} > {output}"


rule extract_and_filter_domains:
    input:
        sequence_domains="output/pfam_domains.csv",
        domain_list_files=config["allowed_domains"]["files"],
    params:
        domain_regex=config["allowed_domains"]["domain_regex"],
        e_val=config["eval"],
    output:
        csv="output/retained_pfam_domains.csv",
        filtered="output/filtered_domains.csv",
        fasta="output/retained_domains.fasta",
    script:
        "scripts/filter_and_merge_domains.py"


rule blast_selected_domains:
    input:
        "output/retained_domains.fasta",
    output:
        "output/blast_out.csv",
    params:
        db=config["blast_db"],
        columns=",".join(
            [
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
            ]
        ),
    shell:
        "echo '{params.columns}' > {output}; "
        "blastp -task blastp-short -query {input} -db {params.db} -outfmt 6 | "
        "sed 's/\\t/,/g' >> {output}"


rule filter_blast_output:
    input:
        blast="output/blast_out.csv",
        anno=config["annos"],
    output:
        csv="output/retained_blast_out.csv",
    script:
        "scripts/filter_blast.py"


from pathlib import Path


rule download_ensembl_sequences:
    input:
        motif_table="data/motif_pwm_df.csv",
    log:
        "logs/ensembl_download.log"
    output:
        directory("output/ensembl_sequences/"),
    script:
        "scripts/get_protein_sequences.R"


ENSEMBL_OUT = Path(config["output"]["ensembl_dir"])


rule combine_ensembl_fasta:
    input:
        ensembl="output/ensembl_sequences/",
        external="data/external_sequences/",
    output:
        ENSEMBL_OUT.joinpath("ensembl_proteins.fa"),
    shell:
        "cat {input.ensembl}/*.fa > {output} "
        "cat {input.external}/*.fa >> {output}"


rule ensembl_sequence_ids:
    input:
        ENSEMBL_OUT.joinpath("ensembl_proteins.fa"),
    output:
        ENSEMBL_OUT.joinpath("protein_ids.txt"),
    shell:
        "grep '>' {input} | sed 's/>//g' | sort | uniq > {output}"


ensembl_config = yaml_to_config_dict(
    config["workflows"]["domain_detection"]["ensembl_config"]
)


module domain_detection:
    snakefile:
        "/projectnb/bradham/workflows/domain_identification/Snakefile"
    config:
        ensembl_config


# run_hmmsearch, input protein to output/ensembl_proteinsfa
use rule * from domain_detection as ensembl_*


use rule run_hmmsearch from domain_detection as ensembl_run_hmmsearch with:
    input:
        hmm=config["hmm"],
        fasta=ENSEMBL_OUT.joinpath("ensembl_proteins.fa"),


use rule missed_sequences from domain_detection as ensembl_missed_sequences with:
    input:
        f1=ENSEMBL_OUT.joinpath("protein_ids.txt"),
        f2=rules.ensembl_missed_sequences.input["f2"],


use rule extract_and_filter_domains from domain_detection as ensembl_extract_and_filter_domains with:
    input:
        tf_fasta=ENSEMBL_OUT.joinpath("ensembl_proteins.fa"),
        sequence_domains=Path(ensembl_config["outdir"]).joinpath("pfam_domains.csv"),
        domain_list_files=ensembl_config["allowed_domains"]["files"],


rule compare_dbds:
    input:
        ensembl_fa=rules.extract_and_filter_domains.output["fasta"],
        tf_fa="output/retained_domains.fasta",
        motif_table="data/motif_pwm_df.csv",
        id2name="all_tf_models.csv",
    output:
        csv="output/dbd_alignment.csv",
    log:
        "logs/dbd_alignment.log",
    script:
        "scripts/compare_dbds.py"


rule select_motifs:
    input:
        id2name="all_tf_models.csv",
        csv="output/dbd_alignment.csv",
        motif_dir="data/cisbg_motifs",
    log:
        "logs/motif_selection.log",
    params:
        pident=0.7,
    output:
        motif_dir=directory("output/retained_motifs"),
        motif_table="output/motif_table.tsv",
        png="output/pident_distributions.png",
    script:
        "scripts/format_and_filter_motifs.py"
