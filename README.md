# dbd_comparison
Repository for comparing DNA Binding Sites for motif selection

## Workflow

This workflow looks to compare DNA binding domains between a query species of specified proteins (e.g. PMC TFs in _lytechinus variegatus_) to possible orthologues in other species in order to determine candidate binding motifs for the query species. The general workflow follows steps outline here:

1. Extract protein amino acids for TFs of interest
2. Annotate each protein sequenc for PFAM domains using HMMER
3. Filter HMMER output to remove overlapped/redundant protein domains, domains with high e-values, and non-DBD domains
4. Match query sequences to UNIPROT orthologues using genome annotations.
5. Download sequences for all orthologues using BioMaRT and ENSEMBL.
6. Repeat 1-3 for orthologous sequences in other species
7. Compare DBD sequences between orthologues in query species (lv) and all other target species (ensembl sequences)
8. Match TF binding motifs to query proteins if a DBD from the given protein has > 70 percent identity with the given orthologous protein.



![](dag.png)


## Installation

To install the pipeline, simply clone this repository. The pipeline requires `snakemake` `hmmer`, `muscle`, `seqtk`, and `pandas` and `Biopython` Python libraries. The current implementation was created for BU's shared computing clusters, and automatically loads `hmmer`, `muscle`, and `seqtk` modules. It also requires the [pfam extraction workflow](https://github.com/BradhamLab/extract_pfam_domains/)

## Configuration

To modify the domains and sequences of interest, modify the [configuration file](config.yaml) as appropriate.
