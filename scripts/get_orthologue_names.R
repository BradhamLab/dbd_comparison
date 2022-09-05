library(biomaRt)
library(dplyr)
library(stringr)

if exists("snakemake") {
    blast_output <- read.csv("output/orthologue_blast.csv")
    scaffold_to_gene_name <- read.csv('gene_model_by_count.csv')
    row.names(scaffold_to_gene_name) <- scaffold_to_gene_name$L_var_ID
    blast_output['urchin_name'] = scaffold_to_gene_name[blast_output$qseqid, 'name']
    ensembl <- useMart('ensembl', dataset="hsapiens_gene_ensembl")
    annotation <- getBM(
        attributes=c("uniprotswissprot", "hgnc_symbol", "uniprot_gn_symbol"),
        filters="uniprotswissprot",
        values=gsub("\\..*$", "", blast_output$sseqid),
        mart=ensembl
    )
    to_hgnc <- annotation$hgnc_symbol
    names(to_hgnc) <- annotation$uniprotswissprot
    best_matches <- blast_output %>%
        dplyr::filter(evalue < 0.001) %>%
        dplyr::mutate(uniprot = stringr::str_replace(sseqid, "\\..*$", "")) %>%
        dplyr::mutate(gene_name = to_hgnc[uniprot]) %>%
        dplyr::filter(!is.na(gene_name)) %>%
        dplyr::group_by(urchin_name) %>%
        dplyr::slice_min(order_by = evalue, n = 5, with_ties=TRUE) %>%
        dplyr::slice_max(order_by = pident, n=1, with_ties=TRUE) %>%
        dplyr::select(
            c('urchin_name', 'gene_name', 'qseqid', 'sseqid', 'pident', 'evalue')
        ) %>%
        dplyr::arrange(-pident, urchin_name)

    missed <- scaffold_to_gene_name %>%
        dplyr::filter(!name %in% best_matches$urchin_name) %>%
        dplyr::group_by(name) %>%
        dplyr::slice_max(order_by = name, n=1, with_ties=FALSE) %>%
        dplyr::select(
            c('name', 'L_var_ID')
        ) %>% dplyr::inner_join(
            blast_output %>%
                dplyr::filter(qseqid %in% missed$L_var_ID) %>%
                dplyr::group_by(qseqid) %>%
                dplyr::slice_min(order_by=evalue, n=5, with_ties=TRUE) %>%
                dplyr::slice_max(order_by=pident, n=1, with_ties=FALSE) %>%
                dplyr::mutate(uniprot = stringr::str_replace(sseqid, "\\..*$", "")) %>%
                dplyr::mutate(hgnc = to_hgnc[uniprot]),
            by=c('L_var_ID' = 'qseqid'),
            keep=FALSE
        ) %>%
        dplyr::select(c('name', 'L_var_ID', 'uniprot', 'hgnc', 'evalue', 'pident'))

    write.csv(
        best_matches,
        "output/human_orthologues.csv",
        quote=FALSE,
        row.names=FALSE
    )

    write.csv(
        missed,
        "output/no_orthologues.csv",
        quote=FALSE,
        row.names=FALSE
    )
}
