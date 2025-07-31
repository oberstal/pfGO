# pfGO 2.2.0

Updates include:

* New helper-functions `reformat_sigGenes()` and `make_enrichRes()` to convert pfGO output to format easily visualized via [enrichplot package functions](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html)
    * See `make_enrichRes()` documentation for pfGO output-to-plot example workflow.\*
* Updated *P. falciparum* GO databases and gene annotations are included (made from PlasmoDB release 68; accessed July 29, 2025):
    * Data objects `Pfal_geneID2GO`, `Pfal_geneID2GO_curated`, and `pf.annot`


\*Ideally would add a full working example of taking pfGO output and turning it into various plots to the documentation for major release. For now I have example code in the `make_enrichRes()` documentation that isn't run/hasn't been well tested. 

# pfGO 2.1

Minor update to remove nonfunctional 'remote' url for dependency that caused some installs to fail

# pfGO 2.0

This release of pfGO includes several updates to increase utility of output files/remove the need for manual parsing before downstream analyses. Several functions are updated to make automated retrieval/incorporation of latest PlasmoDB annotations for gene products and GO terms easy.

* Output files for significant genes in significant terms now include mapped GO terms, their definitions, enrichment values and interest-category
* Formatting fixes to all significant-genes-in-significant-terms output-files for easier downstream parsing.
* Functions pulling data from PlasmoDB.org now require manual URL entry to circumvent inconsistencies in link-structure


# pfGO 1.2

* Formatting fixes to all significant-genes-in-significant-terms output-files for easier downstream parsing.
* Functions pulling data from PlasmoDB.org updated to 'current' release
    * no longer hard-coded to specific version release--should work to access all "Current_Release" links so long as there are no updates to PlasmoDB url structure or file-compression (works with non-compressed files and gzip-compressed files)
    * get.annot function now automatically pulls all 'gene' records from .gff files ("protein_coding_gene","ncRNA_gene","pseudogene")
* Included example-datasets updated from latest PlasmoDB release (version 66)


# pfGO 1.0

* Added a `NEWS.md` file to track changes to the package.
