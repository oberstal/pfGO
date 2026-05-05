
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pfGO

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10019755.svg)](https://doi.org/10.5281/zenodo.10019755)
<!-- badges: end -->

The goal of pfGO is to package handy functions and datasets specifically
enabling *Plasmodium falciparum* functional enrichment analyses. pfGO
acts as a wrapper around the topGO package for much of its enrichment
functionality, while also providing several functions for incorporating
latest gene-ontology and functional annotations.

pfGO enables chaining together many parallel enrichment-analyses at
once, generating thorough logs and outputs supporting reproducible
analyses.

## Installation

You can install pfGO from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("oberstal/pfGO")

# or if authenticating via ssh:
devtools::install_git("https://github.com/oberstal/pfGO")
```

See individual functions/data objects for further documentation. E.g.:

``` r
?run.topGO.meta
?Pfal_geneID2GO

# for all available functions/data:
?pfGO
# then scroll to the bottom and click the "index" link.
```

## Example (quick start)

This is a basic example demonstrating how to run an enrichment analysis
on *piggyBac* pooled phenotypic screening results to identify processes
enabling parasite survival of host fever (using data [as published
previously](https://doi.org/10.1038/s41467-021-24814-1)):

``` r
library(pfGO)
```

``` r
# load included pf GO database and example-data to be tested for functional enrichment
data(Pfal_geneID2GO)
data(exampleMydf)

# run the topGO pipeline on all experimental categories of interest from exampleMydf
run.topGO.meta(mydf = exampleMydf, geneID2GO = Pfal_geneID2GO, pval = 0.05)
```

View example console output [here](#example-console-output)

## Key functions and further documentation

### run.topGO.meta

#### description:

Tests for functional enrichment in gene-categories of interest.
Enrichments are performed by each ontology (molecular function,
biological process, cellular component), for each interest category,
sequentially.

#### parameters:

- `mydf`: data frame with geneIDs in column 1, and interest-category
  classifications in column 2. Additional columns are ignored.
- `geneID2GO`: a list of named vectors of GO IDs–one vector of GO-terms
  for each geneID.
- `pval`: p-value threshold for significance (Fisher test). Defaults to
  0.05.
- `minTermSize`: minimum number of genes that must be mapped to a GO
  term for it to be included in the enrichment analysis. Defaults to 5.
- `algorithm`: algorithm to run for the enrichment analysis. Accepted
  values are c(“classic”, “elim”, “weight”, “weight01”, “lea”,
  “parentchild”). Defaults to “weight01”.
  - See the [topGO package
    documentation](https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO_manual.html#1_Introduction)
    and associated publications for details on algorithms.

#### outputs:

run.topGO.meta creates several output-files, including:

- enrichment results

- significant genes per significant term

  - available in `./Routput/GO/all.combined.sig.genes.per.sig.terms.tsv`

- plots of the GO-term hierarchy relevant to the analysis

- thorough log-files for each gene-category of interest tested against
  the background of all other genes in the analysis

- system logs recording all package versions, etc.

Primary results from run.topGO.meta will be in
`./Routput/GO/all.combined.GO.results.tsv`. Note that run.topGO.meta
will automatically create the `./Routput` directory (and other required
output directories nested in `./Routput`) in your working directory for
you if it does not exist.

#### more function details

The **run.topGO.meta** function:

- defines which genes are “interesting” and which should be defined as
  background for each category specified in mydf,
- makes the GOdata object for topGO,
- tests each category of interest for enriched GO-terms against all the
  other genes included in mydf (the “gene universe”),
- and then outputs results to several tables (.tsv files that can be
  opened in Excel).

Enrichments are performed by each ontology (molecular function,
biological process, cellular component; MF, BP, and CC, respectively)
sequentially on all groups of interest. Results are combined in the
final output-table (`./Routput/GO/all.combined.GO.results.tsv`).

TopGO automatically accounts for genes that cannot be mapped to GO terms
(or are mapped to terms with \< `minTermSize` genes in the analysis)
with “feasible genes” indicated in the topGO.log files in the
`./Routput/GO/` folder.

**Concepts for common use-cases**:

*RNAseq*:

In an RNAseq analysis, common interest-categories might be
“upregulated”, “downregulated”, and “neutral” genes. The gene universe
would consist of all genes expressed above your threshold cutoffs (*not
necessarily all genes in the genome*).

*piggyBac screens*:

In pooled *piggyBac*-mutant screening, common categories might be
“sensitive”, “tolerant”, and “neutral”. The gene universe would consist
of all genes represented in your screened library of mutants (*again,
not all genes in the genome*).

See the included data object *exampleMydf* as an example.

**Using your own custom GO database**:

A correctly formatted geneID2GO object is included for *P. falciparum*
enrichment analyses (*Pfal_geneID2GO*). You may also provide your own,
so long as it is a named character-vector of GO-terms (each vector named
by geneID, with GO terms as each element).

You can use the included **formatGOdb.curated()** function to format a
custom GO database from curated GeneDB annotations for several non-model
organisms (or the **formatGOdb()** function to include all GO
annotations, if you aren’t picky about including automated electronic
annotations). If you’re studying a model organism, several annotations
are already available through the AnnotationDbi bioconductor package
that loads with topGO.

<!-- @seealso [topGO::topGO()] -->

<a id="example-console-output"></a>

## Example console output

Example console output generated running the quick-start example data
(*piggyBac* pooled phenotypic screening results to identify processes
enabling parasite survival of host fever ([similar to as published
previously](https://doi.org/10.1038/s41467-021-24814-1)):

``` r
# load included pf GO database and example-data to be tested for functional enrichment
data(Pfal_geneID2GO_curated)
data(exampleMydf)

# run the topGO pipeline on all experimental categories of interest from exampleMydf
run.topGO.meta(mydf = exampleMydf, geneID2GO = Pfal_geneID2GO_curated, pval = 0.1)
#> 
#> Building most specific GOs .....
#>  ( 283 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 554 GO terms and 715 relations. )
#> 
#> Annotating nodes ...............
#>  ( 334 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 80 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 9:   1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 8:   5 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 7:   9 nodes to be scored    (11 eliminated genes)
#> 
#>   Level 6:   12 nodes to be scored   (33 eliminated genes)
#> 
#>   Level 5:   14 nodes to be scored   (66 eliminated genes)
#> 
#>   Level 4:   15 nodes to be scored   (94 eliminated genes)
#> 
#>   Level 3:   17 nodes to be scored   (127 eliminated genes)
#> 
#>   Level 2:   6 nodes to be scored    (218 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (314 eliminated genes)
#> Loading required package: Rgraphviz
#> Loading required package: grid
#> 
#> Attaching package: 'grid'
#> The following object is masked from 'package:topGO':
#> 
#>     depth
#> 
#> Attaching package: 'Rgraphviz'
#> The following objects are masked from 'package:IRanges':
#> 
#>     from, to
#> The following objects are masked from 'package:S4Vectors':
#> 
#>     from, to
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.A.HS.Sensitive.MF_weight01_10_all  --- no of nodes:  28 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: A.HS.Sensitive
#> Ontology: MF        GO.ID                                        Term Annotated Significant
#> 1  GO:0120543 macromolecular conformation isomerase ac...        14           3
#> 2  GO:0016740                        transferase activity        71          12
#> 3  GO:0003723                                 RNA binding        58          12
#> 4  GO:0003924                             GTPase activity         5           2
#> 5  GO:0005319                  lipid transporter activity         5           2
#> 6  GO:0016301                             kinase activity        26           6
#> 7  GO:0016791                        phosphatase activity        11           3
#> 8  GO:0004518                           nuclease activity         6           2
#> 9  GO:0004674    protein serine/threonine kinase activity        17           4
#> 10 GO:0005198                structural molecule activity        16           4
#>    Expected topGO go.category interest.category algorithm statistic
#> 1      2.14 0.064          MF    A.HS.Sensitive  weight01    fisher
#> 2     10.84 0.076          MF    A.HS.Sensitive  weight01    fisher
#> 3      8.86 0.158          MF    A.HS.Sensitive  weight01    fisher
#> 4      0.76 0.169          MF    A.HS.Sensitive  weight01    fisher
#> 5      0.76 0.169          MF    A.HS.Sensitive  weight01    fisher
#> 6      3.97 0.222          MF    A.HS.Sensitive  weight01    fisher
#> 7      1.68 0.227          MF    A.HS.Sensitive  weight01    fisher
#> 8      0.92 0.229          MF    A.HS.Sensitive  weight01    fisher
#> 9      2.60 0.251          MF    A.HS.Sensitive  weight01    fisher
#> 10     2.44 0.278          MF    A.HS.Sensitive  weight01    fisher
#> 
#> ==============================================================================
#> 
#> interest-category 1 of 4
#> ontology 1 of 3
#> 
#> Building most specific GOs .....
#>  ( 358 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 1125 GO terms and 2077 relations. )
#> 
#> Annotating nodes ...............
#>  ( 354 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 187 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 13:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 12:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 11:  3 nodes to be scored    (11 eliminated genes)
#> 
#>   Level 10:  7 nodes to be scored    (11 eliminated genes)
#> 
#>   Level 9:   12 nodes to be scored   (24 eliminated genes)
#> 
#>   Level 8:   18 nodes to be scored   (51 eliminated genes)
#> 
#>   Level 7:   22 nodes to be scored   (87 eliminated genes)
#> 
#>   Level 6:   33 nodes to be scored   (127 eliminated genes)
#> 
#>   Level 5:   40 nodes to be scored   (200 eliminated genes)
#> 
#>   Level 4:   27 nodes to be scored   (256 eliminated genes)
#> 
#>   Level 3:   17 nodes to be scored   (327 eliminated genes)
#> 
#>   Level 2:   5 nodes to be scored    (343 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (348 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.A.HS.Sensitive.BP_weight01_10_all  --- no of nodes:  57 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: A.HS.Sensitive
#> Ontology: BP        GO.ID                                Term Annotated Significant
#> 1  GO:0010468       regulation of gene expression        27           6
#> 2  GO:0044283 small molecule biosynthetic process         5           3
#> 3  GO:0061024               membrane organization         9           4
#> 4  GO:0035556   intracellular signal transduction        22           7
#> 5  GO:0009410     response to xenobiotic stimulus        28           8
#> 6  GO:0008610          lipid biosynthetic process        12           3
#> 7  GO:0006364                     rRNA processing        14           6
#> 8  GO:0006401               RNA catabolic process         7           3
#> 9  GO:0000460             maturation of 5.8S rRNA         7           3
#> 10 GO:0006163 purine nucleotide metabolic process         8           3
#>    Expected  topGO go.category interest.category algorithm statistic
#> 1      4.27 0.0073          BP    A.HS.Sensitive  weight01    fisher
#> 2      0.79 0.0298          BP    A.HS.Sensitive  weight01    fisher
#> 3      1.42 0.0384          BP    A.HS.Sensitive  weight01    fisher
#> 4      3.48 0.0421          BP    A.HS.Sensitive  weight01    fisher
#> 5      4.43 0.0556          BP    A.HS.Sensitive  weight01    fisher
#> 6      1.90 0.0673          BP    A.HS.Sensitive  weight01    fisher
#> 7      2.21 0.0752          BP    A.HS.Sensitive  weight01    fisher
#> 8      1.11 0.0822          BP    A.HS.Sensitive  weight01    fisher
#> 9      1.11 0.0822          BP    A.HS.Sensitive  weight01    fisher
#> 10     1.27 0.1171          BP    A.HS.Sensitive  weight01    fisher
#> 
#> ==============================================================================
#> 
#> interest-category 1 of 4
#> ontology 2 of 3
#> 
#> Building most specific GOs .....
#>  ( 195 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 367 GO terms and 616 relations. )
#> 
#> Annotating nodes ...............
#>  ( 526 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 99 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   6 nodes to be scored    (16 eliminated genes)
#> 
#>   Level 8:   13 nodes to be scored   (65 eliminated genes)
#> 
#>   Level 7:   12 nodes to be scored   (96 eliminated genes)
#> 
#>   Level 6:   16 nodes to be scored   (173 eliminated genes)
#> 
#>   Level 5:   15 nodes to be scored   (239 eliminated genes)
#> 
#>   Level 4:   14 nodes to be scored   (425 eliminated genes)
#> 
#>   Level 3:   16 nodes to be scored   (459 eliminated genes)
#> 
#>   Level 2:   2 nodes to be scored    (496 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (525 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.A.HS.Sensitive.CC_weight01_10_all  --- no of nodes:  35 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: A.HS.Sensitive
#> Ontology: CC        GO.ID                         Term Annotated Significant Expected topGO
#> 1  GO:0005743 mitochondrial inner membrane         5           3     0.80 0.031
#> 2  GO:1902494            catalytic complex        42           8     6.71 0.032
#> 3  GO:0005739                mitochondrion        51          14     8.14 0.120
#> 4  GO:0005768                     endosome         8           3     1.28 0.121
#> 5  GO:0034399            nuclear periphery         8           3     1.28 0.121
#> 6  GO:0031966       mitochondrial membrane         6           4     0.96 0.155
#> 7  GO:0019866     organelle inner membrane         6           4     0.96 0.155
#> 8  GO:0030139            endocytic vesicle        17           3     2.71 0.161
#> 9  GO:0005929                       cilium         5           2     0.80 0.182
#> 10 GO:0020011                   apicoplast        42           9     6.71 0.211
#>    go.category interest.category algorithm statistic
#> 1           CC    A.HS.Sensitive  weight01    fisher
#> 2           CC    A.HS.Sensitive  weight01    fisher
#> 3           CC    A.HS.Sensitive  weight01    fisher
#> 4           CC    A.HS.Sensitive  weight01    fisher
#> 5           CC    A.HS.Sensitive  weight01    fisher
#> 6           CC    A.HS.Sensitive  weight01    fisher
#> 7           CC    A.HS.Sensitive  weight01    fisher
#> 8           CC    A.HS.Sensitive  weight01    fisher
#> 9           CC    A.HS.Sensitive  weight01    fisher
#> 10          CC    A.HS.Sensitive  weight01    fisher
#> 
#> ==============================================================================
#> 
#> interest-category 1 of 4
#> ontology 3 of 3
#> 
#> Building most specific GOs .....
#>  ( 283 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 554 GO terms and 715 relations. )
#> 
#> Annotating nodes ...............
#>  ( 334 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 52 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 9:   1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 8:   1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 7:   4 nodes to be scored    (11 eliminated genes)
#> 
#>   Level 6:   7 nodes to be scored    (12 eliminated genes)
#> 
#>   Level 5:   11 nodes to be scored   (35 eliminated genes)
#> 
#>   Level 4:   9 nodes to be scored    (57 eliminated genes)
#> 
#>   Level 3:   13 nodes to be scored   (77 eliminated genes)
#> 
#>   Level 2:   5 nodes to be scored    (169 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (299 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.and.CG.phenotype.MF_weight01_10_all  --- no of nodes:  21 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: HS.and.CG.phenotype
#> Ontology: MF        GO.ID                                        Term Annotated Significant
#> 1  GO:0003677                                 DNA binding        26           6
#> 2  GO:0008094       ATP-dependent activity, acting on DNA         8           3
#> 3  GO:0003690                 double-stranded DNA binding         9           3
#> 4  GO:1990837 sequence-specific double-stranded DNA bi...         7           2
#> 5  GO:0016773 phosphotransferase activity, alcohol gro...        25           2
#> 6  GO:0003688              DNA replication origin binding         5           1
#> 7  GO:0003697                 single-stranded DNA binding         5           1
#> 8  GO:0016301                             kinase activity        26           2
#> 9  GO:0003684                         damaged DNA binding         6           1
#> 10 GO:0030674      protein-macromolecule adaptor activity         6           1
#>    Expected  topGO go.category   interest.category algorithm statistic
#> 1      0.93 0.0012          MF HS.and.CG.phenotype  weight01    fisher
#> 2      0.29 0.0018          MF HS.and.CG.phenotype  weight01    fisher
#> 3      0.32 0.0603          MF HS.and.CG.phenotype  weight01    fisher
#> 4      0.25 0.0658          MF HS.and.CG.phenotype  weight01    fisher
#> 5      0.90 0.1643          MF HS.and.CG.phenotype  weight01    fisher
#> 6      0.18 0.1681          MF HS.and.CG.phenotype  weight01    fisher
#> 7      0.18 0.1681          MF HS.and.CG.phenotype  weight01    fisher
#> 8      0.93 0.1940          MF HS.and.CG.phenotype  weight01    fisher
#> 9      0.22 0.1985          MF HS.and.CG.phenotype  weight01    fisher
#> 10     0.22 0.1985          MF HS.and.CG.phenotype  weight01    fisher
#> 
#> ==============================================================================
#> 
#> interest-category 2 of 4
#> ontology 1 of 3
#> 
#> Building most specific GOs .....
#>  ( 358 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 1125 GO terms and 2077 relations. )
#> 
#> Annotating nodes ...............
#>  ( 354 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 84 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 10:  3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 8:   6 nodes to be scored    (24 eliminated genes)
#> 
#>   Level 7:   11 nodes to be scored   (53 eliminated genes)
#> 
#>   Level 6:   16 nodes to be scored   (81 eliminated genes)
#> 
#>   Level 5:   17 nodes to be scored   (132 eliminated genes)
#> 
#>   Level 4:   13 nodes to be scored   (208 eliminated genes)
#> 
#>   Level 3:   9 nodes to be scored    (257 eliminated genes)
#> 
#>   Level 2:   5 nodes to be scored    (287 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (295 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.and.CG.phenotype.BP_weight01_10_all  --- no of nodes:  48 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: HS.and.CG.phenotype
#> Ontology: BP        GO.ID                                    Term Annotated Significant
#> 1  GO:0065004            protein-DNA complex assembly         6           2
#> 2  GO:0006457                         protein folding        11           2
#> 3  GO:0043933 protein-containing complex organization        30           4
#> 4  GO:0006281                              DNA repair        18           2
#> 5  GO:0006270              DNA replication initiation         5           1
#> 6  GO:0006366      transcription by RNA polymerase II         5           1
#> 7  GO:0006891  intra-Golgi vesicle-mediated transport         5           1
#> 8  GO:0065007                   biological regulation        74           2
#> 9  GO:0006338                    chromatin remodeling         6           1
#> 10 GO:0019725                    cellular homeostasis         6           1
#>    Expected topGO go.category   interest.category algorithm statistic
#> 1      0.20 0.015          BP HS.and.CG.phenotype  weight01    fisher
#> 2      0.37 0.049          BP HS.and.CG.phenotype  weight01    fisher
#> 3      1.02 0.055          BP HS.and.CG.phenotype  weight01    fisher
#> 4      0.61 0.119          BP HS.and.CG.phenotype  weight01    fisher
#> 5      0.17 0.159          BP HS.and.CG.phenotype  weight01    fisher
#> 6      0.17 0.159          BP HS.and.CG.phenotype  weight01    fisher
#> 7      0.17 0.159          BP HS.and.CG.phenotype  weight01    fisher
#> 8      2.51 0.180          BP HS.and.CG.phenotype  weight01    fisher
#> 9      0.20 0.188          BP HS.and.CG.phenotype  weight01    fisher
#> 10     0.20 0.188          BP HS.and.CG.phenotype  weight01    fisher
#> 
#> ==============================================================================
#> 
#> interest-category 2 of 4
#> ontology 2 of 3
#> 
#> Building most specific GOs .....
#>  ( 195 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 367 GO terms and 616 relations. )
#> 
#> Annotating nodes ...............
#>  ( 526 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 80 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   5 nodes to be scored    (16 eliminated genes)
#> 
#>   Level 8:   8 nodes to be scored    (64 eliminated genes)
#> 
#>   Level 7:   8 nodes to be scored    (91 eliminated genes)
#> 
#>   Level 6:   14 nodes to be scored   (130 eliminated genes)
#> 
#>   Level 5:   14 nodes to be scored   (215 eliminated genes)
#> 
#>   Level 4:   11 nodes to be scored   (413 eliminated genes)
#> 
#>   Level 3:   13 nodes to be scored   (458 eliminated genes)
#> 
#>   Level 2:   2 nodes to be scored    (481 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (514 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.and.CG.phenotype.CC_weight01_10_all  --- no of nodes:  46 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: HS.and.CG.phenotype
#> Ontology: CC        GO.ID                                     Term Annotated Significant
#> 1  GO:0032993                      protein-DNA complex         5           3
#> 2  GO:0140513       nuclear protein-containing complex        41           4
#> 3  GO:0043231 intracellular membrane-bounded organelle       360          18
#> 4  GO:0005730                                nucleolus        16           2
#> 5  GO:0031982                                  vesicle        54           5
#> 6  GO:1903561                    extracellular vesicle        17           2
#> 7  GO:0005634                                  nucleus       238          13
#> 8  GO:0005743             mitochondrial inner membrane         5           1
#> 9  GO:0030120                             vesicle coat         5           1
#> 10 GO:0032040                 small-subunit processome         5           1
#>    Expected   topGO go.category   interest.category algorithm statistic
#> 1      0.20 0.00052          CC HS.and.CG.phenotype  weight01    fisher
#> 2      1.64 0.07194          CC HS.and.CG.phenotype  weight01    fisher
#> 3     14.37 0.11288          CC HS.and.CG.phenotype  weight01    fisher
#> 4      0.64 0.13010          CC HS.and.CG.phenotype  weight01    fisher
#> 5      2.16 0.13579          CC HS.and.CG.phenotype  weight01    fisher
#> 6      0.68 0.14398          CC HS.and.CG.phenotype  weight01    fisher
#> 7      9.50 0.15948          CC HS.and.CG.phenotype  weight01    fisher
#> 8      0.20 0.18495          CC HS.and.CG.phenotype  weight01    fisher
#> 9      0.20 0.18495          CC HS.and.CG.phenotype  weight01    fisher
#> 10     0.20 0.18495          CC HS.and.CG.phenotype  weight01    fisher
#> 
#> ==============================================================================
#> 
#> interest-category 2 of 4
#> ontology 3 of 3
#> 
#> Building most specific GOs .....
#>  ( 283 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 554 GO terms and 715 relations. )
#> 
#> Annotating nodes ...............
#>  ( 334 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 85 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 8:   4 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 7:   7 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 6:   11 nodes to be scored   (24 eliminated genes)
#> 
#>   Level 5:   16 nodes to be scored   (61 eliminated genes)
#> 
#>   Level 4:   20 nodes to be scored   (91 eliminated genes)
#> 
#>   Level 3:   20 nodes to be scored   (140 eliminated genes)
#> 
#>   Level 2:   6 nodes to be scored    (240 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (318 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.Neutral.MF_weight01_10_all  --- no of nodes:  35 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: HS.Neutral
#> Ontology: MF        GO.ID                                        Term Annotated Significant
#> 1  GO:0050839              cell adhesion molecule binding        17           5
#> 2  GO:0120545 nucleic acid conformation isomerase acti...        11           2
#> 3  GO:0004722 protein serine/threonine phosphatase act...         5           2
#> 4  GO:0003729                                mRNA binding        16           4
#> 5  GO:0140993                  histone modifying activity         6           2
#> 6  GO:0008270                            zinc ion binding         7           2
#> 7  GO:0008757 S-adenosylmethionine-dependent methyltra...         7           2
#> 8  GO:0097367             carbohydrate derivative binding        18           2
#> 9  GO:0016772 transferase activity, transferring phosp...        33           3
#> 10 GO:0008094       ATP-dependent activity, acting on DNA         8           2
#>    Expected topGO go.category interest.category algorithm statistic
#> 1      2.39 0.074          MF        HS.Neutral  weight01    fisher
#> 2      1.55 0.142          MF        HS.Neutral  weight01    fisher
#> 3      0.70 0.147          MF        HS.Neutral  weight01    fisher
#> 4      2.25 0.174          MF        HS.Neutral  weight01    fisher
#> 5      0.84 0.201          MF        HS.Neutral  weight01    fisher
#> 6      0.99 0.257          MF        HS.Neutral  weight01    fisher
#> 7      0.99 0.257          MF        HS.Neutral  weight01    fisher
#> 8      2.53 0.269          MF        HS.Neutral  weight01    fisher
#> 9      4.64 0.281          MF        HS.Neutral  weight01    fisher
#> 10     1.13 0.313          MF        HS.Neutral  weight01    fisher
#> 
#> ==============================================================================
#> 
#> interest-category 3 of 4
#> ontology 1 of 3
#> 
#> Building most specific GOs .....
#>  ( 358 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 1125 GO terms and 2077 relations. )
#> 
#> Annotating nodes ...............
#>  ( 354 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 198 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 13:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 12:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 11:  2 nodes to be scored    (11 eliminated genes)
#> 
#>   Level 10:  9 nodes to be scored    (11 eliminated genes)
#> 
#>   Level 9:   10 nodes to be scored   (17 eliminated genes)
#> 
#>   Level 8:   16 nodes to be scored   (77 eliminated genes)
#> 
#>   Level 7:   32 nodes to be scored   (95 eliminated genes)
#> 
#>   Level 6:   38 nodes to be scored   (140 eliminated genes)
#> 
#>   Level 5:   37 nodes to be scored   (235 eliminated genes)
#> 
#>   Level 4:   26 nodes to be scored   (277 eliminated genes)
#> 
#>   Level 3:   18 nodes to be scored   (318 eliminated genes)
#> 
#>   Level 2:   7 nodes to be scored    (343 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (348 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.Neutral.BP_weight01_10_all  --- no of nodes:  71 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: HS.Neutral
#> Ontology: BP        GO.ID                                        Term Annotated Significant
#> 1  GO:0072594 establishment of protein localization to...         6           4
#> 2  GO:0020033                         antigenic variation        27           8
#> 3  GO:0006355 regulation of DNA-templated transcriptio...        12           5
#> 4  GO:0006605                           protein targeting         8           4
#> 5  GO:0015931    nucleobase-containing compound transport         5           3
#> 6  GO:0033365           protein localization to organelle         5           3
#> 7  GO:0010558 negative regulation of macromolecule bio...         9           3
#> 8  GO:0098609                          cell-cell adhesion        19           5
#> 9  GO:0020013 symbiont-mediated perturbation of host e...        20           5
#> 10 GO:0006886             intracellular protein transport         9           3
#>    Expected  topGO go.category interest.category algorithm statistic
#> 1      0.66 0.0016          BP        HS.Neutral  weight01    fisher
#> 2      2.97 0.0049          BP        HS.Neutral  weight01    fisher
#> 3      1.32 0.0056          BP        HS.Neutral  weight01    fisher
#> 4      0.88 0.0064          BP        HS.Neutral  weight01    fisher
#> 5      0.55 0.0106          BP        HS.Neutral  weight01    fisher
#> 6      0.55 0.0106          BP        HS.Neutral  weight01    fisher
#> 7      0.99 0.0325          BP        HS.Neutral  weight01    fisher
#> 8      2.09 0.0456          BP        HS.Neutral  weight01    fisher
#> 9      2.20 0.0559          BP        HS.Neutral  weight01    fisher
#> 10     0.99 0.0651          BP        HS.Neutral  weight01    fisher
#> 
#> ==============================================================================
#> 
#> interest-category 3 of 4
#> ontology 2 of 3
#> 
#> Building most specific GOs .....
#>  ( 195 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 367 GO terms and 616 relations. )
#> 
#> Annotating nodes ...............
#>  ( 526 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 98 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   4 nodes to be scored    (16 eliminated genes)
#> 
#>   Level 8:   11 nodes to be scored   (75 eliminated genes)
#> 
#>   Level 7:   10 nodes to be scored   (87 eliminated genes)
#> 
#>   Level 6:   17 nodes to be scored   (146 eliminated genes)
#> 
#>   Level 5:   16 nodes to be scored   (233 eliminated genes)
#> 
#>   Level 4:   16 nodes to be scored   (425 eliminated genes)
#> 
#>   Level 3:   17 nodes to be scored   (465 eliminated genes)
#> 
#>   Level 2:   2 nodes to be scored    (499 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (525 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.Neutral.CC_weight01_10_all  --- no of nodes:  35 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: HS.Neutral
#> Ontology: CC        GO.ID                            Term Annotated Significant Expected
#> 1  GO:0020030 infected host cell surface knob        18           5     2.36
#> 2  GO:0005768                        endosome         8           3     1.05
#> 3  GO:0000785                       chromatin         9           3     1.18
#> 4  GO:0043657                       host cell        26           6     3.41
#> 5  GO:0035770       ribonucleoprotein granule         5           2     0.66
#> 6  GO:0020002       host cell plasma membrane        28           6     3.67
#> 7  GO:0020036                  Maurer's cleft        44           8     5.77
#> 8  GO:0016020                        membrane       102          14    13.38
#> 9  GO:0030659    cytoplasmic vesicle membrane         8           2     1.05
#> 10 GO:0030430             host cell cytoplasm        74          14     9.71
#>    topGO go.category interest.category algorithm statistic
#> 1  0.073          CC        HS.Neutral  weight01    fisher
#> 2  0.074          CC        HS.Neutral  weight01    fisher
#> 3  0.085          CC        HS.Neutral  weight01    fisher
#> 4  0.111          CC        HS.Neutral  weight01    fisher
#> 5  0.130          CC        HS.Neutral  weight01    fisher
#> 6  0.147          CC        HS.Neutral  weight01    fisher
#> 7  0.205          CC        HS.Neutral  weight01    fisher
#> 8  0.244          CC        HS.Neutral  weight01    fisher
#> 9  0.245          CC        HS.Neutral  weight01    fisher
#> 10 0.260          CC        HS.Neutral  weight01    fisher
#> 
#> ==============================================================================
#> 
#> interest-category 3 of 4
#> ontology 3 of 3
#> 
#> Building most specific GOs .....
#>  ( 283 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 554 GO terms and 715 relations. )
#> 
#> Annotating nodes ...............
#>  ( 334 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 112 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 9:   1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 8:   5 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 7:   11 nodes to be scored   (11 eliminated genes)
#> 
#>   Level 6:   14 nodes to be scored   (33 eliminated genes)
#> 
#>   Level 5:   24 nodes to be scored   (73 eliminated genes)
#> 
#>   Level 4:   24 nodes to be scored   (104 eliminated genes)
#> 
#>   Level 3:   24 nodes to be scored   (158 eliminated genes)
#> 
#>   Level 2:   8 nodes to be scored    (258 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (328 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.not.classified.MF_weight01_10_all  --- no of nodes:  44 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: not.classified
#> Ontology: MF        GO.ID                                        Term Annotated Significant
#> 1  GO:0004197        cysteine-type endopeptidase activity         6           6
#> 2  GO:0004843       cysteine-type deubiquitinase activity         5           5
#> 3  GO:0022890 inorganic cation transmembrane transport...         5           5
#> 4  GO:0008324 monoatomic cation transmembrane transpor...         5           5
#> 5  GO:0044877          protein-containing complex binding        13          11
#> 6  GO:0008047                   enzyme activator activity         8           7
#> 7  GO:0140097           catalytic activity, acting on DNA        12           7
#> 8  GO:0005524                                 ATP binding        11           9
#> 9  GO:0003735          structural constituent of ribosome        14          11
#> 10 GO:0016747 acyltransferase activity, transferring g...         7           6
#>    Expected topGO go.category interest.category algorithm statistic
#> 1      4.02 0.089          MF    not.classified  weight01    fisher
#> 2      3.35 0.134          MF    not.classified  weight01    fisher
#> 3      3.35 0.134          MF    not.classified  weight01    fisher
#> 4      3.35 0.134          MF    not.classified  weight01    fisher
#> 5      8.72 0.140          MF    not.classified  weight01    fisher
#> 6      5.37 0.198          MF    not.classified  weight01    fisher
#> 7      8.05 0.209          MF    not.classified  weight01    fisher
#> 8      7.38 0.238          MF    not.classified  weight01    fisher
#> 9      9.39 0.266          MF    not.classified  weight01    fisher
#> 10     4.69 0.268          MF    not.classified  weight01    fisher
#> 
#> ==============================================================================
#> 
#> interest-category 4 of 4
#> ontology 1 of 3
#> 
#> Building most specific GOs .....
#>  ( 358 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 1125 GO terms and 2077 relations. )
#> 
#> Annotating nodes ...............
#>  ( 354 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 251 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 13:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 12:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 11:  3 nodes to be scored    (11 eliminated genes)
#> 
#>   Level 10:  12 nodes to be scored   (11 eliminated genes)
#> 
#>   Level 9:   17 nodes to be scored   (24 eliminated genes)
#> 
#>   Level 8:   25 nodes to be scored   (98 eliminated genes)
#> 
#>   Level 7:   39 nodes to be scored   (125 eliminated genes)
#> 
#>   Level 6:   47 nodes to be scored   (174 eliminated genes)
#> 
#>   Level 5:   45 nodes to be scored   (261 eliminated genes)
#> 
#>   Level 4:   31 nodes to be scored   (293 eliminated genes)
#> 
#>   Level 3:   22 nodes to be scored   (329 eliminated genes)
#> 
#>   Level 2:   7 nodes to be scored    (345 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (353 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.not.classified.BP_weight01_10_all  --- no of nodes:  46 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: not.classified
#> Ontology: BP        GO.ID                                        Term Annotated Significant
#> 1  GO:0051604                          protein maturation        22          16
#> 2  GO:0006511 ubiquitin-dependent protein catabolic pr...        10           9
#> 3  GO:0006412                                 translation        26          21
#> 4  GO:0016579                    protein deubiquitination         5           5
#> 5  GO:0051640                      organelle localization         5           5
#> 6  GO:0048583          regulation of response to stimulus         5           5
#> 7  GO:0007165                         signal transduction        27          17
#> 8  GO:0051603 proteolysis involved in protein cataboli...        17          15
#> 9  GO:0070925                          organelle assembly        12          10
#> 10 GO:0044409                    symbiont entry into host        16          13
#>    Expected topGO go.category interest.category algorithm statistic
#> 1     15.35  0.12          BP    not.classified  weight01    fisher
#> 2      6.98  0.14          BP    not.classified  weight01    fisher
#> 3     18.14  0.15          BP    not.classified  weight01    fisher
#> 4      3.49  0.16          BP    not.classified  weight01    fisher
#> 5      3.49  0.16          BP    not.classified  weight01    fisher
#> 6      3.49  0.16          BP    not.classified  weight01    fisher
#> 7     18.84  0.18          BP    not.classified  weight01    fisher
#> 8     11.86  0.23          BP    not.classified  weight01    fisher
#> 9      8.37  0.23          BP    not.classified  weight01    fisher
#> 10    11.16  0.23          BP    not.classified  weight01    fisher
#> 
#> ==============================================================================
#> 
#> interest-category 4 of 4
#> ontology 2 of 3
#> 
#> Building most specific GOs .....
#>  ( 195 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 367 GO terms and 616 relations. )
#> 
#> Annotating nodes ...............
#>  ( 526 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 120 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  5 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   7 nodes to be scored    (16 eliminated genes)
#> 
#>   Level 8:   16 nodes to be scored   (86 eliminated genes)
#> 
#>   Level 7:   14 nodes to be scored   (102 eliminated genes)
#> 
#>   Level 6:   20 nodes to be scored   (183 eliminated genes)
#> 
#>   Level 5:   19 nodes to be scored   (251 eliminated genes)
#> 
#>   Level 4:   17 nodes to be scored   (428 eliminated genes)
#> 
#>   Level 3:   18 nodes to be scored   (465 eliminated genes)
#> 
#>   Level 2:   2 nodes to be scored    (500 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (525 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.not.classified.CC_weight01_10_all  --- no of nodes:  32 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: not.classified
#> Ontology: CC        GO.ID                                 Term Annotated Significant
#> 1  GO:0032991           protein-containing complex       125          88
#> 2  GO:0020003          symbiont-containing vacuole        18          14
#> 3  GO:0009986                         cell surface        24          20
#> 4  GO:0020009                            microneme         6           6
#> 5  GO:1990234                  transferase complex        22          18
#> 6  GO:0020002            host cell plasma membrane        28          22
#> 7  GO:0031410                  cytoplasmic vesicle        37          24
#> 8  GO:0071944                       cell periphery        26          19
#> 9  GO:0044164                    host cell cytosol         7           6
#> 10 GO:0043232 intracellular membraneless organelle        65          46
#>    Expected topGO go.category interest.category algorithm statistic
#> 1     83.65 0.016          CC    not.classified  weight01    fisher
#> 2     12.05 0.027          CC    not.classified  weight01    fisher
#> 3     16.06 0.058          CC    not.classified  weight01    fisher
#> 4      4.02 0.089          CC    not.classified  weight01    fisher
#> 5     14.72 0.103          CC    not.classified  weight01    fisher
#> 6     18.74 0.125          CC    not.classified  weight01    fisher
#> 7     24.76 0.205          CC    not.classified  weight01    fisher
#> 8     17.40 0.265          CC    not.classified  weight01    fisher
#> 9      4.68 0.266          CC    not.classified  weight01    fisher
#> 10    43.50 0.294          CC    not.classified  weight01    fisher
#> 
#> ==============================================================================
#> 
#> interest-category 4 of 4
#> ontology 3 of 3
#> Joining with `by = join_by(GO.ID, go.category, interest.category, algorithm,
#> statistic)`
#> 
#> ==============================================================================
#> Significant terms from final GO enrichment results table (saved to 'Routput/GO/all.combined.GO.results.tsv':
#>         GO.ID                                        Term Annotated Significant
#> 1  GO:0010468               regulation of gene expression        27           6
#> 2  GO:0044283         small molecule biosynthetic process         5           3
#> 3  GO:0061024                       membrane organization         9           4
#> 4  GO:0035556           intracellular signal transduction        22           7
#> 5  GO:0009410             response to xenobiotic stimulus        28           8
#> 6  GO:0008610                  lipid biosynthetic process        12           3
#> 7  GO:0006364                             rRNA processing        14           6
#> 8  GO:0006401                       RNA catabolic process         7           3
#> 9  GO:0000460                     maturation of 5.8S rRNA         7           3
#> 10 GO:0005743                mitochondrial inner membrane         5           3
#> 11 GO:1902494                           catalytic complex        42           8
#> 12 GO:0120543 macromolecular conformation isomerase ac...        14           3
#> 13 GO:0016740                        transferase activity        71          12
#> 14 GO:0072594 establishment of protein localization to...         6           4
#> 15 GO:0020033                         antigenic variation        27           8
#> 16 GO:0006355 regulation of DNA-templated transcriptio...        12           5
#> 17 GO:0006605                           protein targeting         8           4
#> 18 GO:0015931    nucleobase-containing compound transport         5           3
#> 19 GO:0033365           protein localization to organelle         5           3
#> 20 GO:0010558 negative regulation of macromolecule bio...         9           3
#> 21 GO:0098609                          cell-cell adhesion        19           5
#> 22 GO:0020013 symbiont-mediated perturbation of host e...        20           5
#> 23 GO:0006886             intracellular protein transport         9           3
#> 24 GO:0020030             infected host cell surface knob        18           5
#> 25 GO:0005768                                    endosome         8           3
#> 26 GO:0000785                                   chromatin         9           3
#> 27 GO:0050839              cell adhesion molecule binding        17           5
#> 28 GO:0065004                protein-DNA complex assembly         6           2
#> 29 GO:0006457                             protein folding        11           2
#> 30 GO:0043933     protein-containing complex organization        30           4
#> 31 GO:0032993                         protein-DNA complex         5           3
#> 32 GO:0140513          nuclear protein-containing complex        41           4
#> 33 GO:0003677                                 DNA binding        26           6
#> 34 GO:0008094       ATP-dependent activity, acting on DNA         8           3
#> 35 GO:0003690                 double-stranded DNA binding         9           3
#> 36 GO:1990837 sequence-specific double-stranded DNA bi...         7           2
#> 37 GO:0032991                  protein-containing complex       125          88
#> 38 GO:0020003                 symbiont-containing vacuole        18          14
#> 39 GO:0009986                                cell surface        24          20
#> 40 GO:0020009                                   microneme         6           6
#> 41 GO:0004197        cysteine-type endopeptidase activity         6           6
#>    Expected   topGO go.category   interest.category algorithm statistic
#> 1      4.27 0.00730          BP      A.HS.Sensitive  weight01    fisher
#> 2      0.79 0.02980          BP      A.HS.Sensitive  weight01    fisher
#> 3      1.42 0.03840          BP      A.HS.Sensitive  weight01    fisher
#> 4      3.48 0.04210          BP      A.HS.Sensitive  weight01    fisher
#> 5      4.43 0.05560          BP      A.HS.Sensitive  weight01    fisher
#> 6      1.90 0.06730          BP      A.HS.Sensitive  weight01    fisher
#> 7      2.21 0.07520          BP      A.HS.Sensitive  weight01    fisher
#> 8      1.11 0.08220          BP      A.HS.Sensitive  weight01    fisher
#> 9      1.11 0.08220          BP      A.HS.Sensitive  weight01    fisher
#> 10     0.80 0.03100          CC      A.HS.Sensitive  weight01    fisher
#> 11     6.71 0.03200          CC      A.HS.Sensitive  weight01    fisher
#> 12     2.14 0.06400          MF      A.HS.Sensitive  weight01    fisher
#> 13    10.84 0.07600          MF      A.HS.Sensitive  weight01    fisher
#> 14     0.66 0.00160          BP          HS.Neutral  weight01    fisher
#> 15     2.97 0.00490          BP          HS.Neutral  weight01    fisher
#> 16     1.32 0.00560          BP          HS.Neutral  weight01    fisher
#> 17     0.88 0.00640          BP          HS.Neutral  weight01    fisher
#> 18     0.55 0.01060          BP          HS.Neutral  weight01    fisher
#> 19     0.55 0.01060          BP          HS.Neutral  weight01    fisher
#> 20     0.99 0.03250          BP          HS.Neutral  weight01    fisher
#> 21     2.09 0.04560          BP          HS.Neutral  weight01    fisher
#> 22     2.20 0.05590          BP          HS.Neutral  weight01    fisher
#> 23     0.99 0.06510          BP          HS.Neutral  weight01    fisher
#> 24     2.36 0.07300          CC          HS.Neutral  weight01    fisher
#> 25     1.05 0.07400          CC          HS.Neutral  weight01    fisher
#> 26     1.18 0.08500          CC          HS.Neutral  weight01    fisher
#> 27     2.39 0.07400          MF          HS.Neutral  weight01    fisher
#> 28     0.20 0.01500          BP HS.and.CG.phenotype  weight01    fisher
#> 29     0.37 0.04900          BP HS.and.CG.phenotype  weight01    fisher
#> 30     1.02 0.05500          BP HS.and.CG.phenotype  weight01    fisher
#> 31     0.20 0.00052          CC HS.and.CG.phenotype  weight01    fisher
#> 32     1.64 0.07194          CC HS.and.CG.phenotype  weight01    fisher
#> 33     0.93 0.00120          MF HS.and.CG.phenotype  weight01    fisher
#> 34     0.29 0.00180          MF HS.and.CG.phenotype  weight01    fisher
#> 35     0.32 0.06030          MF HS.and.CG.phenotype  weight01    fisher
#> 36     0.25 0.06580          MF HS.and.CG.phenotype  weight01    fisher
#> 37    83.65 0.01600          CC      not.classified  weight01    fisher
#> 38    12.05 0.02700          CC      not.classified  weight01    fisher
#> 39    16.06 0.05800          CC      not.classified  weight01    fisher
#> 40     4.02 0.08900          CC      not.classified  weight01    fisher
#> 41     4.02 0.08900          MF      not.classified  weight01    fisher
#> 
#> ==============================================================================
#> 
#> 
#> 
#> All interesting-gene categories have been tested for GO-term enrichment.
#> 
#> See ALL TOP 30 enriched terms by interesting-gene category in 'Routput/GO/all.combined.GO.results.tsv'.
#> 
#> See log files for topGO-analyses by each interesting-gene category, including all genes in the analysis by GO term in 'Routput/GO/topGO.log.*.txt'.
#> 
#> See all significant genes mapped to all significant GO terms in 'Routput/GO/all.combined.sig.genes.per.sig.terms.tsv'.
#> 
#> ==============================================================================
```

``` r

# if you've run the pipeline, your significant genes in significant terms per category of interest can be loaded:
sig.genes <- read.delim("Routput/GO/all.combined.sig.genes.per.sig.terms.tsv")
```
