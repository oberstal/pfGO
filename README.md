
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
#> No documentation for 'run.topGO.meta' in specified packages and libraries:
#> you could try '??run.topGO.meta'
?Pfal_geneID2GO
#> No documentation for 'Pfal_geneID2GO' in specified packages and libraries:
#> you could try '??Pfal_geneID2GO'

# for all available functions/data:
?pfGO
#> No documentation for 'pfGO' in specified packages and libraries:
#> you could try '??pfGO'
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

#### parameters:

- *mydf*: data frame with geneIDs in column 1, and interest-category
  classifications in column 2.
- *geneID2GO*: a list of named vectors of GO IDs–one vector of GO-terms
  for each geneID.
- *pval*: p-value threshold for significance. Defaults to 0.05.

#### outputs:

run.topGO.meta creates several output-files, including:

- enrichment results

- significant genes per significant term

  - available in “Routput/GO/all.combined.sig.genes.per.sig.terms.tsv”

- plots of the GO-term hierarchy relevant to the analysis

- thorough log-files for each gene-category of interest tested against
  the background of all other genes in the analysis

- system logs recording all package versions, etc.

Primary results from run.topGO.meta will be in
“Routput/GO/all.combined.GO.results.tsv”. Note that run.topGO.meta will
automatically create the Routput directory (and other required output
directories nested in ./Routput) in your working directory for you if it
does not exist.

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
biological process, cellular compartment; MF, BP, and CC, respectively)
sequentially on all groups of interest. Results are combined in the
final output-table (“Routput/GO/all.combined.GO.results.tsv”).

TopGO automatically accounts for genes that cannot be mapped to GO terms
(or are mapped to terms with \< 3 genes in the analysis) with “feasible
genes” indicated in the topGO.log files in the “Routput/GO” folder.

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
#>  ( 128 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 336 GO terms and 437 relations. )
#> 
#> Annotating nodes ...............
#>  ( 193 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 59 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 8:   3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 7:   3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 6:   10 nodes to be scored   (10 eliminated genes)
#> 
#>   Level 5:   11 nodes to be scored   (10 eliminated genes)
#> 
#>   Level 4:   11 nodes to be scored   (42 eliminated genes)
#> 
#>   Level 3:   15 nodes to be scored   (76 eliminated genes)
#> 
#>   Level 2:   5 nodes to be scored    (97 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (183 eliminated genes)
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
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.A.HS.Sensitive.MF_weight01_5_all  --- no of nodes:  16 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: A.HS.Sensitive
#> Ontology: MF        GO.ID                                    Term Annotated Significant
#> 1  GO:0003924                         GTPase activity         3           2
#> 2  GO:0016791                    phosphatase activity         7           3
#> 3  GO:0003735      structural constituent of ribosome         7           3
#> 4  GO:0140096 catalytic activity, acting on a protein        17           4
#> 5  GO:0016887                 ATP hydrolysis activity         6           2
#> 6  GO:0016787                      hydrolase activity        29          10
#> 7  GO:0003729                            mRNA binding        13           3
#> 8  GO:0022857      transmembrane transporter activity         8           2
#> 9  GO:0016301                         kinase activity         8           2
#> 10 GO:0004540                   RNA nuclease activity         3           1
#>    Expected topGO go.category interest.category
#> 1      0.48 0.068          MF    A.HS.Sensitive
#> 2      1.12 0.084          MF    A.HS.Sensitive
#> 3      1.12 0.084          MF    A.HS.Sensitive
#> 4      2.73 0.244          MF    A.HS.Sensitive
#> 5      0.96 0.247          MF    A.HS.Sensitive
#> 6      4.66 0.248          MF    A.HS.Sensitive
#> 7      2.09 0.348          MF    A.HS.Sensitive
#> 8      1.28 0.377          MF    A.HS.Sensitive
#> 9      1.28 0.408          MF    A.HS.Sensitive
#> 10     0.48 0.410          MF    A.HS.Sensitive
#> 
#> ==============================================================================
#> 
#> interest-category 1 of 4
#> ontology 1 of 3
#> 
#> Building most specific GOs .....
#>  ( 114 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 535 GO terms and 1001 relations. )
#> 
#> Annotating nodes ...............
#>  ( 176 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 105 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 9:   1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 8:   4 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 7:   7 nodes to be scored    (4 eliminated genes)
#> 
#>   Level 6:   15 nodes to be scored   (22 eliminated genes)
#> 
#>   Level 5:   27 nodes to be scored   (50 eliminated genes)
#> 
#>   Level 4:   26 nodes to be scored   (93 eliminated genes)
#> 
#>   Level 3:   18 nodes to be scored   (130 eliminated genes)
#> 
#>   Level 2:   6 nodes to be scored    (162 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (165 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.A.HS.Sensitive.BP_weight01_5_all  --- no of nodes:  23 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: A.HS.Sensitive
#> Ontology: BP        GO.ID                                        Term Annotated Significant
#> 1  GO:0016311                           dephosphorylation         3           2
#> 2  GO:0009410             response to xenobiotic stimulus        22           6
#> 3  GO:1901135 carbohydrate derivative metabolic proces...         5           2
#> 4  GO:0015031                           protein transport         5           2
#> 5  GO:0044237                  cellular metabolic process        66          10
#> 6  GO:0006996                      organelle organization         8           2
#> 7  GO:0006468                     protein phosphorylation         9           2
#> 8  GO:0007049                                  cell cycle         3           1
#> 9  GO:0048519 negative regulation of biological proces...         3           1
#> 10 GO:0030522    intracellular receptor signaling pathway         3           1
#>    Expected topGO go.category interest.category
#> 1      0.41 0.049          BP    A.HS.Sensitive
#> 2      3.00 0.056          BP    A.HS.Sensitive
#> 3      0.68 0.134          BP    A.HS.Sensitive
#> 4      0.68 0.134          BP    A.HS.Sensitive
#> 5      9.00 0.251          BP    A.HS.Sensitive
#> 6      1.09 0.299          BP    A.HS.Sensitive
#> 7      1.23 0.354          BP    A.HS.Sensitive
#> 8      0.41 0.358          BP    A.HS.Sensitive
#> 9      0.41 0.358          BP    A.HS.Sensitive
#> 10     0.41 0.358          BP    A.HS.Sensitive
#> 
#> ==============================================================================
#> 
#> interest-category 1 of 4
#> ontology 2 of 3
#> 
#> Building most specific GOs .....
#>  ( 100 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 245 GO terms and 416 relations. )
#> 
#> Annotating nodes ...............
#>  ( 434 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 80 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  2 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  4 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   5 nodes to be scored    (23 eliminated genes)
#> 
#>   Level 8:   8 nodes to be scored    (70 eliminated genes)
#> 
#>   Level 7:   9 nodes to be scored    (93 eliminated genes)
#> 
#>   Level 6:   13 nodes to be scored   (115 eliminated genes)
#> 
#>   Level 5:   10 nodes to be scored   (171 eliminated genes)
#> 
#>   Level 4:   11 nodes to be scored   (342 eliminated genes)
#> 
#>   Level 3:   15 nodes to be scored   (358 eliminated genes)
#> 
#>   Level 2:   2 nodes to be scored    (404 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (432 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.A.HS.Sensitive.CC_weight01_5_all  --- no of nodes:  23 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: A.HS.Sensitive
#> Ontology: CC        GO.ID                               Term Annotated Significant Expected
#> 1  GO:0005739                      mitochondrion        46          13     7.63
#> 2  GO:0031966             mitochondrial membrane         3           2     0.50
#> 3  GO:0140513 nuclear protein-containing complex        12           4     1.99
#> 4  GO:0034399                  nuclear periphery         8           3     1.33
#> 5  GO:0005737                          cytoplasm       176          39    29.20
#> 6  GO:0005840                           ribosome         8           3     1.33
#> 7  GO:0005783              endoplasmic reticulum        19           5     3.15
#> 8  GO:0020011                         apicoplast        42           9     6.97
#> 9  GO:0005794                    Golgi apparatus         6           2     1.00
#> 10 GO:0031981                      nuclear lumen        15           5     2.49
#>    topGO go.category interest.category
#> 1  0.068          CC    A.HS.Sensitive
#> 2  0.073          CC    A.HS.Sensitive
#> 3  0.121          CC    A.HS.Sensitive
#> 4  0.132          CC    A.HS.Sensitive
#> 5  0.163          CC    A.HS.Sensitive
#> 6  0.164          CC    A.HS.Sensitive
#> 7  0.192          CC    A.HS.Sensitive
#> 8  0.245          CC    A.HS.Sensitive
#> 9  0.261          CC    A.HS.Sensitive
#> 10 0.297          CC    A.HS.Sensitive
#> 
#> ==============================================================================
#> 
#> interest-category 1 of 4
#> ontology 3 of 3
#> 
#> Building most specific GOs .....
#>  ( 128 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 336 GO terms and 437 relations. )
#> 
#> Annotating nodes ...............
#>  ( 193 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 16 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 6:   1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 5:   2 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 4:   4 nodes to be scored    (7 eliminated genes)
#> 
#>   Level 3:   5 nodes to be scored    (42 eliminated genes)
#> 
#>   Level 2:   3 nodes to be scored    (64 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (143 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.and.CG.phenotype.MF_weight01_5_all  --- no of nodes:  16 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: HS.and.CG.phenotype
#> Ontology: MF        GO.ID                                        Term Annotated Significant
#> 1  GO:0019899                              enzyme binding         3           1
#> 2  GO:0022804 active transmembrane transporter activit...         3           1
#> 3  GO:0016791                        phosphatase activity         7           1
#> 4  GO:0003723                                 RNA binding        35           1
#> 5  GO:0005515                             protein binding        79           2
#> 6  GO:0003674                          molecular_function       193           4
#> 7  GO:0003676                        nucleic acid binding        47           1
#> 8  GO:0005488                                     binding       128           3
#> 9  GO:0005215                        transporter activity        11           1
#> 10 GO:0042578         phosphoric ester hydrolase activity         7           1
#>    Expected topGO go.category   interest.category
#> 1      0.06 0.061          MF HS.and.CG.phenotype
#> 2      0.06 0.061          MF HS.and.CG.phenotype
#> 3      0.15 0.138          MF HS.and.CG.phenotype
#> 4      0.73 0.554          MF HS.and.CG.phenotype
#> 5      1.64 0.786          MF HS.and.CG.phenotype
#> 6      4.00 1.000          MF HS.and.CG.phenotype
#> 7      0.97 1.000          MF HS.and.CG.phenotype
#> 8      2.65 1.000          MF HS.and.CG.phenotype
#> 9      0.23 1.000          MF HS.and.CG.phenotype
#> 10     0.15 1.000          MF HS.and.CG.phenotype
#> 
#> ==============================================================================
#> 
#> interest-category 2 of 4
#> ontology 1 of 3
#> 
#> Building most specific GOs .....
#>  ( 114 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 535 GO terms and 1001 relations. )
#> 
#> Annotating nodes ...............
#>  ( 176 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 56 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 9:   2 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 8:   4 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 7:   4 nodes to be scored    (9 eliminated genes)
#> 
#>   Level 6:   6 nodes to be scored    (18 eliminated genes)
#> 
#>   Level 5:   12 nodes to be scored   (33 eliminated genes)
#> 
#>   Level 4:   13 nodes to be scored   (50 eliminated genes)
#> 
#>   Level 3:   9 nodes to be scored    (88 eliminated genes)
#> 
#>   Level 2:   5 nodes to be scored    (105 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (113 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.and.CG.phenotype.BP_weight01_5_all  --- no of nodes:  45 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: HS.and.CG.phenotype
#> Ontology: BP        GO.ID                                        Term Annotated Significant
#> 1  GO:0045454                      cell redox homeostasis         3           1
#> 2  GO:0034470                            ncRNA processing         3           1
#> 3  GO:0006457                             protein folding         4           1
#> 4  GO:0006888 endoplasmic reticulum to Golgi vesicle-m...         4           1
#> 5  GO:0065003         protein-containing complex assembly         5           1
#> 6  GO:0042254                         ribosome biogenesis         6           1
#> 7  GO:0006351                 DNA-templated transcription         6           1
#> 8  GO:0006259                       DNA metabolic process         8           1
#> 9  GO:0065007                       biological regulation        21           1
#> 10 GO:1901564 organonitrogen compound metabolic proces...        44           1
#>    Expected topGO go.category   interest.category
#> 1      0.10 0.099          BP HS.and.CG.phenotype
#> 2      0.10 0.099          BP HS.and.CG.phenotype
#> 3      0.14 0.131          BP HS.and.CG.phenotype
#> 4      0.14 0.131          BP HS.and.CG.phenotype
#> 5      0.17 0.161          BP HS.and.CG.phenotype
#> 6      0.20 0.190          BP HS.and.CG.phenotype
#> 7      0.20 0.190          BP HS.and.CG.phenotype
#> 8      0.27 0.247          BP HS.and.CG.phenotype
#> 9      0.72 0.539          BP HS.and.CG.phenotype
#> 10     1.50 1.000          BP HS.and.CG.phenotype
#> 
#> ==============================================================================
#> 
#> interest-category 2 of 4
#> ontology 2 of 3
#> 
#> Building most specific GOs .....
#>  ( 100 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 245 GO terms and 416 relations. )
#> 
#> Annotating nodes ...............
#>  ( 434 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 68 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  2 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  4 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   5 nodes to be scored    (23 eliminated genes)
#> 
#>   Level 8:   7 nodes to be scored    (70 eliminated genes)
#> 
#>   Level 7:   7 nodes to be scored    (90 eliminated genes)
#> 
#>   Level 6:   12 nodes to be scored   (101 eliminated genes)
#> 
#>   Level 5:   8 nodes to be scored    (165 eliminated genes)
#> 
#>   Level 4:   10 nodes to be scored   (335 eliminated genes)
#> 
#>   Level 3:   10 nodes to be scored   (352 eliminated genes)
#> 
#>   Level 2:   2 nodes to be scored    (389 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (416 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.and.CG.phenotype.CC_weight01_5_all  --- no of nodes:  34 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: HS.and.CG.phenotype
#> Ontology: CC        GO.ID                                 Term Annotated Significant
#> 1  GO:0005730                            nucleolus         5           2
#> 2  GO:0032993                  protein-DNA complex         6           2
#> 3  GO:0030684                          preribosome         3           1
#> 4  GO:0030660    Golgi-associated vesicle membrane         3           1
#> 5  GO:1903561                extracellular vesicle        17           2
#> 6  GO:0000785                            chromatin         4           1
#> 7  GO:0030120                         vesicle coat         4           1
#> 8  GO:0005794                      Golgi apparatus         6           1
#> 9  GO:0005634                              nucleus       176           9
#> 10 GO:0020005 symbiont-containing vacuole membrane         9           1
#>    Expected topGO go.category   interest.category
#> 1      0.18 0.012          CC HS.and.CG.phenotype
#> 2      0.22 0.069          CC HS.and.CG.phenotype
#> 3      0.11 0.107          CC HS.and.CG.phenotype
#> 4      0.11 0.107          CC HS.and.CG.phenotype
#> 5      0.63 0.126          CC HS.and.CG.phenotype
#> 6      0.15 0.140          CC HS.and.CG.phenotype
#> 7      0.15 0.140          CC HS.and.CG.phenotype
#> 8      0.22 0.203          CC HS.and.CG.phenotype
#> 9      6.49 0.235          CC HS.and.CG.phenotype
#> 10     0.33 0.289          CC HS.and.CG.phenotype
#> 
#> ==============================================================================
#> 
#> interest-category 2 of 4
#> ontology 3 of 3
#> 
#> Building most specific GOs .....
#>  ( 128 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 336 GO terms and 437 relations. )
#> 
#> Annotating nodes ...............
#>  ( 193 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 49 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 8:   1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 7:   1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 6:   6 nodes to be scored    (6 eliminated genes)
#> 
#>   Level 5:   9 nodes to be scored    (8 eliminated genes)
#> 
#>   Level 4:   14 nodes to be scored   (36 eliminated genes)
#> 
#>   Level 3:   13 nodes to be scored   (71 eliminated genes)
#> 
#>   Level 2:   4 nodes to be scored    (119 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (169 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.Neutral.MF_weight01_5_all  --- no of nodes:  19 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: HS.Neutral
#> Ontology: MF        GO.ID                                        Term Annotated Significant
#> 1  GO:0043565               sequence-specific DNA binding         3           2
#> 2  GO:0008757 S-adenosylmethionine-dependent methyltra...         3           2
#> 3  GO:0016853                          isomerase activity         3           2
#> 4  GO:0008094       ATP-dependent activity, acting on DNA         3           2
#> 5  GO:0016779             nucleotidyltransferase activity         4           2
#> 6  GO:0050839              cell adhesion molecule binding        17           5
#> 7  GO:0003677                                 DNA binding         8           4
#> 8  GO:0046872                           metal ion binding         5           2
#> 9  GO:0003729                                mRNA binding        13           3
#> 10 GO:0016409               palmitoyltransferase activity         3           1
#>    Expected topGO go.category interest.category
#> 1      0.51 0.076          MF        HS.Neutral
#> 2      0.51 0.076          MF        HS.Neutral
#> 3      0.51 0.076          MF        HS.Neutral
#> 4      0.51 0.076          MF        HS.Neutral
#> 5      0.68 0.136          MF        HS.Neutral
#> 6      2.91 0.142          MF        HS.Neutral
#> 7      1.37 0.188          MF        HS.Neutral
#> 8      0.85 0.203          MF        HS.Neutral
#> 9      2.22 0.389          MF        HS.Neutral
#> 10     0.51 0.432          MF        HS.Neutral
#> 
#> ==============================================================================
#> 
#> interest-category 3 of 4
#> ontology 1 of 3
#> 
#> Building most specific GOs .....
#>  ( 114 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 535 GO terms and 1001 relations. )
#> 
#> Annotating nodes ...............
#>  ( 176 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 146 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 10:  3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   7 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 8:   12 nodes to be scored   (33 eliminated genes)
#> 
#>   Level 7:   17 nodes to be scored   (40 eliminated genes)
#> 
#>   Level 6:   22 nodes to be scored   (54 eliminated genes)
#> 
#>   Level 5:   32 nodes to be scored   (82 eliminated genes)
#> 
#>   Level 4:   24 nodes to be scored   (112 eliminated genes)
#> 
#>   Level 3:   21 nodes to be scored   (135 eliminated genes)
#> 
#>   Level 2:   7 nodes to be scored    (162 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (170 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.Neutral.BP_weight01_5_all  --- no of nodes:  77 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: HS.Neutral
#> Ontology: BP        GO.ID                                        Term Annotated Significant
#> 1  GO:0006355 regulation of DNA-templated transcriptio...         5           4
#> 2  GO:0020033                         antigenic variation        27           8
#> 3  GO:0006913                 nucleocytoplasmic transport         3           2
#> 4  GO:0051276                     chromosome organization         3           2
#> 5  GO:0020013 modulation by symbiont of host erythrocy...        20           5
#> 6  GO:0060255 regulation of macromolecule metabolic pr...         8           6
#> 7  GO:0020035    adhesion of symbiont to microvasculature        22           5
#> 8  GO:0010468               regulation of gene expression         7           5
#> 9  GO:0051171 regulation of nitrogen compound metaboli...         7           5
#> 10 GO:0051701 biological process involved in interacti...        53          10
#>    Expected topGO go.category interest.category
#> 1      0.65 0.001          BP        HS.Neutral
#> 2      3.53 0.011          BP        HS.Neutral
#> 3      0.39 0.045          BP        HS.Neutral
#> 4      0.39 0.045          BP        HS.Neutral
#> 5      2.61 0.097          BP        HS.Neutral
#> 6      1.05 0.107          BP        HS.Neutral
#> 7      2.88 0.137          BP        HS.Neutral
#> 8      0.91 0.210          BP        HS.Neutral
#> 9      0.91 0.210          BP        HS.Neutral
#> 10     6.93 0.212          BP        HS.Neutral
#> 
#> ==============================================================================
#> 
#> interest-category 3 of 4
#> ontology 2 of 3
#> 
#> Building most specific GOs .....
#>  ( 100 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 245 GO terms and 416 relations. )
#> 
#> Annotating nodes ...............
#>  ( 434 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 69 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  2 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  4 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   3 nodes to be scored    (23 eliminated genes)
#> 
#>   Level 8:   4 nodes to be scored    (82 eliminated genes)
#> 
#>   Level 7:   7 nodes to be scored    (86 eliminated genes)
#> 
#>   Level 6:   10 nodes to be scored   (100 eliminated genes)
#> 
#>   Level 5:   8 nodes to be scored    (164 eliminated genes)
#> 
#>   Level 4:   12 nodes to be scored   (341 eliminated genes)
#> 
#>   Level 3:   16 nodes to be scored   (364 eliminated genes)
#> 
#>   Level 2:   2 nodes to be scored    (413 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (433 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.Neutral.CC_weight01_5_all  --- no of nodes:  19 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: HS.Neutral
#> Ontology: CC        GO.ID                                     Term Annotated Significant
#> 1  GO:0020030          infected host cell surface knob        18           5
#> 2  GO:0020002                host cell plasma membrane        28           6
#> 3  GO:0030430                      host cell cytoplasm        80          15
#> 4  GO:0020036                           Maurer's cleft        43           8
#> 5  GO:0005634                                  nucleus       176          27
#> 6  GO:0020008                                  rhoptry         8           2
#> 7  GO:0034399                        nuclear periphery         8           2
#> 8  GO:0140535 intracellular protein-containing complex         8           2
#> 9  GO:0043657                                host cell       105          19
#> 10 GO:0020005     symbiont-containing vacuole membrane         9           2
#>    Expected topGO go.category interest.category
#> 1      2.49 0.087          CC        HS.Neutral
#> 2      3.87 0.175          CC        HS.Neutral
#> 3     11.06 0.228          CC        HS.Neutral
#> 4      5.94 0.228          CC        HS.Neutral
#> 5     24.33 0.258          CC        HS.Neutral
#> 6      1.11 0.305          CC        HS.Neutral
#> 7      1.11 0.305          CC        HS.Neutral
#> 8      1.11 0.305          CC        HS.Neutral
#> 9     14.52 0.334          CC        HS.Neutral
#> 10     1.24 0.360          CC        HS.Neutral
#> 
#> ==============================================================================
#> 
#> interest-category 3 of 4
#> ontology 3 of 3
#> 
#> Building most specific GOs .....
#>  ( 128 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 336 GO terms and 437 relations. )
#> 
#> Annotating nodes ...............
#>  ( 193 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 84 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 8:   3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 7:   4 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 6:   12 nodes to be scored   (10 eliminated genes)
#> 
#>   Level 5:   17 nodes to be scored   (14 eliminated genes)
#> 
#>   Level 4:   21 nodes to be scored   (50 eliminated genes)
#> 
#>   Level 3:   19 nodes to be scored   (93 eliminated genes)
#> 
#>   Level 2:   7 nodes to be scored    (136 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (187 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.not.classified.MF_weight01_5_all  --- no of nodes:  11 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: not.classified
#> Ontology: MF        GO.ID                                        Term Annotated Significant
#> 1  GO:0042393                             histone binding         6           6
#> 2  GO:0016874                             ligase activity         4           4
#> 3  GO:0008289                               lipid binding         4           4
#> 4  GO:0140101        catalytic activity, acting on a tRNA         4           4
#> 5  GO:0140030      modification-dependent protein binding         4           4
#> 6  GO:0003743      translation initiation factor activity         4           4
#> 7  GO:0005515                             protein binding        79          55
#> 8  GO:0008233                          peptidase activity         6           5
#> 9  GO:0016765 transferase activity, transferring alkyl...         3           3
#> 10 GO:0140097           catalytic activity, acting on DNA         6           4
#>    Expected topGO go.category interest.category
#> 1      3.89 0.071          MF    not.classified
#> 2      2.59 0.173          MF    not.classified
#> 3      2.59 0.173          MF    not.classified
#> 4      2.59 0.173          MF    not.classified
#> 5      2.59 0.173          MF    not.classified
#> 6      2.59 0.173          MF    not.classified
#> 7     51.17 0.262          MF    not.classified
#> 8      3.89 0.269          MF    not.classified
#> 9      1.94 0.269          MF    not.classified
#> 10     3.89 0.276          MF    not.classified
#> 
#> ==============================================================================
#> 
#> interest-category 4 of 4
#> ontology 1 of 3
#> 
#> Building most specific GOs .....
#>  ( 114 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 535 GO terms and 1001 relations. )
#> 
#> Annotating nodes ...............
#>  ( 176 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 200 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 10:  4 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   11 nodes to be scored   (0 eliminated genes)
#> 
#>   Level 8:   17 nodes to be scored   (36 eliminated genes)
#> 
#>   Level 7:   23 nodes to be scored   (54 eliminated genes)
#> 
#>   Level 6:   31 nodes to be scored   (74 eliminated genes)
#> 
#>   Level 5:   41 nodes to be scored   (105 eliminated genes)
#> 
#>   Level 4:   34 nodes to be scored   (133 eliminated genes)
#> 
#>   Level 3:   29 nodes to be scored   (143 eliminated genes)
#> 
#>   Level 2:   9 nodes to be scored    (169 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (176 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.not.classified.BP_weight01_5_all  --- no of nodes:  34 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: not.classified
#> Ontology: BP        GO.ID                                        Term Annotated Significant
#> 1  GO:0048870                               cell motility         5           5
#> 2  GO:0006413                    translational initiation         4           4
#> 3  GO:0020035    adhesion of symbiont to microvasculature        22          17
#> 4  GO:0044409                             entry into host        11           9
#> 5  GO:0035891                         exit from host cell         7           6
#> 6  GO:0032502                       developmental process         3           3
#> 7  GO:0002377                   immunoglobulin production         3           3
#> 8  GO:0065008            regulation of biological quality         3           3
#> 9  GO:0030003 intracellular monoatomic cation homeosta...         3           3
#> 10 GO:0050776               regulation of immune response         3           3
#>    Expected topGO go.category interest.category
#> 1      3.49  0.16          BP    not.classified
#> 2      2.80  0.24          BP    not.classified
#> 3     15.38  0.29          BP    not.classified
#> 4      7.69  0.30          BP    not.classified
#> 5      4.89  0.32          BP    not.classified
#> 6      2.10  0.34          BP    not.classified
#> 7      2.10  0.34          BP    not.classified
#> 8      2.10  0.34          BP    not.classified
#> 9      2.10  0.34          BP    not.classified
#> 10     2.10  0.34          BP    not.classified
#> 
#> ==============================================================================
#> 
#> interest-category 4 of 4
#> ontology 2 of 3
#> 
#> Building most specific GOs .....
#>  ( 100 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 245 GO terms and 416 relations. )
#> 
#> Annotating nodes ...............
#>  ( 434 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 105 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  2 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  6 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   7 nodes to be scored    (23 eliminated genes)
#> 
#>   Level 8:   13 nodes to be scored   (90 eliminated genes)
#> 
#>   Level 7:   12 nodes to be scored   (97 eliminated genes)
#> 
#>   Level 6:   16 nodes to be scored   (122 eliminated genes)
#> 
#>   Level 5:   13 nodes to be scored   (182 eliminated genes)
#> 
#>   Level 4:   15 nodes to be scored   (346 eliminated genes)
#> 
#>   Level 3:   18 nodes to be scored   (364 eliminated genes)
#> 
#>   Level 2:   2 nodes to be scored    (413 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (433 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.not.classified.CC_weight01_5_all  --- no of nodes:  25 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: not.classified
#> Ontology: CC        GO.ID                                        Term Annotated Significant
#> 1  GO:0020003                 symbiont-containing vacuole        18          14
#> 2  GO:0009986                                cell surface        24          20
#> 3  GO:0020009                                   microneme         6           6
#> 4  GO:0071944                              cell periphery        19          14
#> 5  GO:0020002                   host cell plasma membrane        28          22
#> 6  GO:0032991                  protein-containing complex        55          35
#> 7  GO:0044228                           host cell surface         4           4
#> 8  GO:0031982                                     vesicle        41          29
#> 9  GO:0020026                     merozoite dense granule         3           3
#> 10 GO:0005852 eukaryotic translation initiation factor...         3           3
#>    Expected topGO go.category interest.category
#> 1     11.86 0.023          CC    not.classified
#> 2     15.82 0.046          CC    not.classified
#> 3      3.95 0.080          CC    not.classified
#> 4     12.52 0.081          CC    not.classified
#> 5     18.45 0.102          CC    not.classified
#> 6     36.24 0.130          CC    not.classified
#> 7      2.64 0.187          CC    not.classified
#> 8     27.02 0.282          CC    not.classified
#> 9      1.98 0.285          CC    not.classified
#> 10     1.98 0.285          CC    not.classified
#> 
#> ==============================================================================
#> 
#> interest-category 4 of 4
#> ontology 3 of 3
#> Joining with `by = join_by(GO.ID, go.category, interest.category)`
#> 
#> ==============================================================================
#> Significant terms from final GO enrichment results table (saved to 'Routput/GO/all.combined.GO.results.tsv':
#>         GO.ID                                        Term Annotated Significant
#> 1  GO:0016311                           dephosphorylation         3           2
#> 2  GO:0009410             response to xenobiotic stimulus        22           6
#> 3  GO:0005739                               mitochondrion        46          13
#> 4  GO:0031966                      mitochondrial membrane         3           2
#> 5  GO:0003924                             GTPase activity         3           2
#> 6  GO:0016791                        phosphatase activity         7           3
#> 7  GO:0003735          structural constituent of ribosome         7           3
#> 8  GO:0006355 regulation of DNA-templated transcriptio...         5           4
#> 9  GO:0020033                         antigenic variation        27           8
#> 10 GO:0006913                 nucleocytoplasmic transport         3           2
#> 11 GO:0051276                     chromosome organization         3           2
#> 12 GO:0020013 modulation by symbiont of host erythrocy...        20           5
#> 13 GO:0020030             infected host cell surface knob        18           5
#> 14 GO:0043565               sequence-specific DNA binding         3           2
#> 15 GO:0008757 S-adenosylmethionine-dependent methyltra...         3           2
#> 16 GO:0016853                          isomerase activity         3           2
#> 17 GO:0008094       ATP-dependent activity, acting on DNA         3           2
#> 18 GO:0045454                      cell redox homeostasis         3           1
#> 19 GO:0034470                            ncRNA processing         3           1
#> 20 GO:0005730                                   nucleolus         5           2
#> 21 GO:0032993                         protein-DNA complex         6           2
#> 22 GO:0019899                              enzyme binding         3           1
#> 23 GO:0022804 active transmembrane transporter activit...         3           1
#> 24 GO:0020003                 symbiont-containing vacuole        18          14
#> 25 GO:0009986                                cell surface        24          20
#> 26 GO:0020009                                   microneme         6           6
#> 27 GO:0071944                              cell periphery        19          14
#> 28 GO:0042393                             histone binding         6           6
#>    Expected topGO go.category   interest.category
#> 1      0.41 0.049          BP      A.HS.Sensitive
#> 2      3.00 0.056          BP      A.HS.Sensitive
#> 3      7.63 0.068          CC      A.HS.Sensitive
#> 4      0.50 0.073          CC      A.HS.Sensitive
#> 5      0.48 0.068          MF      A.HS.Sensitive
#> 6      1.12 0.084          MF      A.HS.Sensitive
#> 7      1.12 0.084          MF      A.HS.Sensitive
#> 8      0.65 0.001          BP          HS.Neutral
#> 9      3.53 0.011          BP          HS.Neutral
#> 10     0.39 0.045          BP          HS.Neutral
#> 11     0.39 0.045          BP          HS.Neutral
#> 12     2.61 0.097          BP          HS.Neutral
#> 13     2.49 0.087          CC          HS.Neutral
#> 14     0.51 0.076          MF          HS.Neutral
#> 15     0.51 0.076          MF          HS.Neutral
#> 16     0.51 0.076          MF          HS.Neutral
#> 17     0.51 0.076          MF          HS.Neutral
#> 18     0.10 0.099          BP HS.and.CG.phenotype
#> 19     0.10 0.099          BP HS.and.CG.phenotype
#> 20     0.18 0.012          CC HS.and.CG.phenotype
#> 21     0.22 0.069          CC HS.and.CG.phenotype
#> 22     0.06 0.061          MF HS.and.CG.phenotype
#> 23     0.06 0.061          MF HS.and.CG.phenotype
#> 24    11.86 0.023          CC      not.classified
#> 25    15.82 0.046          CC      not.classified
#> 26     3.95 0.080          CC      not.classified
#> 27    12.52 0.081          CC      not.classified
#> 28     3.89 0.071          MF      not.classified
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
