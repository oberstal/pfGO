
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pfGO

<!-- badges: start -->
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
#>  ( 335 GO terms and 438 relations. )
#> 
#> Annotating nodes ...............
#>  ( 193 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 58 nontrivial nodes
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
#>   Level 4:   10 nodes to be scored   (42 eliminated genes)
#> 
#>   Level 3:   15 nodes to be scored   (76 eliminated genes)
#> 
#>   Level 2:   5 nodes to be scored    (93 eliminated genes)
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
#> 9  GO:0016740                    transferase activity        26           5
#> 10 GO:0016301                         kinase activity         8           2
#>    Expected topGO go.category interest.category
#> 1      0.48 0.068          MF    A.HS.Sensitive
#> 2      1.12 0.084          MF    A.HS.Sensitive
#> 3      1.12 0.084          MF    A.HS.Sensitive
#> 4      2.73 0.244          MF    A.HS.Sensitive
#> 5      0.96 0.247          MF    A.HS.Sensitive
#> 6      4.66 0.248          MF    A.HS.Sensitive
#> 7      2.09 0.348          MF    A.HS.Sensitive
#> 8      1.28 0.377          MF    A.HS.Sensitive
#> 9      4.18 0.395          MF    A.HS.Sensitive
#> 10     1.28 0.408          MF    A.HS.Sensitive
#> 
#> ==============================================================================
#> 
#> interest-category 1 of 4
#> ontology 1 of 3
#> 
#> Building most specific GOs .....
#>  ( 115 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 551 GO terms and 1041 relations. )
#> 
#> Annotating nodes ...............
#>  ( 177 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 107 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 9:   1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 8:   3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 7:   6 nodes to be scored    (4 eliminated genes)
#> 
#>   Level 6:   16 nodes to be scored   (19 eliminated genes)
#> 
#>   Level 5:   29 nodes to be scored   (47 eliminated genes)
#> 
#>   Level 4:   25 nodes to be scored   (89 eliminated genes)
#> 
#>   Level 3:   19 nodes to be scored   (127 eliminated genes)
#> 
#>   Level 2:   7 nodes to be scored    (159 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (166 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.A.HS.Sensitive.BP_weight01_5_all  --- no of nodes:  23 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: A.HS.Sensitive
#> Ontology: BP        GO.ID                                        Term Annotated Significant
#> 1  GO:0016311                           dephosphorylation         3           2
#> 2  GO:0009410             response to xenobiotic stimulus        22           6
#> 3  GO:1901135 carbohydrate derivative metabolic proces...         5           2
#> 4  GO:0015031                           protein transport         5           2
#> 5  GO:0044237                  cellular metabolic process        63          10
#> 6  GO:0006996                      organelle organization         8           2
#> 7  GO:0007049                                  cell cycle         3           1
#> 8  GO:0048519 negative regulation of biological proces...         3           1
#> 9  GO:0030522    intracellular receptor signaling pathway         3           1
#> 10 GO:0044265    cellular macromolecule catabolic process         3           1
#>    Expected topGO go.category interest.category
#> 1      0.42 0.053          BP    A.HS.Sensitive
#> 2      3.11 0.066          BP    A.HS.Sensitive
#> 3      0.71 0.139          BP    A.HS.Sensitive
#> 4      0.71 0.139          BP    A.HS.Sensitive
#> 5      8.90 0.258          BP    A.HS.Sensitive
#> 6      1.13 0.315          BP    A.HS.Sensitive
#> 7      0.42 0.368          BP    A.HS.Sensitive
#> 8      0.42 0.368          BP    A.HS.Sensitive
#> 9      0.42 0.368          BP    A.HS.Sensitive
#> 10     0.42 0.368          BP    A.HS.Sensitive
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
#>  ( 246 GO terms and 418 relations. )
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
#>  ( 335 GO terms and 438 relations. )
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
#>  ( 115 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 551 GO terms and 1041 relations. )
#> 
#> Annotating nodes ...............
#>  ( 177 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 55 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 10:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   2 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 8:   3 nodes to be scored    (6 eliminated genes)
#> 
#>   Level 7:   3 nodes to be scored    (9 eliminated genes)
#> 
#>   Level 6:   5 nodes to be scored    (14 eliminated genes)
#> 
#>   Level 5:   13 nodes to be scored   (29 eliminated genes)
#> 
#>   Level 4:   12 nodes to be scored   (37 eliminated genes)
#> 
#>   Level 3:   10 nodes to be scored   (70 eliminated genes)
#> 
#>   Level 2:   5 nodes to be scored    (101 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (114 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.and.CG.phenotype.BP_weight01_5_all  --- no of nodes:  38 
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
#> 9  GO:0065007                       biological regulation        22           1
#> 10 GO:0016070                       RNA metabolic process        17           2
#>    Expected topGO go.category   interest.category
#> 1      0.10 0.099          BP HS.and.CG.phenotype
#> 2      0.10 0.099          BP HS.and.CG.phenotype
#> 3      0.14 0.130          BP HS.and.CG.phenotype
#> 4      0.14 0.130          BP HS.and.CG.phenotype
#> 5      0.17 0.160          BP HS.and.CG.phenotype
#> 6      0.20 0.189          BP HS.and.CG.phenotype
#> 7      0.20 0.189          BP HS.and.CG.phenotype
#> 8      0.27 0.245          BP HS.and.CG.phenotype
#> 9      0.75 0.555          BP HS.and.CG.phenotype
#> 10     0.58 1.000          BP HS.and.CG.phenotype
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
#>  ( 246 GO terms and 418 relations. )
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
#> 1  GO:0032993                  protein-DNA complex         3           2
#> 2  GO:0005730                            nucleolus         5           2
#> 3  GO:0030684                          preribosome         3           1
#> 4  GO:0030660    Golgi-associated vesicle membrane         3           1
#> 5  GO:1903561                extracellular vesicle        17           2
#> 6  GO:0000785                            chromatin         4           1
#> 7  GO:0030120                         vesicle coat         4           1
#> 8  GO:0005794                      Golgi apparatus         6           1
#> 9  GO:0005634                              nucleus       176           9
#> 10 GO:0020005 symbiont-containing vacuole membrane         9           1
#>    Expected  topGO go.category   interest.category
#> 1      0.11 0.0037          CC HS.and.CG.phenotype
#> 2      0.18 0.0120          CC HS.and.CG.phenotype
#> 3      0.11 0.1068          CC HS.and.CG.phenotype
#> 4      0.11 0.1068          CC HS.and.CG.phenotype
#> 5      0.63 0.1255          CC HS.and.CG.phenotype
#> 6      0.15 0.1400          CC HS.and.CG.phenotype
#> 7      0.15 0.1400          CC HS.and.CG.phenotype
#> 8      0.22 0.2029          CC HS.and.CG.phenotype
#> 9      6.49 0.2353          CC HS.and.CG.phenotype
#> 10     0.33 0.2891          CC HS.and.CG.phenotype
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
#>  ( 335 GO terms and 438 relations. )
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
#>   Level 2:   4 nodes to be scored    (118 eliminated genes)
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
#>  ( 115 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 551 GO terms and 1041 relations. )
#> 
#> Annotating nodes ...............
#>  ( 177 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 149 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  4 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   7 nodes to be scored    (5 eliminated genes)
#> 
#>   Level 8:   10 nodes to be scored   (34 eliminated genes)
#> 
#>   Level 7:   16 nodes to be scored   (40 eliminated genes)
#> 
#>   Level 6:   25 nodes to be scored   (53 eliminated genes)
#> 
#>   Level 5:   33 nodes to be scored   (82 eliminated genes)
#> 
#>   Level 4:   24 nodes to be scored   (107 eliminated genes)
#> 
#>   Level 3:   21 nodes to be scored   (130 eliminated genes)
#> 
#>   Level 2:   7 nodes to be scored    (159 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (167 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.Neutral.BP_weight01_5_all  --- no of nodes:  79 
#> 
#> ==============================================================================
#> GO enrichment results for interest category: HS.Neutral
#> Ontology: BP        GO.ID                                        Term Annotated Significant
#> 1  GO:0006355 regulation of DNA-templated transcriptio...         5           4
#> 2  GO:0020033                         antigenic variation        27           8
#> 3  GO:0006913                 nucleocytoplasmic transport         3           2
#> 4  GO:0051276                     chromosome organization         3           2
#> 5  GO:0020013 modulation by symbiont of host erythrocy...        20           5
#> 6  GO:0020035    adhesion of symbiont to microvasculature        22           5
#> 7  GO:0051701 biological process involved in interacti...        53          10
#> 8  GO:0071705                 nitrogen compound transport         7           2
#> 9  GO:0006259                       DNA metabolic process         8           2
#> 10 GO:0006508                                 proteolysis         8           2
#>    Expected topGO go.category interest.category
#> 1      0.65 0.001          BP        HS.Neutral
#> 2      3.51 0.011          BP        HS.Neutral
#> 3      0.39 0.045          BP        HS.Neutral
#> 4      0.39 0.045          BP        HS.Neutral
#> 5      2.60 0.095          BP        HS.Neutral
#> 6      2.86 0.134          BP        HS.Neutral
#> 7      6.89 0.211          BP        HS.Neutral
#> 8      0.91 0.240          BP        HS.Neutral
#> 9      1.04 0.241          BP        HS.Neutral
#> 10     1.04 0.278          BP        HS.Neutral
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
#>  ( 246 GO terms and 418 relations. )
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
#>   Level 3:   15 nodes to be scored   (364 eliminated genes)
#> 
#>   Level 2:   2 nodes to be scored    (413 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (432 eliminated genes)
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
#>  ( 335 GO terms and 438 relations. )
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
#>   Level 2:   7 nodes to be scored    (135 eliminated genes)
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
#>  ( 115 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 551 GO terms and 1041 relations. )
#> 
#> Annotating nodes ...............
#>  ( 177 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 203 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  5 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   11 nodes to be scored   (5 eliminated genes)
#> 
#>   Level 8:   13 nodes to be scored   (37 eliminated genes)
#> 
#>   Level 7:   21 nodes to be scored   (54 eliminated genes)
#> 
#>   Level 6:   36 nodes to be scored   (68 eliminated genes)
#> 
#>   Level 5:   43 nodes to be scored   (101 eliminated genes)
#> 
#>   Level 4:   33 nodes to be scored   (130 eliminated genes)
#> 
#>   Level 3:   29 nodes to be scored   (140 eliminated genes)
#> 
#>   Level 2:   10 nodes to be scored   (166 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (177 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.not.classified.BP_weight01_5_all  --- no of nodes:  36 
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
#> 8  GO:0030003 intracellular monoatomic cation homeosta...         3           3
#> 9  GO:0050776               regulation of immune response         3           3
#> 10 GO:0051301                               cell division         3           3
#>    Expected topGO go.category interest.category
#> 1      3.47  0.16          BP    not.classified
#> 2      2.78  0.23          BP    not.classified
#> 3     15.29  0.28          BP    not.classified
#> 4      7.64  0.29          BP    not.classified
#> 5      4.86  0.31          BP    not.classified
#> 6      2.08  0.33          BP    not.classified
#> 7      2.08  0.33          BP    not.classified
#> 8      2.08  0.33          BP    not.classified
#> 9      2.08  0.33          BP    not.classified
#> 10     2.08  0.33          BP    not.classified
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
#>  ( 246 GO terms and 418 relations. )
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
#> 6  GO:0032991                  protein-containing complex        54          35
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
#> 6     35.59 0.128          CC    not.classified
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
#> Final GO enrichment results table (saved to 'Routput/GO/all.combined.GO.results.tsv':
#>          GO.ID                                        Term Annotated
#> 1   GO:0016311                           dephosphorylation         3
#> 2   GO:0009410             response to xenobiotic stimulus        22
#> 3   GO:1901135 carbohydrate derivative metabolic proces...         5
#> 4   GO:0015031                           protein transport         5
#> 5   GO:0044237                  cellular metabolic process        63
#> 6   GO:0006996                      organelle organization         8
#> 7   GO:0007049                                  cell cycle         3
#> 8   GO:0048519 negative regulation of biological proces...         3
#> 9   GO:0030522    intracellular receptor signaling pathway         3
#> 10  GO:0044265    cellular macromolecule catabolic process         3
#> 11  GO:0005739                               mitochondrion        46
#> 12  GO:0031966                      mitochondrial membrane         3
#> 13  GO:0140513          nuclear protein-containing complex        12
#> 14  GO:0034399                           nuclear periphery         8
#> 15  GO:0005737                                   cytoplasm       176
#> 16  GO:0005840                                    ribosome         8
#> 17  GO:0005783                       endoplasmic reticulum        19
#> 18  GO:0020011                                  apicoplast        42
#> 19  GO:0005794                             Golgi apparatus         6
#> 20  GO:0031981                               nuclear lumen        15
#> 21  GO:0003924                             GTPase activity         3
#> 22  GO:0016791                        phosphatase activity         7
#> 23  GO:0003735          structural constituent of ribosome         7
#> 24  GO:0140096     catalytic activity, acting on a protein        17
#> 25  GO:0016887                     ATP hydrolysis activity         6
#> 26  GO:0016787                          hydrolase activity        29
#> 27  GO:0003729                                mRNA binding        13
#> 28  GO:0022857          transmembrane transporter activity         8
#> 29  GO:0016740                        transferase activity        26
#> 30  GO:0016301                             kinase activity         8
#> 31  GO:0006355 regulation of DNA-templated transcriptio...         5
#> 32  GO:0020033                         antigenic variation        27
#> 33  GO:0006913                 nucleocytoplasmic transport         3
#> 34  GO:0051276                     chromosome organization         3
#> 35  GO:0020013 modulation by symbiont of host erythrocy...        20
#> 36  GO:0020035    adhesion of symbiont to microvasculature        22
#> 37  GO:0051701 biological process involved in interacti...        53
#> 38  GO:0071705                 nitrogen compound transport         7
#> 39  GO:0006259                       DNA metabolic process         8
#> 40  GO:0006508                                 proteolysis         8
#> 41  GO:0020030             infected host cell surface knob        18
#> 42  GO:0020002                   host cell plasma membrane        28
#> 43  GO:0030430                         host cell cytoplasm        80
#> 44  GO:0020036                              Maurer's cleft        43
#> 45  GO:0005634                                     nucleus       176
#> 46  GO:0020008                                     rhoptry         8
#> 47  GO:0034399                           nuclear periphery         8
#> 48  GO:0140535    intracellular protein-containing complex         8
#> 49  GO:0043657                                   host cell       105
#> 50  GO:0020005        symbiont-containing vacuole membrane         9
#> 51  GO:0043565               sequence-specific DNA binding         3
#> 52  GO:0008757 S-adenosylmethionine-dependent methyltra...         3
#> 53  GO:0016853                          isomerase activity         3
#> 54  GO:0008094       ATP-dependent activity, acting on DNA         3
#> 55  GO:0016779             nucleotidyltransferase activity         4
#> 56  GO:0050839              cell adhesion molecule binding        17
#> 57  GO:0003677                                 DNA binding         8
#> 58  GO:0046872                           metal ion binding         5
#> 59  GO:0003729                                mRNA binding        13
#> 60  GO:0016409               palmitoyltransferase activity         3
#> 61  GO:0045454                      cell redox homeostasis         3
#> 62  GO:0034470                            ncRNA processing         3
#> 63  GO:0006457                             protein folding         4
#> 64  GO:0006888 endoplasmic reticulum to Golgi vesicle-m...         4
#> 65  GO:0065003         protein-containing complex assembly         5
#> 66  GO:0042254                         ribosome biogenesis         6
#> 67  GO:0006351                 DNA-templated transcription         6
#> 68  GO:0006259                       DNA metabolic process         8
#> 69  GO:0065007                       biological regulation        22
#> 70  GO:0016070                       RNA metabolic process        17
#> 71  GO:0032993                         protein-DNA complex         3
#> 72  GO:0005730                                   nucleolus         5
#> 73  GO:0030684                                 preribosome         3
#> 74  GO:0030660           Golgi-associated vesicle membrane         3
#> 75  GO:1903561                       extracellular vesicle        17
#> 76  GO:0000785                                   chromatin         4
#> 77  GO:0030120                                vesicle coat         4
#> 78  GO:0005794                             Golgi apparatus         6
#> 79  GO:0005634                                     nucleus       176
#> 80  GO:0020005        symbiont-containing vacuole membrane         9
#> 81  GO:0019899                              enzyme binding         3
#> 82  GO:0022804 active transmembrane transporter activit...         3
#> 83  GO:0016791                        phosphatase activity         7
#> 84  GO:0003723                                 RNA binding        35
#> 85  GO:0005515                             protein binding        79
#> 86  GO:0003674                          molecular_function       193
#> 87  GO:0003676                        nucleic acid binding        47
#> 88  GO:0005488                                     binding       128
#> 89  GO:0005215                        transporter activity        11
#> 90  GO:0042578         phosphoric ester hydrolase activity         7
#> 91  GO:0048870                               cell motility         5
#> 92  GO:0006413                    translational initiation         4
#> 93  GO:0020035    adhesion of symbiont to microvasculature        22
#> 94  GO:0044409                             entry into host        11
#> 95  GO:0035891                         exit from host cell         7
#> 96  GO:0032502                       developmental process         3
#> 97  GO:0002377                   immunoglobulin production         3
#> 98  GO:0030003 intracellular monoatomic cation homeosta...         3
#> 99  GO:0050776               regulation of immune response         3
#> 100 GO:0051301                               cell division         3
#> 101 GO:0020003                 symbiont-containing vacuole        18
#> 102 GO:0009986                                cell surface        24
#> 103 GO:0020009                                   microneme         6
#> 104 GO:0071944                              cell periphery        19
#> 105 GO:0020002                   host cell plasma membrane        28
#> 106 GO:0032991                  protein-containing complex        54
#> 107 GO:0044228                           host cell surface         4
#> 108 GO:0031982                                     vesicle        41
#> 109 GO:0020026                     merozoite dense granule         3
#> 110 GO:0005852 eukaryotic translation initiation factor...         3
#> 111 GO:0042393                             histone binding         6
#> 112 GO:0016874                             ligase activity         4
#> 113 GO:0008289                               lipid binding         4
#> 114 GO:0140101        catalytic activity, acting on a tRNA         4
#> 115 GO:0140030      modification-dependent protein binding         4
#> 116 GO:0003743      translation initiation factor activity         4
#> 117 GO:0005515                             protein binding        79
#> 118 GO:0008233                          peptidase activity         6
#> 119 GO:0016765 transferase activity, transferring alkyl...         3
#> 120 GO:0140097           catalytic activity, acting on DNA         6
#>     Significant Expected  topGO go.category   interest.category
#> 1             2     0.42 0.0530          BP      A.HS.Sensitive
#> 2             6     3.11 0.0660          BP      A.HS.Sensitive
#> 3             2     0.71 0.1390          BP      A.HS.Sensitive
#> 4             2     0.71 0.1390          BP      A.HS.Sensitive
#> 5            10     8.90 0.2580          BP      A.HS.Sensitive
#> 6             2     1.13 0.3150          BP      A.HS.Sensitive
#> 7             1     0.42 0.3680          BP      A.HS.Sensitive
#> 8             1     0.42 0.3680          BP      A.HS.Sensitive
#> 9             1     0.42 0.3680          BP      A.HS.Sensitive
#> 10            1     0.42 0.3680          BP      A.HS.Sensitive
#> 11           13     7.63 0.0680          CC      A.HS.Sensitive
#> 12            2     0.50 0.0730          CC      A.HS.Sensitive
#> 13            4     1.99 0.1210          CC      A.HS.Sensitive
#> 14            3     1.33 0.1320          CC      A.HS.Sensitive
#> 15           39    29.20 0.1630          CC      A.HS.Sensitive
#> 16            3     1.33 0.1640          CC      A.HS.Sensitive
#> 17            5     3.15 0.1920          CC      A.HS.Sensitive
#> 18            9     6.97 0.2450          CC      A.HS.Sensitive
#> 19            2     1.00 0.2610          CC      A.HS.Sensitive
#> 20            5     2.49 0.2970          CC      A.HS.Sensitive
#> 21            2     0.48 0.0680          MF      A.HS.Sensitive
#> 22            3     1.12 0.0840          MF      A.HS.Sensitive
#> 23            3     1.12 0.0840          MF      A.HS.Sensitive
#> 24            4     2.73 0.2440          MF      A.HS.Sensitive
#> 25            2     0.96 0.2470          MF      A.HS.Sensitive
#> 26           10     4.66 0.2480          MF      A.HS.Sensitive
#> 27            3     2.09 0.3480          MF      A.HS.Sensitive
#> 28            2     1.28 0.3770          MF      A.HS.Sensitive
#> 29            5     4.18 0.3950          MF      A.HS.Sensitive
#> 30            2     1.28 0.4080          MF      A.HS.Sensitive
#> 31            4     0.65 0.0010          BP          HS.Neutral
#> 32            8     3.51 0.0110          BP          HS.Neutral
#> 33            2     0.39 0.0450          BP          HS.Neutral
#> 34            2     0.39 0.0450          BP          HS.Neutral
#> 35            5     2.60 0.0950          BP          HS.Neutral
#> 36            5     2.86 0.1340          BP          HS.Neutral
#> 37           10     6.89 0.2110          BP          HS.Neutral
#> 38            2     0.91 0.2400          BP          HS.Neutral
#> 39            2     1.04 0.2410          BP          HS.Neutral
#> 40            2     1.04 0.2780          BP          HS.Neutral
#> 41            5     2.49 0.0870          CC          HS.Neutral
#> 42            6     3.87 0.1750          CC          HS.Neutral
#> 43           15    11.06 0.2280          CC          HS.Neutral
#> 44            8     5.94 0.2280          CC          HS.Neutral
#> 45           27    24.33 0.2580          CC          HS.Neutral
#> 46            2     1.11 0.3050          CC          HS.Neutral
#> 47            2     1.11 0.3050          CC          HS.Neutral
#> 48            2     1.11 0.3050          CC          HS.Neutral
#> 49           19    14.52 0.3340          CC          HS.Neutral
#> 50            2     1.24 0.3600          CC          HS.Neutral
#> 51            2     0.51 0.0760          MF          HS.Neutral
#> 52            2     0.51 0.0760          MF          HS.Neutral
#> 53            2     0.51 0.0760          MF          HS.Neutral
#> 54            2     0.51 0.0760          MF          HS.Neutral
#> 55            2     0.68 0.1360          MF          HS.Neutral
#> 56            5     2.91 0.1420          MF          HS.Neutral
#> 57            4     1.37 0.1880          MF          HS.Neutral
#> 58            2     0.85 0.2030          MF          HS.Neutral
#> 59            3     2.22 0.3890          MF          HS.Neutral
#> 60            1     0.51 0.4320          MF          HS.Neutral
#> 61            1     0.10 0.0990          BP HS.and.CG.phenotype
#> 62            1     0.10 0.0990          BP HS.and.CG.phenotype
#> 63            1     0.14 0.1300          BP HS.and.CG.phenotype
#> 64            1     0.14 0.1300          BP HS.and.CG.phenotype
#> 65            1     0.17 0.1600          BP HS.and.CG.phenotype
#> 66            1     0.20 0.1890          BP HS.and.CG.phenotype
#> 67            1     0.20 0.1890          BP HS.and.CG.phenotype
#> 68            1     0.27 0.2450          BP HS.and.CG.phenotype
#> 69            1     0.75 0.5550          BP HS.and.CG.phenotype
#> 70            2     0.58 1.0000          BP HS.and.CG.phenotype
#> 71            2     0.11 0.0037          CC HS.and.CG.phenotype
#> 72            2     0.18 0.0120          CC HS.and.CG.phenotype
#> 73            1     0.11 0.1068          CC HS.and.CG.phenotype
#> 74            1     0.11 0.1068          CC HS.and.CG.phenotype
#> 75            2     0.63 0.1255          CC HS.and.CG.phenotype
#> 76            1     0.15 0.1400          CC HS.and.CG.phenotype
#> 77            1     0.15 0.1400          CC HS.and.CG.phenotype
#> 78            1     0.22 0.2029          CC HS.and.CG.phenotype
#> 79            9     6.49 0.2353          CC HS.and.CG.phenotype
#> 80            1     0.33 0.2891          CC HS.and.CG.phenotype
#> 81            1     0.06 0.0610          MF HS.and.CG.phenotype
#> 82            1     0.06 0.0610          MF HS.and.CG.phenotype
#> 83            1     0.15 0.1380          MF HS.and.CG.phenotype
#> 84            1     0.73 0.5540          MF HS.and.CG.phenotype
#> 85            2     1.64 0.7860          MF HS.and.CG.phenotype
#> 86            4     4.00 1.0000          MF HS.and.CG.phenotype
#> 87            1     0.97 1.0000          MF HS.and.CG.phenotype
#> 88            3     2.65 1.0000          MF HS.and.CG.phenotype
#> 89            1     0.23 1.0000          MF HS.and.CG.phenotype
#> 90            1     0.15 1.0000          MF HS.and.CG.phenotype
#> 91            5     3.47 0.1600          BP      not.classified
#> 92            4     2.78 0.2300          BP      not.classified
#> 93           17    15.29 0.2800          BP      not.classified
#> 94            9     7.64 0.2900          BP      not.classified
#> 95            6     4.86 0.3100          BP      not.classified
#> 96            3     2.08 0.3300          BP      not.classified
#> 97            3     2.08 0.3300          BP      not.classified
#> 98            3     2.08 0.3300          BP      not.classified
#> 99            3     2.08 0.3300          BP      not.classified
#> 100           3     2.08 0.3300          BP      not.classified
#> 101          14    11.86 0.0230          CC      not.classified
#> 102          20    15.82 0.0460          CC      not.classified
#> 103           6     3.95 0.0800          CC      not.classified
#> 104          14    12.52 0.0810          CC      not.classified
#> 105          22    18.45 0.1020          CC      not.classified
#> 106          35    35.59 0.1280          CC      not.classified
#> 107           4     2.64 0.1870          CC      not.classified
#> 108          29    27.02 0.2820          CC      not.classified
#> 109           3     1.98 0.2850          CC      not.classified
#> 110           3     1.98 0.2850          CC      not.classified
#> 111           6     3.89 0.0710          MF      not.classified
#> 112           4     2.59 0.1730          MF      not.classified
#> 113           4     2.59 0.1730          MF      not.classified
#> 114           4     2.59 0.1730          MF      not.classified
#> 115           4     2.59 0.1730          MF      not.classified
#> 116           4     2.59 0.1730          MF      not.classified
#> 117          55    51.17 0.2620          MF      not.classified
#> 118           5     3.89 0.2690          MF      not.classified
#> 119           3     1.94 0.2690          MF      not.classified
#> 120           4     3.89 0.2760          MF      not.classified
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
