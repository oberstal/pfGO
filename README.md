
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

You can install the development version of pfGO from
[GitHub](https://github.com/) with:

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

## Example

This is a basic example demonstrating how to run an enrichment analysis:

``` r
library(pfGO)
```

``` r
# load included pf GO database and example-data to be tested for functional enrichment
data("Pfal_geneID2GO")
data("exampleMydf")

# run the topGO pipeline on all experimental categories of interest from exampleMydf
run.topGO.meta(mydf = exampleMydf, geneID2GO = Pfal_geneID2GO, pval = 0.05)
#> 
#> Building most specific GOs .....
#>  ( 286 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 558 GO terms and 720 relations. )
#> 
#> Annotating nodes ...............
#>  ( 339 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 119 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  2 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   6 nodes to be scored    (3 eliminated genes)
#> 
#>   Level 8:   9 nodes to be scored    (8 eliminated genes)
#> 
#>   Level 7:   12 nodes to be scored   (32 eliminated genes)
#> 
#>   Level 6:   20 nodes to be scored   (56 eliminated genes)
#> 
#>   Level 5:   21 nodes to be scored   (86 eliminated genes)
#> 
#>   Level 4:   23 nodes to be scored   (130 eliminated genes)
#> 
#>   Level 3:   18 nodes to be scored   (195 eliminated genes)
#> 
#>   Level 2:   6 nodes to be scored    (239 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (324 eliminated genes)
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
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.A.HS.Sensitive.MF_weight01_5_all  --- no of nodes:  21 
#> 'data.frame':    0 obs. of  0 variables
#> NULL
#> data frame with 0 columns and 0 rows
#> 
#> Building most specific GOs .....
#>  ( 384 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 1342 GO terms and 2675 relations. )
#> 
#> Annotating nodes ...............
#>  ( 347 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 313 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 12:  2 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 11:  7 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  11 nodes to be scored   (14 eliminated genes)
#> 
#>   Level 9:   26 nodes to be scored   (48 eliminated genes)
#> 
#>   Level 8:   33 nodes to be scored   (81 eliminated genes)
#> 
#>   Level 7:   39 nodes to be scored   (129 eliminated genes)
#> 
#>   Level 6:   60 nodes to be scored   (199 eliminated genes)
#> 
#>   Level 5:   61 nodes to be scored   (234 eliminated genes)
#> 
#>   Level 4:   39 nodes to be scored   (278 eliminated genes)
#> 
#>   Level 3:   27 nodes to be scored   (314 eliminated genes)
#> 
#>   Level 2:   7 nodes to be scored    (336 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (344 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.A.HS.Sensitive.BP_weight01_5_all  --- no of nodes:  57 
#> 'data.frame':    6 obs. of  2 variables:
#>  $ GO.ID: Factor w/ 6 levels "PF3D7_0623100",..: 1 2 3 4 5 6
#>  $ 1    : chr  "PF3D7_0623100" "PF3D7_0804000" "PF3D7_1314900" "PF3D7_1316000" ...
#> NULL
#>           GO.ID             1
#> 1 PF3D7_0623100 PF3D7_0623100
#> 2 PF3D7_0804000 PF3D7_0804000
#> 3 PF3D7_1314900 PF3D7_1314900
#> 4 PF3D7_1316000 PF3D7_1316000
#> 5 PF3D7_1445500 PF3D7_1445500
#> 6 PF3D7_1468200 PF3D7_1468200
#> 
#> Building most specific GOs .....
#>  ( 200 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 380 GO terms and 646 relations. )
#> 
#> Annotating nodes ...............
#>  ( 518 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 134 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  6 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   11 nodes to be scored   (26 eliminated genes)
#> 
#>   Level 8:   17 nodes to be scored   (80 eliminated genes)
#> 
#>   Level 7:   18 nodes to be scored   (122 eliminated genes)
#> 
#>   Level 6:   21 nodes to be scored   (184 eliminated genes)
#> 
#>   Level 5:   21 nodes to be scored   (247 eliminated genes)
#> 
#>   Level 4:   16 nodes to be scored   (419 eliminated genes)
#> 
#>   Level 3:   18 nodes to be scored   (452 eliminated genes)
#> 
#>   Level 2:   2 nodes to be scored    (497 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (518 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.A.HS.Sensitive.CC_weight01_5_all  --- no of nodes:  36 
#> 'data.frame':    0 obs. of  0 variables
#> NULL
#> data frame with 0 columns and 0 rows
#> 
#> Building most specific GOs .....
#>  ( 286 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 558 GO terms and 720 relations. )
#> 
#> Annotating nodes ...............
#>  ( 339 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 59 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 9:   3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 8:   3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 7:   5 nodes to be scored    (22 eliminated genes)
#> 
#>   Level 6:   11 nodes to be scored   (38 eliminated genes)
#> 
#>   Level 5:   11 nodes to be scored   (62 eliminated genes)
#> 
#>   Level 4:   10 nodes to be scored   (87 eliminated genes)
#> 
#>   Level 3:   12 nodes to be scored   (157 eliminated genes)
#> 
#>   Level 2:   3 nodes to be scored    (181 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (301 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.and.CG.phenotype.MF_weight01_5_all  --- no of nodes:  19 
#> 'data.frame':    11 obs. of  4 variables:
#>  $ GO.ID     : Factor w/ 3 levels "GO:0008094","GO:0003677",..: 1 1 1 2 2 2 2 2 2 3 ...
#>  $ GO:0008094: chr  "PF3D7_0706700" "PF3D7_0818700" "PF3D7_1211300" NA ...
#>  $ GO:0003677: chr  NA NA NA "PF3D7_0110800" ...
#>  $ GO:1990837: chr  NA NA NA NA ...
#> NULL
#>         GO.ID    GO:0008094    GO:0003677    GO:1990837
#> 1  GO:0008094 PF3D7_0706700          <NA>          <NA>
#> 2  GO:0008094 PF3D7_0818700          <NA>          <NA>
#> 3  GO:0008094 PF3D7_1211300          <NA>          <NA>
#> 4  GO:0003677          <NA> PF3D7_0110800          <NA>
#> 5  GO:0003677          <NA> PF3D7_0706700          <NA>
#> 6  GO:0003677          <NA> PF3D7_0714000          <NA>
#> 7  GO:0003677          <NA> PF3D7_0818700          <NA>
#> 8  GO:0003677          <NA> PF3D7_1211300          <NA>
#> 9  GO:0003677          <NA> PF3D7_1363200          <NA>
#> 10 GO:1990837          <NA>          <NA> PF3D7_0110800
#> 11 GO:1990837          <NA>          <NA> PF3D7_1211300
#> 
#> Building most specific GOs .....
#>  ( 384 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 1342 GO terms and 2675 relations. )
#> 
#> Annotating nodes ...............
#>  ( 347 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 127 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   6 nodes to be scored    (5 eliminated genes)
#> 
#>   Level 8:   13 nodes to be scored   (32 eliminated genes)
#> 
#>   Level 7:   14 nodes to be scored   (65 eliminated genes)
#> 
#>   Level 6:   19 nodes to be scored   (121 eliminated genes)
#> 
#>   Level 5:   30 nodes to be scored   (185 eliminated genes)
#> 
#>   Level 4:   20 nodes to be scored   (229 eliminated genes)
#> 
#>   Level 3:   15 nodes to be scored   (266 eliminated genes)
#> 
#>   Level 2:   5 nodes to be scored    (282 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (295 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.and.CG.phenotype.BP_weight01_5_all  --- no of nodes:  22 
#> 'data.frame':    0 obs. of  0 variables
#> NULL
#> data frame with 0 columns and 0 rows
#> 
#> Building most specific GOs .....
#>  ( 200 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 380 GO terms and 646 relations. )
#> 
#> Annotating nodes ...............
#>  ( 518 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 96 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  5 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   7 nodes to be scored    (26 eliminated genes)
#> 
#>   Level 8:   11 nodes to be scored   (71 eliminated genes)
#> 
#>   Level 7:   10 nodes to be scored   (98 eliminated genes)
#> 
#>   Level 6:   14 nodes to be scored   (136 eliminated genes)
#> 
#>   Level 5:   15 nodes to be scored   (215 eliminated genes)
#> 
#>   Level 4:   12 nodes to be scored   (396 eliminated genes)
#> 
#>   Level 3:   16 nodes to be scored   (448 eliminated genes)
#> 
#>   Level 2:   2 nodes to be scored    (489 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (509 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.and.CG.phenotype.CC_weight01_5_all  --- no of nodes:  15 
#> 'data.frame':    3 obs. of  2 variables:
#>  $ GO.ID: Factor w/ 3 levels "PF3D7_0110800",..: 1 2 3
#>  $ 1    : chr  "PF3D7_0110800" "PF3D7_0714000" "PF3D7_1211300"
#> NULL
#>           GO.ID             1
#> 1 PF3D7_0110800 PF3D7_0110800
#> 2 PF3D7_0714000 PF3D7_0714000
#> 3 PF3D7_1211300 PF3D7_1211300
#> 
#> Building most specific GOs .....
#>  ( 286 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 558 GO terms and 720 relations. )
#> 
#> Annotating nodes ...............
#>  ( 339 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 117 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 10:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   6 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 8:   6 nodes to be scored    (3 eliminated genes)
#> 
#>   Level 7:   9 nodes to be scored    (29 eliminated genes)
#> 
#>   Level 6:   20 nodes to be scored   (50 eliminated genes)
#> 
#>   Level 5:   22 nodes to be scored   (75 eliminated genes)
#> 
#>   Level 4:   25 nodes to be scored   (129 eliminated genes)
#> 
#>   Level 3:   22 nodes to be scored   (201 eliminated genes)
#> 
#>   Level 2:   5 nodes to be scored    (251 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (324 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.Neutral.MF_weight01_5_all  --- no of nodes:  22 
#> 'data.frame':    2 obs. of  2 variables:
#>  $ GO.ID: Factor w/ 2 levels "PF3D7_1221000",..: 1 2
#>  $ 1    : chr  "PF3D7_1221000" "PF3D7_1322100"
#> NULL
#>           GO.ID             1
#> 1 PF3D7_1221000 PF3D7_1221000
#> 2 PF3D7_1322100 PF3D7_1322100
#> 
#> Building most specific GOs .....
#>  ( 384 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 1342 GO terms and 2675 relations. )
#> 
#> Annotating nodes ...............
#>  ( 347 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 357 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 12:  4 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 11:  11 nodes to be scored   (0 eliminated genes)
#> 
#>   Level 10:  17 nodes to be scored   (25 eliminated genes)
#> 
#>   Level 9:   27 nodes to be scored   (49 eliminated genes)
#> 
#>   Level 8:   34 nodes to be scored   (100 eliminated genes)
#> 
#>   Level 7:   42 nodes to be scored   (137 eliminated genes)
#> 
#>   Level 6:   68 nodes to be scored   (218 eliminated genes)
#> 
#>   Level 5:   69 nodes to be scored   (251 eliminated genes)
#> 
#>   Level 4:   44 nodes to be scored   (294 eliminated genes)
#> 
#>   Level 3:   31 nodes to be scored   (321 eliminated genes)
#> 
#>   Level 2:   9 nodes to be scored    (338 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (340 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.Neutral.BP_weight01_5_all  --- no of nodes:  59 
#> 'data.frame':    17 obs. of  5 variables:
#>  $ GO.ID     : Factor w/ 4 levels "GO:0072594","GO:0020033",..: 1 1 1 1 2 2 2 2 2 2 ...
#>  $ GO:0072594: chr  "PF3D7_0524000" "PF3D7_0808100" "PF3D7_1136400" "PF3D7_1448600" ...
#>  $ GO:0020033: chr  NA NA NA NA ...
#>  $ GO:0032268: chr  NA NA NA NA ...
#>  $ GO:0007034: chr  NA NA NA NA ...
#> NULL
#>         GO.ID    GO:0072594    GO:0020033    GO:0032268    GO:0007034
#> 1  GO:0072594 PF3D7_0524000          <NA>          <NA>          <NA>
#> 2  GO:0072594 PF3D7_0808100          <NA>          <NA>          <NA>
#> 3  GO:0072594 PF3D7_1136400          <NA>          <NA>          <NA>
#> 4  GO:0072594 PF3D7_1448600          <NA>          <NA>          <NA>
#> 5  GO:0020033          <NA> PF3D7_0500100          <NA>          <NA>
#> 6  GO:0020033          <NA> PF3D7_0632800          <NA>          <NA>
#> 7  GO:0020033          <NA> PF3D7_0808700          <NA>          <NA>
#> 8  GO:0020033          <NA> PF3D7_1040800          <NA>          <NA>
#> 9  GO:0020033          <NA> PF3D7_1200600          <NA>          <NA>
#> 10 GO:0020033          <NA> PF3D7_1219300          <NA>          <NA>
#> 11 GO:0020033          <NA> PF3D7_1221000          <NA>          <NA>
#> 12 GO:0020033          <NA> PF3D7_1322100          <NA>          <NA>
#> 13 GO:0032268          <NA>          <NA> PF3D7_0811300          <NA>
#> 14 GO:0032268          <NA>          <NA> PF3D7_0907700          <NA>
#> 15 GO:0032268          <NA>          <NA> PF3D7_1108000          <NA>
#> 16 GO:0007034          <NA>          <NA>          <NA> PF3D7_0808100
#> 17 GO:0007034          <NA>          <NA>          <NA> PF3D7_1448600
#> 
#> Building most specific GOs .....
#>  ( 200 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 380 GO terms and 646 relations. )
#> 
#> Annotating nodes ...............
#>  ( 518 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 110 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  2 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  4 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   7 nodes to be scored    (23 eliminated genes)
#> 
#>   Level 8:   12 nodes to be scored   (82 eliminated genes)
#> 
#>   Level 7:   12 nodes to be scored   (100 eliminated genes)
#> 
#>   Level 6:   18 nodes to be scored   (149 eliminated genes)
#> 
#>   Level 5:   18 nodes to be scored   (234 eliminated genes)
#> 
#>   Level 4:   17 nodes to be scored   (420 eliminated genes)
#> 
#>   Level 3:   17 nodes to be scored   (457 eliminated genes)
#> 
#>   Level 2:   2 nodes to be scored    (500 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (517 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.HS.Neutral.CC_weight01_5_all  --- no of nodes:  38 
#> 'data.frame':    2 obs. of  2 variables:
#>  $ GO.ID: Factor w/ 2 levels "PF3D7_0808100",..: 1 2
#>  $ 1    : chr  "PF3D7_0808100" "PF3D7_1448600"
#> NULL
#>           GO.ID             1
#> 1 PF3D7_0808100 PF3D7_0808100
#> 2 PF3D7_1448600 PF3D7_1448600
#> 
#> Building most specific GOs .....
#>  ( 286 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 558 GO terms and 720 relations. )
#> 
#> Annotating nodes ...............
#>  ( 339 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 180 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  3 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   7 nodes to be scored    (3 eliminated genes)
#> 
#>   Level 8:   9 nodes to be scored    (11 eliminated genes)
#> 
#>   Level 7:   13 nodes to be scored   (36 eliminated genes)
#> 
#>   Level 6:   33 nodes to be scored   (56 eliminated genes)
#> 
#>   Level 5:   39 nodes to be scored   (91 eliminated genes)
#> 
#>   Level 4:   41 nodes to be scored   (165 eliminated genes)
#> 
#>   Level 3:   25 nodes to be scored   (235 eliminated genes)
#> 
#>   Level 2:   8 nodes to be scored    (281 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (335 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.not.classified.MF_weight01_5_all  --- no of nodes:  16 
#> 'data.frame':    0 obs. of  0 variables
#> NULL
#> data frame with 0 columns and 0 rows
#> 
#> Building most specific GOs .....
#>  ( 384 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 1342 GO terms and 2675 relations. )
#> 
#> Annotating nodes ...............
#>  ( 347 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 484 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 14:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 13:  1 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 12:  5 nodes to be scored    (3 eliminated genes)
#> 
#>   Level 11:  15 nodes to be scored   (3 eliminated genes)
#> 
#>   Level 10:  28 nodes to be scored   (28 eliminated genes)
#> 
#>   Level 9:   43 nodes to be scored   (64 eliminated genes)
#> 
#>   Level 8:   56 nodes to be scored   (137 eliminated genes)
#> 
#>   Level 7:   66 nodes to be scored   (175 eliminated genes)
#> 
#>   Level 6:   88 nodes to be scored   (254 eliminated genes)
#> 
#>   Level 5:   82 nodes to be scored   (283 eliminated genes)
#> 
#>   Level 4:   52 nodes to be scored   (314 eliminated genes)
#> 
#>   Level 3:   36 nodes to be scored   (326 eliminated genes)
#> 
#>   Level 2:   10 nodes to be scored   (340 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (346 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.not.classified.BP_weight01_5_all  --- no of nodes:  47 
#> 'data.frame':    0 obs. of  0 variables
#> NULL
#> data frame with 0 columns and 0 rows
#> 
#> Building most specific GOs .....
#>  ( 200 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 380 GO terms and 646 relations. )
#> 
#> Annotating nodes ...............
#>  ( 518 genes annotated to the GO terms. )
#> 
#>           -- Weight01 Algorithm -- 
#> 
#>       the algorithm is scoring 175 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#> 
#>   Level 11:  4 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  9 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 9:   16 nodes to be scored   (29 eliminated genes)
#> 
#>   Level 8:   25 nodes to be scored   (101 eliminated genes)
#> 
#>   Level 7:   23 nodes to be scored   (133 eliminated genes)
#> 
#>   Level 6:   27 nodes to be scored   (200 eliminated genes)
#> 
#>   Level 5:   26 nodes to be scored   (261 eliminated genes)
#> 
#>   Level 4:   20 nodes to be scored   (424 eliminated genes)
#> 
#>   Level 3:   22 nodes to be scored   (458 eliminated genes)
#> 
#>   Level 2:   2 nodes to be scored    (501 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (518 eliminated genes)
#> /Users/Jenna/Library/CloudStorage/Box-Box/Genomics/omicshub/projects/pfGO/Routput/GO/hierarchy.plots/tGO.not.classified.CC_weight01_5_all  --- no of nodes:  25 
#> 'data.frame':    14 obs. of  2 variables:
#>  $ GO.ID: Factor w/ 14 levels "PF3D7_0104100",..: 1 2 3 4 5 6 7 8 9 10 ...
#>  $ 1    : chr  "PF3D7_0104100" "PF3D7_0104200" "PF3D7_0207500" "PF3D7_0207600" ...
#> NULL
#>            GO.ID             1
#> 1  PF3D7_0104100 PF3D7_0104100
#> 2  PF3D7_0104200 PF3D7_0104200
#> 3  PF3D7_0207500 PF3D7_0207500
#> 4  PF3D7_0207600 PF3D7_0207600
#> 5  PF3D7_0511300 PF3D7_0511300
#> 6  PF3D7_0732500 PF3D7_0732500
#> 7  PF3D7_1021800 PF3D7_1021800
#> 8  PF3D7_1041000 PF3D7_1041000
#> 9  PF3D7_1101100 PF3D7_1101100
#> 10 PF3D7_1116800 PF3D7_1116800
#> 11 PF3D7_1135400 PF3D7_1135400
#> 12 PF3D7_1300200 PF3D7_1300200
#> 13 PF3D7_1334800 PF3D7_1334800
#> 14 PF3D7_1418100 PF3D7_1418100
#> 
#> All interesting-gene categories have been tested for GO-term enrichment.
#> 
#> See ALL TOP 30 enriched terms by interesting-gene category in 'Routput/GO/results*.tab.txt' and 'Routput/GO/all.combined.GO.results.tab.txt'.
#> 
#> See log files for topGO-analyses by each interesting-gene category, including all genes in the analysis by GO term in 'Routput/GO/genes_by_GOterm.*.tab.txt'.
#>          GO.ID                                        Term Annotated
#> 1   GO:0042802                   identical protein binding         3
#> 2   GO:0008569 ATP-dependent microtubule motor activity...         3
#> 3   GO:0045505           dynein intermediate chain binding         3
#> 4   GO:0051959     dynein light intermediate chain binding         4
#> 5   GO:0008170                N-methyltransferase activity         4
#> 6   GO:0008173              RNA methyltransferase activity         4
#> 7   GO:0036094                      small molecule binding        18
#> 8   GO:0004540                       ribonuclease activity         5
#> 9   GO:0003924                             GTPase activity         5
#> 10  GO:0005319                  lipid transporter activity         5
#> 11  GO:0010468               regulation of gene expression        27
#> 12  GO:0006401                       RNA catabolic process         8
#> 13  GO:0019682 glyceraldehyde-3-phosphate metabolic pro...         3
#> 14  GO:0060627    regulation of vesicle-mediated transport         3
#> 15  GO:0046394        carboxylic acid biosynthetic process         3
#> 16  GO:0044262     cellular carbohydrate metabolic process         3
#> 17  GO:0019751                    polyol metabolic process         3
#> 18  GO:0007007 inner mitochondrial membrane organizatio...         3
#> 19  GO:0009132    nucleoside diphosphate metabolic process         3
#> 20  GO:0140053               mitochondrial gene expression         3
#> 21  GO:0005930                                     axoneme         3
#> 22  GO:0030286                              dynein complex         3
#> 23  GO:0005768                                    endosome         8
#> 24  GO:0019866                    organelle inner membrane        10
#> 25  GO:0031966                      mitochondrial membrane        10
#> 26  GO:0005856                                cytoskeleton         9
#> 27  GO:0012505                         endomembrane system        51
#> 28  GO:0030139                           endocytic vesicle        15
#> 29  GO:0020011                                  apicoplast        42
#> 30  GO:0016021              integral component of membrane        47
#> 31  GO:0008094               DNA-dependent ATPase activity         7
#> 32  GO:0003677                                 DNA binding        25
#> 33  GO:1990837 sequence-specific double-stranded DNA bi...         6
#> 34  GO:0003690                 double-stranded DNA binding         8
#> 35  GO:0031072                  heat shock protein binding         3
#> 36  GO:0022853 active ion transmembrane transporter act...         3
#> 37  GO:1901505 carbohydrate derivative transmembrane tr...         3
#> 38  GO:0015932 nucleobase-containing compound transmemb...         3
#> 39  GO:0015291 secondary active transmembrane transport...         4
#> 40  GO:0019899                              enzyme binding         4
#> 41  GO:0006457                             protein folding        11
#> 42  GO:0043933 protein-containing complex subunit organ...        31
#> 43  GO:0065004                protein-DNA complex assembly         6
#> 44  GO:0006898               receptor-mediated endocytosis         3
#> 45  GO:0051336            regulation of hydrolase activity         3
#> 46  GO:0006334                         nucleosome assembly         3
#> 47  GO:0045454                      cell redox homeostasis         3
#> 48  GO:0046939                  nucleotide phosphorylation         3
#> 49  GO:0022411              cellular component disassembly         3
#> 50  GO:0051276                     chromosome organization        21
#> 51  GO:0032993                         protein-DNA complex         5
#> 52  GO:0031982                                     vesicle        50
#> 53  GO:0005667             transcription regulator complex         3
#> 54  GO:0000785                                   chromatin         3
#> 55  GO:0042555                                 MCM complex         3
#> 56  GO:0030126                           COPI vesicle coat         3
#> 57  GO:0000139                              Golgi membrane         3
#> 58  GO:0005634                                     nucleus       237
#> 59  GO:0005730                                   nucleolus        16
#> 60  GO:1903561                       extracellular vesicle        17
#> 61  GO:0018024 histone-lysine N-methyltransferase activ...         3
#> 62  GO:0050839              cell adhesion molecule binding        17
#> 63  GO:0016779             nucleotidyltransferase activity         4
#> 64  GO:0005048                     signal sequence binding         4
#> 65  GO:0016853                          isomerase activity         4
#> 66  GO:1901681                     sulfur compound binding         4
#> 67  GO:0042277                             peptide binding         5
#> 68  GO:0004722 protein serine/threonine phosphatase act...         5
#> 69  GO:0003729                                mRNA binding        16
#> 70  GO:0008270                            zinc ion binding         7
#> 71  GO:0072594 establishment of protein localization to...         6
#> 72  GO:0020033                         antigenic variation        27
#> 73  GO:0032268 regulation of cellular protein metabolic...        10
#> 74  GO:0007034                          vacuolar transport         3
#> 75  GO:0006605                           protein targeting         8
#> 76  GO:0020013 modulation by symbiont of host erythrocy...        20
#> 77  GO:0045892 negative regulation of transcription, DN...         4
#> 78  GO:0006892       post-Golgi vesicle-mediated transport         4
#> 79  GO:0034968                  histone lysine methylation         4
#> 80  GO:0006612               protein targeting to membrane         4
#> 81  GO:0010008                           endosome membrane         3
#> 82  GO:0020030             infected host cell surface knob        18
#> 83  GO:0036464       cytoplasmic ribonucleoprotein granule         4
#> 84  GO:0020002                   host cell plasma membrane        28
#> 85  GO:0005654                                 nucleoplasm        11
#> 86  GO:0043657                                   host cell       105
#> 87  GO:0020036                              Maurer's cleft        43
#> 88  GO:0016020                                    membrane       108
#> 89  GO:0005783                       endoplasmic reticulum        31
#> 90  GO:0020008                                     rhoptry         8
#> 91  GO:0004197        cysteine-type endopeptidase activity         6
#> 92  GO:0042393                             histone binding         8
#> 93  GO:0004843 thiol-dependent ubiquitin-specific prote...         5
#> 94  GO:0140096     catalytic activity, acting on a protein        62
#> 95  GO:0016298                             lipase activity         4
#> 96  GO:0008080                N-acetyltransferase activity         4
#> 97  GO:0008047                   enzyme activator activity         7
#> 98  GO:0005524                                 ATP binding        11
#> 99  GO:0003743      translation initiation factor activity         7
#> 100 GO:0003735          structural constituent of ribosome        14
#> 101 GO:0035891                         exit from host cell         6
#> 102 GO:0016579                    protein deubiquitination         5
#> 103 GO:0007165                         signal transduction        23
#> 104 GO:0051603 proteolysis involved in cellular protein...        17
#> 105 GO:0043161 proteasome-mediated ubiquitin-dependent ...         4
#> 106 GO:1903047                  mitotic cell cycle process         9
#> 107 GO:0000724 double-strand break repair via homologou...         8
#> 108 GO:0006412                                 translation        27
#> 109 GO:0020035 cytoadherence to microvasculature, media...        22
#> 110 GO:0046474    glycerophospholipid biosynthetic process         7
#> 111 GO:0020003                 symbiont-containing vacuole        18
#> 112 GO:0020009                                   microneme         6
#> 113 GO:0009986                                cell surface        25
#> 114 GO:0020002                   host cell plasma membrane        28
#> 115 GO:0032991                  protein-containing complex       121
#> 116 GO:0030176 integral component of endoplasmic reticu...         8
#> 117 GO:0071944                              cell periphery        27
#> 118 GO:0044228                           host cell surface         4
#> 119 GO:0098796                    membrane protein complex        20
#> 120 GO:0000151                    ubiquitin ligase complex         3
#>     Significant Expected   topGO go.category   interest.category
#> 1             2     0.46 0.06300          MF      A.HS.Sensitive
#> 2             2     0.46 0.06300          MF      A.HS.Sensitive
#> 3             2     0.46 0.06300          MF      A.HS.Sensitive
#> 4             2     0.61 0.11300          MF      A.HS.Sensitive
#> 5             2     0.61 0.15200          MF      A.HS.Sensitive
#> 6             2     0.61 0.15200          MF      A.HS.Sensitive
#> 7             3     2.76 0.15500          MF      A.HS.Sensitive
#> 8             2     0.77 0.17000          MF      A.HS.Sensitive
#> 9             2     0.77 0.17000          MF      A.HS.Sensitive
#> 10            2     0.77 0.17000          MF      A.HS.Sensitive
#> 11            6     4.36 0.00790          BP      A.HS.Sensitive
#> 12            3     1.29 0.06850          BP      A.HS.Sensitive
#> 13            2     0.48 0.06890          BP      A.HS.Sensitive
#> 14            2     0.48 0.06890          BP      A.HS.Sensitive
#> 15            2     0.48 0.06890          BP      A.HS.Sensitive
#> 16            2     0.48 0.06890          BP      A.HS.Sensitive
#> 17            2     0.48 0.06890          BP      A.HS.Sensitive
#> 18            2     0.48 0.06890          BP      A.HS.Sensitive
#> 19            2     0.48 0.06890          BP      A.HS.Sensitive
#> 20            2     0.48 0.06890          BP      A.HS.Sensitive
#> 21            2     0.48 0.06800          CC      A.HS.Sensitive
#> 22            2     0.48 0.06800          CC      A.HS.Sensitive
#> 23            3     1.28 0.12100          CC      A.HS.Sensitive
#> 24            4     1.60 0.15700          CC      A.HS.Sensitive
#> 25            4     1.60 0.15700          CC      A.HS.Sensitive
#> 26            3     1.44 0.15900          CC      A.HS.Sensitive
#> 27            9     8.17 0.16000          CC      A.HS.Sensitive
#> 28            3     2.40 0.16100          CC      A.HS.Sensitive
#> 29            9     6.73 0.21400          CC      A.HS.Sensitive
#> 30            9     7.53 0.21900          CC      A.HS.Sensitive
#> 31            3     0.25 0.00110          MF HS.and.CG.phenotype
#> 32            6     0.88 0.00120          MF HS.and.CG.phenotype
#> 33            2     0.21 0.03290          MF HS.and.CG.phenotype
#> 34            3     0.28 0.05920          MF HS.and.CG.phenotype
#> 35            1     0.11 0.10280          MF HS.and.CG.phenotype
#> 36            1     0.11 0.10280          MF HS.and.CG.phenotype
#> 37            1     0.11 0.10280          MF HS.and.CG.phenotype
#> 38            1     0.11 0.10280          MF HS.and.CG.phenotype
#> 39            1     0.14 0.13480          MF HS.and.CG.phenotype
#> 40            1     0.14 0.13480          MF HS.and.CG.phenotype
#> 41            2     0.41 0.05900          BP HS.and.CG.phenotype
#> 42            4     1.16 0.06200          BP HS.and.CG.phenotype
#> 43            2     0.22 0.10100          BP HS.and.CG.phenotype
#> 44            1     0.11 0.10900          BP HS.and.CG.phenotype
#> 45            1     0.11 0.10900          BP HS.and.CG.phenotype
#> 46            1     0.11 0.10900          BP HS.and.CG.phenotype
#> 47            1     0.11 0.10900          BP HS.and.CG.phenotype
#> 48            1     0.11 0.10900          BP HS.and.CG.phenotype
#> 49            1     0.11 0.10900          BP HS.and.CG.phenotype
#> 50            3     0.79 0.12700          BP HS.and.CG.phenotype
#> 51            3     0.19 0.00047          CC HS.and.CG.phenotype
#> 52            5     1.93 0.09869          CC HS.and.CG.phenotype
#> 53            1     0.12 0.11162          CC HS.and.CG.phenotype
#> 54            1     0.12 0.11162          CC HS.and.CG.phenotype
#> 55            1     0.12 0.11162          CC HS.and.CG.phenotype
#> 56            1     0.12 0.11162          CC HS.and.CG.phenotype
#> 57            1     0.12 0.11162          CC HS.and.CG.phenotype
#> 58           13     9.15 0.11629          CC HS.and.CG.phenotype
#> 59            2     0.62 0.12292          CC HS.and.CG.phenotype
#> 60            2     0.66 0.13615          CC HS.and.CG.phenotype
#> 61            2     0.41 0.04900          MF          HS.Neutral
#> 62            5     2.31 0.06500          MF          HS.Neutral
#> 63            2     0.54 0.09000          MF          HS.Neutral
#> 64            2     0.54 0.09000          MF          HS.Neutral
#> 65            2     0.54 0.09000          MF          HS.Neutral
#> 66            2     0.54 0.09000          MF          HS.Neutral
#> 67            3     0.68 0.13100          MF          HS.Neutral
#> 68            2     0.68 0.13400          MF          HS.Neutral
#> 69            4     2.17 0.15800          MF          HS.Neutral
#> 70            2     0.95 0.24300          MF          HS.Neutral
#> 71            4     0.67 0.00180          BP          HS.Neutral
#> 72            8     3.03 0.00560          BP          HS.Neutral
#> 73            3     1.12 0.03400          BP          HS.Neutral
#> 74            2     0.34 0.03440          BP          HS.Neutral
#> 75            4     0.90 0.05920          BP          HS.Neutral
#> 76            5     2.25 0.06010          BP          HS.Neutral
#> 77            2     0.45 0.06390          BP          HS.Neutral
#> 78            2     0.45 0.06390          BP          HS.Neutral
#> 79            2     0.45 0.06390          BP          HS.Neutral
#> 80            2     0.45 0.06390          BP          HS.Neutral
#> 81            2     0.39 0.04700          CC          HS.Neutral
#> 82            5     2.36 0.07300          CC          HS.Neutral
#> 83            2     0.53 0.08500          CC          HS.Neutral
#> 84            6     3.68 0.14700          CC          HS.Neutral
#> 85            3     1.44 0.16500          CC          HS.Neutral
#> 86           19    13.78 0.16800          CC          HS.Neutral
#> 87            8     5.64 0.18700          CC          HS.Neutral
#> 88           14    14.18 0.19000          CC          HS.Neutral
#> 89            5     4.07 0.23500          CC          HS.Neutral
#> 90            2     1.05 0.28300          CC          HS.Neutral
#> 91            6     4.05 0.09300          MF      not.classified
#> 92            8     5.40 0.13600          MF      not.classified
#> 93            5     3.38 0.13900          MF      not.classified
#> 94           45    41.88 0.19700          MF      not.classified
#> 95            4     2.70 0.20600          MF      not.classified
#> 96            4     2.70 0.20600          MF      not.classified
#> 97            6     4.73 0.20700          MF      not.classified
#> 98            9     7.43 0.24900          MF      not.classified
#> 99            6     4.73 0.27700          MF      not.classified
#> 100          11     9.46 0.28000          MF      not.classified
#> 101           6     4.13 0.10000          BP      not.classified
#> 102           5     3.44 0.15000          BP      not.classified
#> 103          15    15.84 0.16000          BP      not.classified
#> 104          15    11.71 0.22000          BP      not.classified
#> 105           4     2.76 0.22000          BP      not.classified
#> 106           6     6.20 0.23000          BP      not.classified
#> 107           7     5.51 0.23000          BP      not.classified
#> 108          22    18.60 0.25000          BP      not.classified
#> 109          17    15.15 0.27000          BP      not.classified
#> 110           6     4.82 0.30000          BP      not.classified
#> 111          14    12.06 0.02700          CC      not.classified
#> 112           6     4.02 0.08900          CC      not.classified
#> 113          20    16.75 0.09700          CC      not.classified
#> 114          22    18.76 0.12700          CC      not.classified
#> 115          84    81.06 0.13000          CC      not.classified
#> 116           7     5.36 0.13400          CC      not.classified
#> 117          20    18.09 0.19900          CC      not.classified
#> 118           4     2.68 0.20000          CC      not.classified
#> 119          14    13.40 0.20200          CC      not.classified
#> 120           3     2.01 0.30000          CC      not.classified
```
