## code to prepare `DATASET` dataset goes here

usethis::use_data(geneID2GO, overwrite = TRUE)
# curated P. falciparum GOdb (accessed from GeneDB Dec 08, 2020)
input <- "data-raw/Pf_curatedDec082020_GOdb.txt"
geneID2GO <- readMappings(input)
Pfal_geneID2GO <- geneID2GO

# example *piggyBac* mutant classification-data from pooled 1k heat shock screen (to use as "mydf" parameter for run.topGO.meta())
usethis::use_data(mydf, overwrite = TRUE)
input <- "data-raw/run.topGO.input.df.txt"
mydf <- read.delim(input, header = TRUE)
