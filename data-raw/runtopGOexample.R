# curated P. falciparum GOdb (accessed from GeneDB Oct 28/Nov 5, 2021)
#input <- "data-raw/Pf_curatedOct282021_GOdb.out"
#Pfal_geneID2GO <- readMappings(input)


# load included Pfal geneID2GO object
data("Pfal_geneID2GO")
# read in example *piggyBac* mutant classification-data from pooled 1k heat shock screen (to use as "mydf" parameter for run.topGO.meta())
input <- "data-raw/run.topGO.input.df.txt"
mydf <- read.delim(input, header = TRUE)
# run topGO pipeline
run.topGO.meta(mydf,Pfal_geneID2GO)
