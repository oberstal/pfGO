## code to prepare `mydf` dataset goes here
#library(read.table)

# example *piggyBac* mutant classification-data from pooled 1k heat shock screen (to use as "mydf" parameter for run.topGO.meta())
input <- "data-raw/run.topGO.input.df.txt"
exampleMydf <- read.delim(input, header = TRUE)
usethis::use_data(exampleMydf, overwrite = TRUE)
