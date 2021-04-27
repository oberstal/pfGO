## 4/26/2021: editing GO pipeline to make output more efficient

# read in test mydf for running the function
input <- "Rdata/run.TopGO.input.df. TEST.txt"
mydf.test <- read.delim(input, header = TRUE)

# example output for troubleshooting too
input <- "Rdata/all.combined.GO.results.tab.txt"
res <- read.delim(input, header = TRUE)

# read in GO stuff
input <- "/Users/Jenna/Box/Genomics/Research/for_piggyBac_mysql/tables/GeneDB/GO/Pf_curatedDec082020_GOdb.out"
geneID2GO <- readMappings(input)

mydf <- mydf.test

# test enrichment for EACH interesting-genes category against all the other background genes
geneids = mydf[, 1]
interesting_genes_raw = as.factor(mydf[, 2])
interest.categories = levels(interesting_genes_raw)


geneList.out <- sapply(interest.categories, get.interestingGenes, mydf = mydf, simplify = FALSE)
  # using sapply to keep the names attribute, which we need to keep the lists named after the interest category. "simplify=FALSE" is required so output doesn't get turned into a matrix (gets kept as a list of lists,binary scores for each geneID by interest category)

test.sig.genes <- lapply(geneList.out, get.resbyOntology, o = "BP")



# try vectorizing instead of looping by interest.category over mydf (so the vector would be interesting_genes)
get.interestingGenes <- function(interest.category, mydf = "mydf"){
  myInterestingGenes = mydf[mydf[,2] == interest.category,1]
  # make a named factor classifying each gene as interesting or not based on myInterestingGenes (0 = not interesting, 1 = interesting)
  geneList = factor(as.integer(geneids %in% myInterestingGenes))
  names(geneList) = geneids
  geneList
}


test.sig.genes <- lapply(geneList.out, get.resbyOntology, o = "BP")
  # siggenes is just "significant" geneIDs by interest-category--does not provide the enriched term that made them significant

test.out <- lapply(geneList.out, get.resbyOntology2, o = "BP")

test.bp.out2 <- lapply(geneList.out, get.resbyOntology3, o = "BP")

# getting closer to the format i want the sig genes in terms in
  #unlist2 is a function from annotationdbi that "doesn't mangle the names"
unlist.test <- unlist2(test.bp.out2$UP_NF_only$sig.genes.by.sig.term)
unlist.test2 <- unlist2(test.bp.out2[2])

#alright time to forget all this. Optimize this stupid set of functions later. Just use what works and see if I can get a basic package made ffs!!!

GOdata <- new(
  "topGOdata",
  ontology = "BP",
  allGenes = geneList.out[[7]],
  annot = annFUN.gene2GO,
  gene2GO = geneID2GO,
  nodeSize = 3,
  description = 'GO analysis of genes comprising each art-R interest-category against all other genes in the comparison'
)



ont = c("BP", "MF", "CC")


# vectorize the ontology part over the output from get.interestingGenes
get.resbyOntology3 <- function(geneList="geneList", o = "BP"){

  GOdata = new(
    "topGOdata",
    ontology = o,
    allGenes = geneList,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO,
    nodeSize = 3,
    description = 'GO analysis of genes comprising each art-R interest-category against all other genes in the comparison'
  )

  # Enrichment analyses: could use any number of statistical tests; the "weight01" algorithm is the default

  ### resultFisher <- runTest(GOdata,algorithm="classic", statistic="Fisher")
  resultTopgo = runTest(GOdata, algorithm = "weight01", statistic =
                          "Fisher")
  ### resultElim <- runTest(GOdata,algorithm="elim", statistic="Fisher")

  # output genes in significant GO terms
  sig.genes = sigGenes(GOdata)
  # what does sig.genes look like? a list of sig genes named after each interest-category

  ## make a results table with ALL enriched GO terms
  res = GenTable(
    GOdata,
    topGO = resultTopgo,
    orderBy = "topGO",
    ranksOf = "fisher"
    #topNodes = 50
  )

  # add columns to results-table for go-category and interest category
  res$go.category = o
#  res$interest.category = i

  res$topGO[which(res$topGO == "<1e-30")] <-
    "1e-30"
  # next convert the topGO column to numeric (this way no "NA's" will be introduced)
  res$topGO = as.numeric(res$topGO)
  res.significant = res[which(res$topGO <= 0.1),]


  # goresults.genes works when tested step by step below ----
  goresults.genes = sapply(res.significant$GO.ID, function(x) {
    # get all genes for each GO term (will be list of geneIDs named for GO id)
    genes = genesInTerm(GOdata, x)
    # for every gene in the GO-id list,
    sapply(genes, function(x) {
      genes.in.term.test = x[x %in% sig.genes]
      genes.in.term.test
    } )
  })
  # to return multiple outputs (enrichment results, then significant genes, for each interest category)
  outputs = list("enrichment.results" = res, "significant.genes" = sig.genes, "sig.genes.by.sig.term" = goresults.genes)

  return(outputs)
}




# original crappy working version below ----
interesting.category.counter = 0
for (i in interesting_genes) {
  interesting.category.counter = interesting.category.counter + 1

  myInterestingGenes = mydf[which(mydf[, 2] == i), 1]

  # make a named factor classifying each gene as interesting or not based on myInterestingGenes (0 = not interesting, 1 = interesting)
  geneList = factor(as.integer(geneids %in% myInterestingGenes))
  names(geneList) = geneids

  ## iterate through each essentiality category by ontology ##
  ontology = c("MF", "BP", "CC")
  ontology.counter = 0
  for (o in ontology) {
  #make GOdata object if you're using a custom GO database (ie, Plasmodium)
  GOdata = new(
    "topGOdata",
    ontology = o,
    allGenes = geneList,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO,
    nodeSize = 3,
    description = 'GO analysis of genes comprising each art-R interest-category against all other genes in the comparison'
  )

  # Enrichment analyses: could use any number of statistical tests; the "weight01" algorithm is the default

  ### resultFisher <- runTest(GOdata,algorithm="classic", statistic="Fisher")
  resultTopgo = runTest(GOdata, algorithm = "weight01", statistic =
                          "Fisher")
  ### resultElim <- runTest(GOdata,algorithm="elim", statistic="Fisher")

  # output genes in significant GO terms
  sig.genes = sigGenes(GOdata)

  # generate a plot of the GO hierarchy highlighting the most significant terms (can help demonstrate how terms were collapsed) and output to .pdf
  printGraph(
    GOdata,
    resultTopgo,
    firstSigNodes = 5,
    fn.prefix = paste("Routput/GO/hierarchy.plots/tGO", i, o, sep = "."),
    useInfo = "all",
    pdfSW = TRUE
  )

  ## make a results table with ALL enriched GO terms
  res = GenTable(
    GOdata,
    topGO = resultTopgo,
    orderBy = "topGO",
    ranksOf = "fisher"
    #topNodes = 50
  )

  # add columns to results-table for go-category and interest category
  res$go.category = o
  res$interest.category = i

  res$topGO[which(res$topGO == "<1e-30")] <-
    "1e-30"
  # next convert the topGO column to numeric (this way no "NA's" will be introduced)
  res$topGO = as.numeric(res$topGO)
  res.significant = res[which(res$topGO <= 0.1),]


  if (ontology.counter != 1) {
    combined.GO.output = rbind.data.frame(combined.GO.output, res)
    combined.significant.GO.output = rbind.data.frame(combined.significant.GO.output, res.significant)
    # try using bind_rows instead to get around problem of joining dataframes with uneven number of columns. IT WORKS!
    combined.sig.per.term.output = dplyr::bind_rows(combined.sig.per.term.output, genes.in.terms.df)

  } else {
    combined.GO.output = res
    combined.significant.GO.output = res.significant
    combined.sig.per.term.output = genes.in.terms.df
  }
  } ### ONTOLOGY LOOP ENDS HERE




# # to output ONLY the significant genes from enriched GO terms, not every gene in a significant GO term in the gene universe(this seems to be the part where things go wrong and pipeline breaks--when it comes to converting to a df):

## testing goresults.genes ----
# steps:  (1) get all genes for each GO term (will be list of geneIDs named for GO id); (2) # for every gene in the GO-id list, pull out the ones that are significant


 = sapply(res.significant$GO.ID, ...?)

get.goresults.genes <- function(GOdata,x){
  # get all genes for each GO term (will be list of geneIDs named for GO id)
  genes = genesInTerm(GOdata, x)
  # for every gene in the GO-id list, pull out the ones that are significant
  genes.in.term.test = sapply(genes, get.genesInTerm)
}

genes = sapply(GOdata, genesInTerm)

  # get all genes for each GO term (will be list of geneIDs named for GO id)

  genes = genesInTerm(GOdata, x)
  # for every gene in the GO-id list,
  goresults.genes = get.genesInTerm(x)
  goresults.genes
  }


goresults.genes = sapply(res.significant$GO.ID, function(x) {

  # get all genes for each GO term (will be list of geneIDs named for GO id)
  genes = genesInTerm(GOdata, x)
  # for every gene in the GO-id list, pull out the ones that are significant
  get.genesInTerm <- function(x) {
    genes.in.term.test = x[x %in% sig.genes]
    genes.in.term.test
  }
  genes.in.term.test = sapply(genes, get.genesInTerm)
  } )
  # can I convert genes.in.term.test from

})


# goresults.genes works when tested step by step below ----
goresults.genes = sapply(res.significant$GO.ID, function(x) {
  # get all genes for each GO term (will be list of geneIDs named for GO id)
  genes = genesInTerm(GOdata, x)
  # for every gene in the GO-id list,
  sapply(genes, function(x) {
    genes.in.term.test = x[x %in% sig.genes]
    genes.in.term.test
  } )
})
