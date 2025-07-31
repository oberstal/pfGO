# 2025/03/03
# Author: J. Oberstaller

## make_enrichRes----
#' @title make_enrichRes
#' @description
#' Function to turn output of [reformat_sigGenes()] into an enrichResult object, which can then be visualized in many different ways with functions from the enrichplot package.
#'
#'[See an enrichplot visualization overview here](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html).
#' @param my.cpInput  a reformatted dataframe from pfGO output file all.combined.sig.genes.per.sig.term.tsv (reformatted using the [reformat_sigGenes()] function)
#' @param interestCategory  character string of your category of interest (should match a value from the interest.category column of my.cpInput)
#' @param ontology  one of "BP", "CC" or "MF"
#' @param pval_cutoff significance cutoff. Default is 0.05.
#' @returns An enrichResult object (which can be visualized in various ways using functions from the enrichplot package). The enrichResult object class is described in the DOSE package.
#' @examples
#' # example pfGO-output-to-visualization workflow
#' \dontrun{
#' my.sigGenes <- read.delim("Routput/GO/all.combined.sig.genes.per.sig.term.tsv")
#' my.cp.df <- reformat_sigGenes(my.sigGenes)
#' my.enrichRes <- make_enrichRes(my.cp.df, interestCategory = "sensitive", ontology = "BP")
#'
#' ## from this point, enrichment results can be visualized many different ways with functions included in the enrichplot package.
#'
#' barplot(my.enrichRes)
#' cnetplot(my.enrichRes)
#' dotplot(my.enrichRes)
#'
#' }
#' @export
#'
make_enrichRes <- function(my.cpInput, interestCategory, ontology = "BP", pval_cutoff = 0.05){

  #require(DOSE)
  #require(enrichplot)

  res_df = my.cpInput[my.cpInput$interest.category==interestCategory & my.cpInput$go.category==ontology,]

  # make vector out of geneIDs (actually geneNames) for gene universe
  geneNames = unlist(stringr::str_split(res_df$geneID, pattern = "/"), use.names = FALSE)

  ## create enrichResult object (enrichResult object is a class defined from 'DOSE' package)
  enrichRes = methods::new("enrichResult",
                           readable = TRUE,
                           result = res_df,
                           pvalueCutoff = pval_cutoff,
                           qvalueCutoff = 1,
                           pAdjustMethod = "none",
                           organism = "P.falciparum",
                           ontology = ontology,
                           keytype = "UNKNOWN",
                           gene = res_df$geneID,
                           universe = unique(geneNames[my.cpInput$interest.category!=interestCategory])
  )

  # fill in term similarity slot in enrichRes object (required for cnetplots)
  termSim = enrichplot::pairwise_termsim(enrichRes)

  # return the termSim enrichRes object
  termSim

}
