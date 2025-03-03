## reformat_sigGenes----
#' @title
#' Reformat sig genes output for clusterProfiler
#' @description
#' Reformats significant genes in significant terms output to something more compatible with clusterProfiler/enrichPlot for visualization. Output is a dataframe with the appropriate columns that can be turned into an enrichRes object (see the DOSE package), which can then be plotted with enrichPlot functions.
#' @param my.sigGenes   a dataframe of pfGO output all.combined.sig.genes.per.sig.term.tsv
#'
#' @export
#'
#' @examples
#' # example code
#'
#'\dontrun{
#' my.sigGenes <- read.delim("Routput/GO/all.combined.sig.genes.per.sig.term.tsv")
#' my.cp.df <- reformat_sigGenes(my.sigGenes)
#' }
#'
reformat_sigGenes <- function(my.sigGenes){
  my.cpInput <- my.sigGenes %>%
    dplyr::group_by(interest.category,
             GO.ID,
             Annotated,
             Significant,
             topGO) %>%
    dplyr::reframe(ID=Term,
              Description = Term,
              GeneRatio = paste(Significant,
                                "/",
                                Annotated,
                                sep = ""),
              BgRatio = paste(as.character(Annotated-Significant),
                              "/",
                              Annotated),
              pvalue = topGO,
              p.adjust = topGO,
              geneID = paste(geneName,
                             sep = "",
                             collapse = "/"),
              Count = Significant,
              go.category = go.category) %>%
    dplyr::distinct()

  my.cpInput <- as.data.frame(my.cpInput)

  # remove old columns
  my.cpInput <- my.cpInput[,-c(3:5)]
  # add column for qvalue
  my.cpInput$qvalue <- qvalue::qvalue(my.cpInput$pvalue, fdr.level = NULL, lambda = 0)$qvalues
  # add rownames
  rownames(my.cpInput) <- my.cpInput$ID

  my.cpInput
}
