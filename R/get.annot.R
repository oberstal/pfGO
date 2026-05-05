## get.annot ----
#' @title
#' Extracts and formats annotations from a .gff file.
#'
#' @description
#' Extracts and formats all gene annotations (biotypes "protein_coding_gene", ncRNA_gene", and "pseudogene") from a .gff file. Not required presently for the GO enrichment pipeline, but provides useful context for results. Can be used as-is with a provided gff file as input, or is called by get_pfannot to get the gff file from PlasmoDB.
#'
#' @param x input gff file
#' @family gene product annotation functions


#' *notes on gff format*
#' The gff file should be in tabular format with 9 columns, one for each annotated feature associated with a geneID. No formatting is necessary when using the file at the provided url.
#'
#' an annotation created from PlasmoDB's latest \emph{P. falciparum} .gff file (from PlasmoDB v66; accessed September 28, 2023) pre-formatted using this function is included in this package (pf.annot).
#'
#' @seealso [formatGOdb()]
#' @export
get.annot <- function(x) {
  # "type" column is standard third column from .gff file
  # keep only entries for "protein_coding_gene","ncRNA_gene" and "pseudogene" types
  x = x %>% dplyr::filter(type == "protein_coding_gene" | type == "ncRNA_gene" | type == "pseudogene")
  # Keep only informative columns
  x = x[, c(1:5, 7, 9)]

  # format attributes column to get annotations
  ## geneID will be first field (semicolon-delimited) in attributes column. remove all characters except the geneID
  x$geneID = stringr::str_replace(x$attributes, ";.+", "")
  x$geneID = stringr::str_replace(x$geneID, "ID=", "")

  # description is always last or next-to last field in attributes column
  x$description = stringr::str_replace(x$attributes, ".+description=", "")
  # remove any trailing fields
  x$description = stringr::str_replace(x$description, ";.+", "")
  # and replace the "%2C" ascii formatting with comma
  x$description = stringr::str_replace_all(x$description, "\\%2C", ",")


  # some entries have "Name=" entries; others only have "description=". use names for the genes that have a genesymbol, and for the others set name = to what's left (the geneID)
  x$geneName = "blank"
  # remove all characters leading up to the gene symbol
  x$geneName = stringr::str_replace(x$attributes, ".+Name=", "")
  # then remove all characters after the gene symbol (will start with semicolon)
  x$geneName = stringr::str_replace(x$geneName, ";.+", "")
  x$geneName = stringr::str_replace(x$geneName, "ID=", "")
  # now drop "attributes" column
  x$attributes = NULL
  x
}
