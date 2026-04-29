## get.pfannot ----
#' @title
#' Extracts and formats annotations from a gff file from PlasmoDB.
#'
#' @description
#' Extracts and formats annotations from a gff file from PlasmoDB. Not required presently for the GO enrichment pipeline, but provides useful context for results. Opens a connection to the .gff file from PlasmoDB without downloading it, then calls get_annot() to extract and format the annotation.
#'
#' @param gff_url url or filepath to a .gff file. Defaults to \url{https://plasmodb.org/common/downloads/release-68/Pfalciparum3D7/gff/data/PlasmoDB-68_Pfalciparum3D7.gff}
#' @family gene product annotation functions
#'
#' @details # \strong{Notes on gff format}
#' The .gff file should be in tabular format with 9 columns, one for each annotated feature associated with a geneID. No formatting is necessary when using the provided url.
#'
#' an annotation created from PlasmoDB's \emph{P. falciparum} 3D7 .gff file (beta release 69, accessed April 2026) pre-formatted using this function and ready for run.topGO.meta is included in this package [pf.annot].
#'
#' You \emph{will} need to update the gff url accordingly to the latest version when PlasmoDB is updated.
#'
#' @seealso [get.annot()]
#' @export
get.pfannot <-
  function(gff_url = "https://plasmodb.org/common/downloads/release-68/Pfalciparum3D7/gff/data/PlasmoDB-68_Pfalciparum3D7.gff") {
    # make connection to gff file without downloading it, then read it in.
    if (grepl("^https?://", gff_url)) {
      con = gzcon(url(gff_url))
    } else {
      con = gzcon(file(gff_url))
    }
    input = readLines(con)

    x = utils::read.delim(textConnection(input),
                   sep = "\t",
                   comment.char = "#",
                   header = FALSE,
                   stringsAsFactors = TRUE)
    # set column names for standard gff file
    colnames(x) = c("seqid","source","type","geneStart","geneEnd","score","strand","phase","attributes")

    annot = get.annot(x) # then the get.annot function formats the gff file for non-redundancy and readability
    # get.annot keeps columns 1:5,7,9 and output is slightly different. rearrange columns to make it easier to get in bed format should that be needed downstream
    annot = annot[,c(1,4:6,3,7:9,2)]
    # final column order will be c("seqid","geneStart","geneEnd","strand","type","geneID","description","geneName","source")
    annot
  }

## including some call to Rgraphviz to nullify build-error
ignore_unused_imports <- function() {
  Rgraphviz::aaa_fun
}
