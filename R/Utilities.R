#' Creating directory structure
#'
#' Creates an output-directory structure for the pfGO topGO pipeline. All the included makeDir functions evaluate to the newly created path, or to the existing path if it already exists.
#'
#' param ... not required. Defaults to creating top-level Routput directory in current working directory.
#'
#' @name makeDirs
#' @noRd
#' @keywords internal
NULL
#> NULL

#' @noRd
#' @rdname makeDirs
#' @keywords internal
#' @export
makeRoutput.dir <- function () {
  mainDir = getwd()
  subDir = "/Routput"
  newDir = paste(mainDir,
                 subDir,
                 sep = "")
  dir.create(newDir,
             showWarnings = FALSE)
  newDir
}

#' @noRd
#' @rdname makeDirs
#' @keywords internal
#' @export
makeGOoutput.dir <- function() {
  file.path = getwd()
  mainDir = "/Routput/"
  subDir = "/GO"
  newDir = paste(file.path,
                 mainDir,
                 subDir,
                 sep = "")
  dir.create(newDir, showWarnings = FALSE)
  newDir
}


#' @noRd
#' @rdname makeDirs
#' @keywords internal
#' @export
makeGOhierarchy.dir <- function() {
  file.path = getwd()
  mainDir = "/Routput/GO"
  subDir = "/hierarchy.plots"
  newDir = paste(file.path,
                 mainDir,
                 subDir,
                 sep = "")
  dir.create(newDir, showWarnings = FALSE)
  newDir
}



## get.GOdef -----
#' @title
#' Get GO terms and definitions
#'
#' @description
#' Extracts the term, term-definition, and ontology for all GOIDs (e.g. GO:0000027) in a geneID2GO object into a dataframe. Output also includes a column for all geneIDs mapping to each term as a comma-separated list. Not required presently for the GO enrichment pipeline, but provides useful context for results--very handy to use as a lookup-table for both GO terms and individual geneIDs.
#'
#'
#' @param geneID2GO  geneID2GO object. Defaults to Pfal_geneID2GO_curated (if using the default, you must first load the data-object \link{Pfal_geneID2GO_curated} ).
#'
#' @export
#'
#' @examples
#' # example code
#' data(Pfal_geneID2GO_curated)  ## load Pfal_geneID2GO_curated object
#' go2genesLookup.df <- get.GOdef(Pfal_geneID2GO_curated)
#'
get.GOdef <- function(geneID2GO = Pfal_geneID2GO_curated){
  # get go2genes
  go2genes = topGO::inverseList(geneID2GO)
  go2genes.df = data.frame(sapply(go2genes,stringr::str_flatten_comma))
  go2genes.df$GOID = rownames(go2genes.df)
  rownames(go2genes.df) = NULL
  colnames(go2genes.df) = c("geneIDsList","GOID")
  # grab all GO terms in GOdb
  terms = names(go2genes)
  # extract all relevant term definitions from GO.db package
  GOdef.df = AnnotationDbi::select(x = GO.db::GO.db, keys = terms, columns = as.character(AnnotationDbi::columns(GO.db::GO.db)), keytype = "GOID")
  # join GO term definitions with all genes mapping to each term (geneIDs in a comma-separated string for each GOterm)
  GOdef.df = dplyr::left_join(GOdef.df,go2genes.df)
  GOdef.df
}






# ########## Other useful functions #############
# ## to.pdf: Plot the output of any function to pdf ##
# # usage: to.pdf(expr = plotFPKM(clusteredbyFPKM), filename = paste("PfPbTg_FPKM_by_binID_and_cluster.pdf", sep = ""), width = 10, height = 4 )
# # need to add a "print" line to the end of the function whose output you want to print if it doesn't already have one, or .pdf will be blank
#
# to.pdf <- function(expr, filename = "Rfigs/plot.pdf", width = 10, height = 4, ..., verbose=TRUE) {
#   if ( verbose )
#     cat(sprintf("Creating %s\n", filename))
#   pdf(filename, ...)
#   on.exit(dev.off())
#   eval.parent(substitute(expr))
# }

