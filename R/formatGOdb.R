## formatGOdb ----
#' @title
#' Creates a GO database compatible with topGO from .gaf file retrieved from PlasmoDB.
#'
#' @description
#' Fetches \emph{P. falciparum} annotations from PlasmoDB, from which it creates a GO database compatible with topGO. The output can be used as the geneID2GO parameter for run.topGO.meta. A generally useful tool for keeping GO analyses up-to-date.
#'
#' @param url  url or filepath to a .gaf or gz-compressed gaf file. E.g. current latest \emph{Pf} consortium .gaf annotation file hosted at PlasmoDB: \url{https://plasmodb.org/common/downloads/release-68/Pfalciparum3D7/gaf/PlasmoDB-68_Pfalciparum3D7_GO.gaf.gz}.
#' Note that zip archives (extension: .zip) won't work directly with this function. Zip files need to be downloaded and unzipped, then the 'url' argument should be updated to the local unzipped file path.
#' @param organism  optional string to include in your output file-name. Defaults to "PF3D7".
#' @param plasmoDB_version  optional string to include the PlasmoDB release version in your output filename (e.g. "v68").
#' @family formatGOdb functions
#'
#' @details
#'Outputs a dated file named "*_GOdb.out" to the ./Routput folder. ./Routput will be created if it doesn't already exist. Note that the PlasmoDB release version isn't recorded automatically in the filename, so it is good practice to include the current PlasmoDB release number in the filename via the 'organism' parameter for your records.
#'
#' You will need to run `topGO::readmappings()` on the file generated with formatGOdb before using it as the geneID2GO parameter.
#'
#' Example: geneID2GO <- topGO::readmappings("./Routput/PF3D7_v68_20260422_GOdb.out")
#'
#'
#' @details # \strong{Notes on .gaf format}
#' The .gaf file should be in tabular format with 17 columns, one row for each GO term associated with a geneID. No formatting is necessary when using the provided url.
#'
#' Retrieves GO annotations assigned by \strong{all} evidence-codes. 
#'
#' a GOdb from PlasmoDB's \emph{P. falciparum} annotation (from PlasmoDB beta release 69; accessed April 2026) pre-formatted using the this function and ready for run.topGO.meta is included in this package [Pfal_geneID2GO].
#' 
#' A version of this function that weeds out any non-curated, Inferred-from-Electronic-Annotation assignments, Inferred from Biological aspect of Ancestor, or Nontraceable Author Statements is also included in this package (evidence codes IEA,IBA,NAS; see \code{\link{formatGOdb.curated}}).
#'
#' @seealso [formatGOdb.curated()]
#'
#' @export
formatGOdb <- function(url,
           organism = "PF3D7",
           plasmoDB_version = "") {
  # updating to require manual input of URL as PlasmoDB link structure is not yet consistent. Example that will work: https://plasmodb.org/common/downloads/release-68/Pfalciparum3D7/gaf/PlasmoDB-68_Pfalciparum3D7_GO.gaf.gz

  # old pre-v65 plasmoDB link: "https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gaf/PlasmoDB-56_Pfalciparum3D7_GO.gaf"
    # make connection to gaf file without downloading it, then read it in.
    if (grepl("^https?://", url)) {
      con = gzcon(url(url))
    } else {
      con = gzcon(file(url))
    }
    input = readLines(con)
    x = utils::read.delim(textConnection(input),
                   sep = "\t",
                   comment.char = "!",
                   header = FALSE,
                   stringsAsFactors = TRUE)
    # Keep only the 2nd and 5th columns (geneID, GO)
    x = x[, c(2, 5)]

    # keep only unique rows (in original file--a geneID will be listed multiple times with the same GO.ID if there are different lines of evidence supporting it)
    x = unique.data.frame(x)

    # replace all the "/.[0-9]" in the geneID column so no transcript-names
    x[,1] = gsub(pattern = "\\.[0-9]",
                       replacement = "",
                       x[,1],
                       perl = TRUE)

    # keep only unique rows again because now multi-transcript genes will have multiple entries named the same thing, though they don't add GO info (go for one isoform is same for another)
    x = unique.data.frame(x)

    # get unique geneIDs
    id_list = as.vector(unique(x[, 1]))
#    print(head(id_list))
#    print(str(id_list))

    # get date for output filename
    today = Sys.Date() # e.g., "2026-04-22"
    today = format(today, format = "%Y%m%d") # e.g., "20260422"
#    today = format(today, format = "%b%d%Y")

    # create output directory if it doesn't exist
    dir.create("Routput", showWarnings = FALSE)

    # create output file
    GOdb.file = paste("Routput/",
                      organism,
                      "_",
                      plasmoDB_version,
                      "_",
                      today,
                      "_GOdb.out",
                      sep = "")

    cat("creating output file", GOdb.file, "(this may take a couple minutes) . . .\n\n")

    # ensure we overwrite any existing file
    if (file.exists(GOdb.file)) file.remove(GOdb.file)

    for (i in id_list) {
      cat(paste(i),
          "\t",
          sep = "",
          file = GOdb.file,
          append = TRUE)

      GO = droplevels(x[which(x[, 1] == i), 2])

      for (g in GO) {
        # if g is not the last GO term for the geneID, paste GO.ID followed by a comma. Otherwise paste GO.ID alone (no trailing commas)

        if (g != GO[length(GO)]) {
          cat(paste(g),
              ",",
              sep = "",
              file = GOdb.file,
              append = TRUE)
        } else {
          cat(paste(g),
              "\n",
              sep = "",
              file = GOdb.file,
              append = TRUE)
        }
      }
    }
    GO.db = utils::read.delim(GOdb.file, header = FALSE)
    cat("File is complete. To create the geneID2GO object for running a GO enrichment, run 'topGO::readMappings(\"", GOdb.file, "\")'\n\n", sep = "")
    return(invisible(GO.db))
  }

