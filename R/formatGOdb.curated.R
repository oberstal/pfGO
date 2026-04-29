## formatGOdb.curated----
#' @title
#' Generates new GO database from curated evidence codes only for functional enrichment using run.topGO.meta.
#'
#' @description
#' Generates new GO database from curated evidence codes only for functional enrichment using run.topGO.meta.
#'
#' @param url url or filepath to a .gaf or gz-compressed gaf file. This should be the complete GO.gaf file (not one that has "Curated" in the PlasmoDB filename).
#'
#' @param organism optional string to include in your output filename. Defaults to "PF3D7".
#' @param plasmoDB_version optional string to include the PlasmoDB release version in your output filename (e.g. "v68").
#'
#' @family formatGOdb functions
#'
#' @details
#' Outputs a dated file named "*_GOdb.out" to the ./Routput folder. ./Routput will be created if it doesn't already exist.
#'
#' You will need to run `topGO::readmappings()` on the file generated with the formatGOdb functions to use it as the geneID2GO parameter.
#'
#' Example: geneID2GO <- topGO::readmappings("./Routput/PF3D7_v68_20260422_GOdb.out")
#'
#' @details # \strong{Notes on gaf.gz format}
#' The .gaf or .gaf.gz file should be in tabular format with 17 columns, one row for each GO term associated with a geneID. No formatting is necessary when using the provided url. Note that zip archives (extension: .zip) won't work directly with this function. Zip files need to be downloaded and unzipped, then the 'url' argument should be updated to the local unzipped file path.
#'
#' A GOdb from Plasmo DB's \emph{P. falciparum} annotation (from PlasmoDB beta release 69; accessed April 2026) pre-formatted using this function and ready for run.topGO.meta is included in this package ([Pfal_geneID2GO_curated]) and is used as the default geneID2GO input to run.topGO.meta(). 
#'
#' @seealso [formatGOdb()]
#' @export
formatGOdb.curated <-
  function(url,
           organism = "PF3D7",
           plasmoDB_version = ""){

    # old plasmoDB gaf link structure (up until release 65): "https://plasmodb.org/common/downloads/release-56/Pfalciparum3D7/gaf/PlasmoDB-56_Pfalciparum3D7_GO.gaf"
      ## NOTE: seems plasmoDB is not following a consistent convention yet after release 64; v65 uses zip archives while v66 switches to gzip files (up until v65 files were not compressed). The method here will work to extract gzip and noncompressed files. Zip archives are more complicated and won't work with this function. Zip files need to be downloaded and unzipped, then update the 'url' argument to the local unzipped file path.

    # make connection to .gaf or .gaf.gz file without downloading it, then read it in.
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

    # weed out any electronically-inferred evidence-codes (only want curated GO terms)
    x = x[x[,7]!="IEA",]
    ## NOTE: don't need to filter manually anymore, PlasmoDB offers their own curated gaf file with IEA,NAS and IBA codes filtered

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
                      "_curated_",
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
