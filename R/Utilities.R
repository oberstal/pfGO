library(topGO)
#library(ggplot2) # loads with tidyverse
library(plyr) # does not load with tidyverse
library(dplyr) # use tidyverse instead?
library(scales)
library(formatR)
library(stringr) # loads with tidyverse

#' Creating directory-structure
#'
#' Creates the output-directory structure I use for my topGO pipeline. All the included makeDir functions evaluate to the newly created path, or to the existing path if it already exists.
#'
#' @param ... not required. Defaults to creating top-level Routput directory in current working directory.
#'
#' @name makeDirs
#' @export
NULL
#> NULL

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


#' @rdname makeDirs
#' @keywords internal
#' @export
makeGOsig.genes.dir <- function() {
  file.path = getwd()
  mainDir = "/Routput/GO"
  subDir = "/sig.genes.by.term"
  newDir = paste(file.path, mainDir, subDir, sep = "")
  dir.create(newDir, showWarnings = FALSE)
  newDir
}

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








############### run.topGO.meta ----
# script slightly modified from my run.topGO script especially to test for GO enrichment in art-R meta-analysis gene-categories of interest. STILL CLUNKY BUT IT WORKS

#' @title
#' Run JO's topGO pipeline
#'
#' @description
#' This function tests for functional enrichment in gene-categories of interest.
#'
#' @param mydf data frame with geneIDs in column 1, and interest-category classifications in column 2
#' @param geneID2GO A list of named vectors of GO IDs--one vector of GO-terms for each geneID.
#'
#' @section **outputs**:
#' run.topGO.meta creates several output-files, including:
#' * enrichment results,
#' * significant genes per significant term,
#' * plots of the GO-term hierarchy relevant to the analysis, and
#' * thorough log-files for each gene-category of interest tested against the background of all other genes in the analysis.
#'
#' Primary results from run.topGO.meta will be in "Routput/GO/all.combined.GO.results.tab.txt".
#'
#'
#' @details
#' The **run.topGO.meta** function:
#' * defines which genes are "interesting" and which should be defined as background for each category specified in mydf,
#' * makes the GOdata object for topGO,
#' * tests each category of interest for enriched GO-terms against all the other genes included in mydf (the "gene universe"),
#' * and then outputs results to several tables (tab.txt files that can be opened in Excel).
#'
#' Enrichments are performed by each ontology (molecular function, biological process, cellular compartment) sequentially on all groups of interest. Results are combined in the final output-table ("Routput/GO/all.combined.GO.results.tab.txt").

#'
#'
#' TopGO automatically accounts for genes that cannot be mapped to GO terms (or are mapped to terms with < 3 genes in the analysis) with "feasible genes" indicated in the topGO.log files in the "Routput/GO" folder.
#'
#'
#' @section **Concepts for common use-cases**:
#'
#' *RNAseq*:
#' In an RNAseq analysis, common categories might be "upregulated", "downregulated", and "neutral". The gene universe would consist of all genes detected above your threshold cutoffs (*not necessarily all genes in the genome*).
#'
#' *piggyBac screens*:
#' In pooled *piggyBac*-mutant screening, common categories might be "sensitive", "tolerant", and "neutral". The gene universe would consist of all genes represented in your screened library of mutants (*again, not all genes in the genome*).
### See the included [mydf] as an example. (doesn't work, add in later)
#'
#'
#' @section **Using your own custom GO database**:
#'
#' A correctly formatted geneID2GO object is included for P. falciparum enrichment analyses ([Pfal_geneID2GO]). You may also provide your own, so long as it is a named character-vector of GO-terms (each vector named by geneID, with GO terms as each element).
#'
#' You can use the included [formatGOdb.curated()] function to format a custom GO database from curated GeneDB annotations for several non-model organisms (or the [formatGOdb()] function to include all GO annotations, if you aren't picky about quality of automated electronic annotations). If you're studying a model organism, several annotations are already available through the AnnotationDbi bioconductor package that loads with topGO.
#'
#' @seealso [topGO::topGO()]
#'
#' @examples
#' run.topGO.meta(mydf,Pfal_geneID2GO)
#'
#' @export
run.topGO.meta <- function(mydf = "mydf", geneID2GO = "geneID2GO", pval = 0.05) {
  require(topGO)
  require(tidyverse)
  require(plyr)

  # make required directories for output if they don't exist. each one evaluates to that path, so can save the paths as variables
  R.dir = makeRoutput.dir()
  GO.dir = makeGOoutput.dir()
  sigGO.dir = makeGOsig.genes.dir()
  hier.dir = makeGOhierarchy.dir()

  # go ahead and write the input-df to file in the run-folder so you don't have to wonder later:
  write.table(
    mydf,
    file = "Routput/GO/run.TopGO.input.df.txt",
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )


  # test enrichment for EACH interesting-genes category against all the other background genes
  geneids = mydf[, 1]
  interesting_genes_raw = as.factor(mydf[, 2])
  interesting_genes = levels(interesting_genes_raw)

  interesting.category.counter = 0
  for (i in interesting_genes) {
    interesting.category.counter = interesting.category.counter + 1

    logfile = paste(GO.dir,"/topGO.log.", i, ".txt", sep = "")
    cat(
      paste(
        "***Genes by GO-term for art-R meta-analysis interest-category ",
        i,
        "***\n",
        sep = ""
      ),
      file = logfile,
      append = TRUE
    )

    myInterestingGenes = mydf[which(mydf[, 2] == i), 1]

    # make a named factor classifying each gene as interesting or not based on myInterestingGenes (0 = not interesting, 1 = interesting)
    geneList = factor(as.integer(geneids %in% myInterestingGenes))
    names(geneList) = geneids

    ## iterate through each essentiality category by ontology ##
    ontology = c("MF", "BP", "CC")
    ontology.counter = 0
    for (o in ontology) {
      # print some notes to logfile
      cat(
        paste(
          "\n\n\n====================================================================================\n"
        ),
        file = logfile,
        append = TRUE
      )
      cat(
        paste("Gene-category of interest: ", i, sep = ""),
        file = logfile,
        append = TRUE
      )
      cat(paste("\nOntology: ", o, sep = ""),
          file = logfile,
          append = TRUE)
      cat(
        paste(
          "\n====================================================================================\n"
        ),
        file = logfile,
        append = TRUE
      )
      cat(
        paste("\n*** Begin topGO results summary ***\n"),
        file = logfile,
        append = TRUE
      )

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
      capture.output(GOdata, file = logfile, append = TRUE)

      # Enrichment analyses: could use any number of statistical tests; the "weight01" algorithm is the default

      ### resultFisher <- runTest(GOdata,algorithm="classic", statistic="Fisher")
      resultTopgo = topGO::runTest(GOdata, algorithm = "weight01", statistic =
                              "Fisher")
      ### resultElim <- runTest(GOdata,algorithm="elim", statistic="Fisher")

      # output genes in significant GO terms to logfile
      sig.genes = topGO::sigGenes(GOdata)
      cat(
        paste("All genes in significant GO terms in gene universe:\n\n\n"),
        file = logfile,
        append = TRUE
      )
      capture.output(sig.genes, file = logfile, append = TRUE)

      cat(
        paste(
          "\n==============================================================================\n"
        ),
        file = logfile,
        append = TRUE
      )

      # generate a plot of the GO hierarchy highlighting the most significant terms (can help demonstrate how terms were collapsed) and output to .pdf
      topGO::printGraph(
        GOdata,
        resultTopgo,
        firstSigNodes = 5,
        fn.prefix = paste(hier.dir,
                          "/tGO.",
                          i,
                          ".",
                          o,
                          sep = ""),
        useInfo = "all",
        pdfSW = TRUE
      )

      ## make a results table with ALL enriched GO terms
      res = topGO::GenTable(
        GOdata,
        topGO = resultTopgo,
        orderBy = "topGO",
        ranksOf = "fisher"
        #topNodes = 50
      )

      # add columns to results-table for go-category and interest category
      res$go.category = o
      res$interest.category = i
      # change any p.values that are way too low (close to 0--looks like cutoff for p-val assignment is 1e-30) and have the "<" symbol which makes this column be considered as a factor, not as numbers, to a minimum value so they will be kept for FDR-adjustment:

      res$topGO[which(res$topGO == "<1e-30")] <-
        "1e-30"
      # next convert the topGO column to numeric (this way no "NA's" will be introduced)
      res$topGO = as.numeric(res$topGO)

      # do a p. adjust for FDR on the entire dataset (NOT NECESSARY for topGO method):
      #        res$FDR = p.adjust(res$topGO, method = "fdr")
      res.significant = res[which(res$topGO <= pval),]
      #                                    & res$FDR <= 0.05),]

      # # to output ONLY the significant genes from enriched GO terms, not every gene in a significant GO term in the gene universe(this seems to be the part where things go wrong and pipeline breaks--when it comes to converting to a df):
      # goresults.genes works when tested step by step
      goresults.genes = sapply(res.significant$GO.ID, function(x) {
        # get all genes for each GO term (will be list of geneIDs named for GO id)
        genes = topGO::genesInTerm(GOdata, x)
        # for every gene in the GO-id list,
        sapply(genes, function(x) {
          genes.in.term.test = x[x %in% sig.genes]
          genes.in.term.test
        } )
      })

      cat(
        paste(
          "\n*** SIGNIFICANT genes in significant GO terms: ***\n\n\n"
        ),
        file = logfile,
        append = TRUE
      )
      capture.output(goresults.genes, file = logfile, append = TRUE)
      cat(
        paste(
          "\n\n--------------------------------------------------------------------------------\n\n\n"
        ),
        file = logfile,
        append = TRUE
      )

      genes.in.term.lists = sapply(goresults.genes,function(x){
        new.row = unlist(x, recursive = FALSE, use.names = TRUE)
      })

      # for testing
      genes.in.terms.df = plyr::ldply(genes.in.term.lists, rbind, .id = "GO.ID")
      # current output is is a df with unnecessary empty NA values in columns to preserve spacing (prob from writing from a list of lists to a table row by row?? )
        # output should ideally just be a long list--column 1 = GO ids, column 2 = geneID mapped to that term. One geneID/GO pair per row. Then for each GO term, add column for term definition, and for each geneID, add column for annotation.

      # THIS MAGICAL ONE-LINER TURNS THE WHOLE LIST OF LISTS OF UNEVEN SIZES INTO A BEAUTIFUL DF (this one works)
#      genes.in.terms.df = plyr::ldply(genes.in.term.lists, rbind)
      ## 4/26/2021: this is where I would fix the format of the output genes-in-terms file
      print(str(genes.in.terms.df))
      print(genes.in.terms.df)


      # # goresults.genes is a list of lists, one for each GO term. Reformat to df and write to file:

      cat(
        paste(
          "\n*** Significant genes in significant GO-terms (table output) ***\n\n\n"
        ),
        file = logfile,
        append = TRUE
      )
      capture.output(genes.in.terms.df, file = logfile, append = TRUE)
      ontology.counter = ontology.counter + 1
      # print status messages to logfile
      cat(
        paste(
          "\nontology-category counter is ",
          ontology.counter,
          " of 3",
          sep = ""
        ),
        file = logfile,
        append = TRUE
      )
      cat(
        paste(
          "\ninteresting-category counter is ",
          interesting.category.counter,
          " of ",
          length(interesting_genes),
          sep = ""
        ),
        file = logfile,
        append = TRUE
      )

      # build on to the results-table for each ontology (3 loops for each interest category)

      if (ontology.counter != 1) {
        combined.GO.output = rbind.data.frame(combined.GO.output, res)
        combined.significant.GO.output = rbind.data.frame(combined.significant.GO.output, res.significant)
        #        combined.sig.per.term.output = rbind.data.frame(combined.sig.per.term.output, genes.in.terms.df)
        # try using bind_rows instead to get around problem of joining dataframes with uneven number of columns. IT WORKS!
        combined.sig.per.term.output = dplyr::bind_rows(combined.sig.per.term.output, genes.in.terms.df)

      } else {
        combined.GO.output = res
        combined.significant.GO.output = res.significant
        combined.sig.per.term.output = genes.in.terms.df
      }
    } ### ONTOLOGY LOOP ENDS HERE

    # print status-messages to indicate end of interest-category log-file
    cat(
      paste(
        "\n--------------------------------------------------------------------------------\n"
      ),
      file = logfile,
      append = TRUE
    )
    cat(
      paste(
        "************************************************************************************\n"
      ),
      file = logfile,
      append = TRUE
    )
    cat(
      paste("***END category ", i, " enrichment results***\n"),
      file = logfile,
      append = TRUE
    )
    cat(
      paste(
        "************************************************************************************\n"
      ),
      file = logfile,
      append = TRUE
    )

    ##### output a table with all GO (MF, BP, CC) significant-genes-per-term analyses in one for each interest category
    write.table(
      combined.sig.per.term.output,
      file = paste(
        "Routput/GO/sig.genes.by.term/",
        i,
        "_sig.genes.per.term.txt",
        sep = ""
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    #    build on to the results-tables for each interest-category (1 loop for each of the interest categories)
    if (interesting.category.counter != 1) {
      all.bin.combined.GO.output = rbind.data.frame(all.bin.combined.GO.output, combined.GO.output)

    } else {
      all.bin.combined.GO.output = combined.GO.output
    }
  }### THIS ENDS THE LOOP FOR EACH INTERESTING-GENE CATEGORY (counter incremented at beginning of loop)

  # output a table with ALL interest-category GO analyses in one with an added column for interest-category
  write.table(
    all.bin.combined.GO.output,
    paste("Routput/GO/all.combined.GO.results.tab.txt",
          sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  # print some progress-messages to screen
  cat("\nAll interesting-gene categories have been tested for GO-term enrichment.")
  cat(
    "\n\nSee SIGNIFICANT enrichment by interesting-gene category in 'Routput/GO/significant.GO.results*.tab.txt' and 'Routput/GO/all.combined.significant.GO.results.tab.txt'."
  )
  cat(
    "\n\nSee ALL TOP 30 enrichment by interesting-gene category in 'Routput/GO/results*.tab.txt' and 'Routput/GO/all.combined.GO.results.tab.txt'."
  )
  cat(
    "\n\nSee log files for topGO-analyses by each interesting-gene category, including all genes in the analysis by GO term in 'Routput/GO/genes_by_GOterm.*.tab.txt'."
  )

  # sig genes per term doesn't work with creating runGOdata object yet. But option for the future to store multiple outputs from single run. returns major output as R object (does not interfere with created output files)
#  setClass("runGOdata", slots = c(enrichmentResult = "ANY", sigGenes = "ANY"))

#  myrunGOdata <- new("runGOdata", enrichmentResult = all.bin.combined.GO.output, sigGenes = combined.sig.per.term.output)

#  myrunGOdata
  all.bin.combined.GO.output
}

#' @title get.value
#'
#'@description
#' a lookup function to match an ID and return the matching values
#' @param id a character-string you want to match on
#' @param lookupvector a named vector you will search for your id
#'
#' @details
#' This one needs more documentation. I don't remember what it returns as I wrote it a long time ago. Could be a vector of matching values, or could be a vector of positions you can then filter the lookupvector on to return the values.
#'
#' @export
get.value <- function(id, lookupvector = named.vector){
  value = unname(lookupvector[id])
  return(value)
}


#' @title
#' formatGOdb
#'
#' @description
#' This function fetches P. falciparum annotations from GeneDB, from which it creates a GO database compatible with topGO. The output can be used as the geneID2GO parameter for run.topGO.meta. A generally useful tool for keeping GO analyses up-to-date.
#'
#' @param gaf.gz_url = url to a gaf.gz file. Defaults to latest Pf consortium .gaf annotation file from GeneDB, hosted here: ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/Pfalciparum/Pfalciparum.gaf.gz
#' @param organism = optional string to include in your output file-name. Defaults to "Pf".
#'
#' @details
#'
#' You may need to run `topGO::readmappings()` on the file generated with formatGOdb before using it as the geneID2GO parameter.
#'
#' *notes on gaf.gz format*
#' The gaf.gz file should be in tabular format with 17 columns, one row for each GO term associated with a geneID. No formatting is necessary when using the provided url.
#'
#' Retrieves GO annotations assigned by *all* evidence-codes. A version that weeds out any non-curated, inferred-from-electronic-annotation assignments is also included in this package (evidence code IEA; see format.curated.GOdb))
#'
#' a GOdb from Sanger's latest P. falciparum annotation (accessed December 8, 2020) pre-formatted using the *curated* version of this function and ready for run.topGO.meta is included in this package (Pfal_geneID2GO).
#'
#' @seealso [formatGOdb.curated()]
#'
#' @export
formatGOdb <- function(gaf.gz_url = "ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/Pfalciparum/Pfalciparum.gaf.gz",
           organism = "Pf") {
    # make connection to gaf.gz file without downloading it, then read it in.
    con = gzcon(url(gaf.gz_url))
    input = readLines(con)
    x = read.delim(textConnection(input),
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
    today = Sys.Date()
    today = format(today, format = "%b%d%Y")

    # create output file
    GOdb.file = paste("Routput/", organism, "_", today, "_GOdb.out", sep = "")

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
    GO.db = read.delim(GOdb.file, header = FALSE)
    return(GO.db)
  }

#' @title
#' formatGOdb.curated
#'
#' @description
#' Generates new GO database from curated evidence-codes only for functional enrichment using run.topGO.meta.
#'
#' @param gaf.gz_url url to a gaf.gz file. Defaults to latest Pf consortium .gaf annotation file from GeneDB.
#'@param organism optional string to include in your output file-name. Defaults to "Pf".
#'
#' @details
#' You will need to run `topGO::readmappings()` on the file generated with formatGOdb before using it as the geneID2GO parameter.
#'
#' *notes on gaf.gz format*
#' The gaf.gz file should be in tabular format with 17 columns, one row for each GO term associated with a geneID. No formatting is necessary when using the provided url.
#'
#' a GOdb from Sanger's latest P. falciparum annotation (accessed December 8, 2020) pre-formatted using this function and ready for run.topGO.meta is included in this package (Pfal_geneID2GO).
#'
#' @seealso [formatGOdb()]
#' @export
formatGOdb.curated <-
  function(gaf.gz_url = "https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gaf/PlasmoDB-54_Pfalciparum3D7_GO.gaf",
           organism = "Pf"){

#  function(gaf.gz_url = "ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/Pfalciparum/Pfalciparum.gaf.gz",
#           organism = "Pf") {
    # make connection to gaf.gz file without downloading it, then read it in.

    con = gzcon(url(gaf.gz_url))

    input = readLines(con)
    x = read.delim(textConnection(input),
                   sep = "\t",
                   comment.char = "!",
                   header = FALSE,
                   stringsAsFactors = TRUE)

    # weed out any electronically-inferred evidence-codes (only want curated GO terms)
    x = x[x[,7]!="IEA",]

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
    today = Sys.Date()
    today = format(today, format = "%b%d%Y")

    # create output file
    GOdb.file = paste("Routput/", organism, "_curated", today, "_GOdb.out", sep = "")

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
    GO.db = read.delim(GOdb.file, header = FALSE)
    return(GO.db)
    }


#' @title
#' get_annot
#'
#' @description
#' Extracts and formats annotations from a gff file. Not required presently for the GO enrichment pipeline, but provides useful context for results. Can be used as-is with a provided gff file as input, or is called by get_pfannot to get the gff file from plasmoDB.
#'
#' @param x input gff file


#' *notes on gff format*
#' The gff file should be in tabular format with 9 columns, one for each annotated feature associated with a geneID. No formatting is necessary when using the file at the provided url.
#'
#' an annotation created from PlasmoDB's latest P. falciparum gff file (accessed November 1, 2021) pre-formatted using this function and ready for run.topGO.meta is included in this package (pf.annot).
#'
#' @seealso [formatGOdb()]
#' @export
get_annot <- function(x) {
  require(tidyverse)
  # keep only entries for "protein_coding_gene" or "ncRNA_gene" types
  x = x %>% filter(type == "protein_coding_gene" | type == "ncRNA_gene")
  # Keep only informative columns
  x = x[, c(1:5, 7, 9)]

  # format attributes column to get annotations
  ## geneID will be first field (semicolon-delimited) in attributes column. remove all characters except the geneID
  x$geneID = stringr::str_replace(x$attributes, ";.+", "")
  x$geneID = stringr::str_replace(x$geneID, "ID=", "")

  # description is always last field in attributes column
  x$description = stringr::str_replace(x$attributes, ".+description=", "")
  # and replace the "%2C" misformattings with comma
  x$description = stringr::str_replace(x$description, "\\%2C", ",")

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


#' @title
#' get_pfannot
#'
#' @description
#' Extracts and formats annotations from a gff file from plasmoDB. Not required presently for the GO enrichment pipeline, but provides useful context for results. Opens a connection to the gff file from plasmoDB without downloading it, then calls get_annot() to extract and format the annotation.
#'
#' @param gff_url connection to gff file. Defaults to "https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gff/data/PlasmoDB-54_Pfalciparum3D7.gff"


#' *notes on gff format*
#' The gff file should be in tabular format with 9 columns, one for each annotated feature associated with a geneID. No formatting is necessary when using the provided url.
#'
#' an annotation created from PlasmoDB's latest P. falciparum gff file (accessed November 1, 2021) pre-formatted using this function and ready for run.topGO.meta is included in this package (pf.annot).
#'
#' @seealso [get_annot()]
#' @export
get_pfannot <-
  function(gff_url = "https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gff/data/PlasmoDB-54_Pfalciparum3D7.gff",
           organism = "Pf") {
    # make connection to gff file without downloading it, then read it in.
    con = gzcon(url(gff_url))
    input = readLines(con)

    x = read.delim(textConnection(input),
                   sep = "\t",
                   comment.char = "#",
                   header = FALSE,
                   stringsAsFactors = TRUE)
    # set column names for standard gff file
    colnames(x) = c("seqid","source","type","feature_start","feature_end","score","strand","phase","attributes")

    annot = get_annot(x) # then the get_annot function formats the gff file for non-redundancy and readability
    annot
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

