#library(topGO)
#library(ggplot2) # loads with tidyverse
#library(plyr) # does not load with tidyverse
#library(dplyr) # use tidyverse instead?
#library(scales)
#library(formatR)
#library(stringr) # loads with tidyverse

#' Creating directory-structure
#'
#' Creates the output-directory structure I use for my topGO pipeline. All the included makeDir functions evaluate to the newly created path, or to the existing path if it already exists.
#'
#' @param ... not required. Defaults to creating top-level Routput directory in current working directory.
#'
#' @name makeDirs
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
#' Tests for functional enrichment in gene-categories of interest.
#'
#' @param mydf data frame with geneIDs in column 1, and interest-category classifications in column 2.
#' @param geneID2GO A list of named vectors of GO IDs--one vector of GO-terms for each geneID.
#' @param pval pvalue threshold for significance. Defaults to 0.05.
#'
#' @section **outputs**:
#' run.topGO.meta creates several output-files, including:
#' * enrichment results,
#' * significant genes per significant term,
#' * plots of the GO-term hierarchy relevant to the analysis, and
#' * thorough log-files for each gene-category of interest tested against the background of all other genes in the analysis.
#'
#' Primary results from run.topGO.meta will be in "Routput/GO/all.combined.GO.results.tsv".
#'
#'
#' @details
#' The **run.topGO.meta** function:
#' * defines which genes are "interesting" and which should be defined as background for each category specified in mydf,
#' * makes the GOdata object for topGO,
#' * tests each category of interest for enriched GO-terms against all the other genes included in mydf (the "gene universe"),
#' * and then outputs results to table (.tsv files that can be opened in Excel).
#'
#' Enrichments are performed by each ontology (molecular function, biological process, cellular compartment) sequentially on all groups of interest. Results are combined in the final output-table ("Routput/GO/all.combined.GO.results.tsv").

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
#' See the included [exampleMydf] as an example.
#'
#'
#' @section **Using your own custom GO database**:
#'
#' A correctly formatted geneID2GO object is included for P. falciparum enrichment analyses ([Pfal_geneID2GO]). You may also provide your own, so long as it is a named character-vector of GO-terms (each vector named by geneID, with GO terms as each element).
#'
#' You can use the included [formatGOdb.curated()] function to format a custom GO database from curated GeneDB/PlasmoDB annotations for several non-model organisms (or the [formatGOdb()] function to include all GO annotations, if you aren't picky about quality of automated electronic annotations). If you're studying a model organism, several annotations are already available and can be downloaded through the AnnotationDbi bioconductor package that loads with topGO.
#'
#' @seealso [topGO::topGO()]
#'
#' @examples
#'
#' run.topGO.meta(exampleMydf,Pfal_geneID2GO)
#'
#' @export
run.topGO.meta <- function(mydf = "mydf", geneID2GO = "Pfal_geneID2GO", pval = 0.05) {
#  require(topGO)
#  require(tidyverse)
 # require(plyr)

  # make required directories for output if they don't exist. each one evaluates to that path, so can save the paths as variables
  R.dir = makeRoutput.dir()
  GO.dir = makeGOoutput.dir()
  hier.dir = makeGOhierarchy.dir()

  # go ahead and write the input-df to file in the run-folder so you don't have to wonder later:
  utils::write.table(
    mydf,
    file = "Routput/GO/run.TopGO.input.df.txt",
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )

  # create logfile with system info, all loaded packages and version numbers for reproducibility:
  syslog = "Routput/syslog.txt"
  utils::capture.output(utils::sessionInfo(), file = syslog)


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
        "***Genes by GO-term per interest-category ",
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
          "\n\n\n====================================================================================\n",
          "Gene-category of interest: ",
          i,
          "\nOntology: ",
          o,
          "\n====================================================================================\n",
          "\n*** Begin topGO results summary ***\n",
          sep = ""
        ),
        file = logfile,
        append = TRUE
      )

      #make GOdata object if you're using a custom GO database (ie, Plasmodium)
      GOdata = methods::new(
        "topGOdata",
        ontology = o,
        allGenes = geneList,
        annot = annFUN.gene2GO,
        gene2GO = geneID2GO,
        nodeSize = 3,
        description = 'GO analysis of genes comprising each interest-category against all other genes in the comparison'
      )
      utils::capture.output(GOdata, file = logfile, append = TRUE)

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
      utils::capture.output(sig.genes, file = logfile, append = TRUE)

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

      res$go.category = o
      res$interest.category = i

      #### print results for go category to screen
      cat(
        paste(
          "\n==============================================================================\n",
          "GO enrichment results for interest category: ",
          i,
          "\nOntology: ",
          o,
          sep = ""
          )
      )

      print(res)

      cat(
        paste(
          "\n==============================================================================\n",
          "\ninterest-category ",
          interesting.category.counter,
          " of ",
          length(interesting_genes),
          "\nontology ",
          ontology.counter +1,
          " of 3\n",
          sep = ""
          )
      )
      ####


      # do a p. adjust for FDR on the entire dataset (NOT NECESSARY for topGO method):
      #        res$FDR = p.adjust(res$topGO, method = "fdr")
      res.significant = res[which(res$topGO <= pval),]
      #                                    & res$FDR <= 0.05),]



      # # to output ONLY the significant genes from enriched GO terms, not every gene in a significant GO term in the gene universe:
      goresults.genes = sapply(res.significant$GO.ID, function(x) {
        # get all genes for each GO term (will be list of geneIDs named for GO id)
        genes = topGO::genesInTerm(GOdata, x)
        # for every gene in the GO-id list,
        sapply(genes, function(x) {
          genes.in.term.test = x[x %in% sig.genes]
          genes.in.term.test
        }, simplify = FALSE )
      }, simplify = FALSE ) # do not simplify as some will be returned as vector and others as matrix

      ## grab gene annotations to add to output tables later
      pf.annot = get.pfannot()

      if(length(goresults.genes)>0){
      genes.in.term.lists = AnnotationDbi::unlist2(goresults.genes, use.names = TRUE)
      genes.in.terms.df = data.frame(GO.ID=names(genes.in.term.lists),geneID=genes.in.term.lists)
      genes.in.terms.df$go.category = o
      genes.in.terms.df$interest.category = i
      genes.in.terms.df = suppressMessages(dplyr::left_join(genes.in.terms.df,
                                           pf.annot,
                                           na_matches ="never"))
      }

      if(length(goresults.genes)==0){
        genes.in.term.lists = NULL
        genes.in.terms.df = NULL
      }

        # output (genes.in.terms.df) is a long df--column 1 = GO ids, column 2 = geneID mapped to that term. One geneID/GO pair per row. Then for each GO term, add column for term definition, and for each geneID, add column for annotation.

      ## write sig genes in sig terms df to logfile:

      cat(
        paste(
          "\n==============================================================================\n",
          "\n*** SIGNIFICANT genes in significant (i.e. enriched) GO-terms ***\n",
          "pval parameter: ", pval,
          "\n\n"
        ),
        file = logfile,
        append = TRUE
      )
      utils::capture.output(genes.in.terms.df, file = logfile, append = TRUE)
      ontology.counter = ontology.counter + 1
      # print status messages to logfile
      cat(
        paste(
          "\nontology-category counter is ",
          ontology.counter,
          " of 3",
          "\ninteresting-category counter is ",
          interesting.category.counter,
          " of ",
          length(interesting_genes),
          "\n",
          sep = ""
        ),
        file = logfile,
        append = TRUE
      )

      # build on to the results-table for each ontology (3 loops for each interest category)

      if (ontology.counter != 1) {
        combined.GO.output = rbind.data.frame(combined.GO.output, res)
        combined.significant.GO.output = rbind.data.frame(combined.significant.GO.output, res.significant)
        combined.sig.per.term.output = rbind.data.frame(combined.sig.per.term.output, genes.in.terms.df)
#        sig.per.term.output = genes.in.terms.df

      } else {
        combined.GO.output = res
        combined.significant.GO.output = res.significant
        combined.sig.per.term.output = genes.in.terms.df
#        sig.per.term.output = genes.in.terms.df
      }

        } ### ONTOLOGY LOOP ENDS HERE

    # print status-messages to indicate end of interest-category log-file
    cat(
      paste(
        "\n--------------------------------------------------------------------------------\n",
        "************************************************************************************\n",
        "***END category ", i, " enrichment results***\n",
        "************************************************************************************\n",
        "\n\nSee TOP 30 enriched terms for all interesting-gene categories in 'Routput/GO/all.combined.GO.results.tsv'.",
        "\n\nSee all significant genes mapped to all enriched GO terms in 'Routput/GO/all.combined.sig.genes.per.sig.terms.tsv'.\n",
        sep = ""
        ),
      file = logfile,
      append = TRUE
      )

    #    build on to the results-tables and sig genes per sig term table for each interest-category (1 loop for each of the interest categories)
    if (interesting.category.counter != 1) {
      all.bin.combined.GO.output = rbind.data.frame(all.bin.combined.GO.output, combined.GO.output)
      all.bin.combined.sig.per.term.output = rbind.data.frame(all.bin.combined.sig.per.term.output, combined.sig.per.term.output)

    } else {
      all.bin.combined.GO.output = combined.GO.output
      all.bin.combined.sig.per.term.output = combined.sig.per.term.output
    }
  }### THIS ENDS THE LOOP FOR EACH INTERESTING-GENE CATEGORY (counter incremented at beginning of loop)

  ## sort final output table by interest category, go category, and pvalue (from most to least significant)
  all.bin.combined.GO.output = all.bin.combined.GO.output %>% dplyr::arrange(interest.category,go.category,topGO)

  # output a table with ALL interest-category GO analyses in one with an added column for interest-category
  utils::write.table(
    all.bin.combined.GO.output,
    paste("Routput/GO/all.combined.GO.results.tsv",
          sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  # and output a table with ALL interest-category sig genes per sig term results in one with gene annotations and GOterm annotations:
  all.bin.combined.sig.per.term.output = dplyr::right_join(all.bin.combined.GO.output,all.bin.combined.sig.per.term.output,na_matches ="never")

  ## also sort final output table by interest category, go category, GO.ID, and pvalue (from most to least significant)
  all.bin.combined.sig.per.term.output = all.bin.combined.sig.per.term.output %>% dplyr::arrange(interest.category,go.category,GO.ID,topGO)

  utils::write.table(
    all.bin.combined.sig.per.term.output,
    paste("Routput/GO/all.combined.sig.genes.per.sig.terms.tsv",
          sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  # print final GO results table to screen
  cat("\n==============================================================================\n")
  cat("Final GO enrichment results table (saved to 'Routput/GO/all.combined.GO.results.tsv':\n")
  print(all.bin.combined.GO.output)

  # print some progress-messages to screen
  cat("\n==============================================================================\n\n\n")
  cat("\nAll interesting-gene categories have been tested for GO-term enrichment.")
  cat(
    "\n\nSee ALL TOP 30 enriched terms by interesting-gene category in 'Routput/GO/all.combined.GO.results.tsv'."
  )
  cat(
    "\n\nSee log files for topGO-analyses by each interesting-gene category, including all genes in the analysis by GO term in 'Routput/GO/topGO.log.*.txt'."
  )
  cat(
    "\n\nSee all significant genes mapped to all significant GO terms in 'Routput/GO/all.combined.sig.genes.per.sig.terms.tsv'."
  )
  cat("\n\n==============================================================================\n\n")
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
get.value <- function(id, lookupvector){
  value = unname(lookupvector[id])
  return(value)
}

## formatGOdb ----
#' @title
#' Creates a GO database compatible with topGO from .gaf file retrieved from PlasmoDB.
#'
#' @description
#' Fetches \emph{P. falciparum} annotations from PlasmoDB, from which it creates a GO database compatible with topGO. The output can be used as the geneID2GO parameter for run.topGO.meta. A generally useful tool for keeping GO analyses up-to-date.
#'
#' @param url = url or filepath to a .gaf or gz-compressed gaf file. Defaults to latest \emph{Pf} consortium .gaf annotation file hosted at PlasmoDB: /url{https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gaf/PlasmoDB-CURRENT_Pfalciparum3D7_GO.gaf.gzip}. Note that zip archives won't work directly with this function. Zip files need to be downloaded and unzipped, then the 'url' argument should be updated to the local unzipped file path.
#' @param organism = optional string to include in your output file-name. Defaults to "Pf".
#' @param plasmoDB_version = optional string to include the PlasmoDB release version in your output filename (e.g. "v66").
#'
#' @details
#'Outputs a dated file named "*_GOdb.out" to the ./Routput folder. ./Routput will be created if it doesn't already exist. Note that the PlasmoDB release version isn't recorded automatically in the filename, so it is good practice to include the current PlasmoDB release number in the filename via the 'organism' parameter for your records.
#'
#' You will need to run `topGO::readmappings()` on the file generated with formatGOdb before using it as the geneID2GO parameter.
#'
#' Example: geneID2GO <- topGO::readmappings("./Routput/Pf_Mar022022_GOdb.out")
#'
#'
#' @details # \strong{Notes on .gaf format}
#' The .gaf file should be in tabular format with 17 columns, one row for each GO term associated with a geneID. No formatting is necessary when using the provided url.
#'
#' Retrieves GO annotations assigned by \strong{all} evidence-codes. A version that weeds out any non-curated, Inferred-from-Electronic-Annotation assignments, Inferred from Biological aspect of Ancestor, or Nontraceable Author Statements is also included in this package (evidence codes IEA,IBA,NAS; see \code{\link{formatGOdb.curated}})
#'
#' a GOdb from PlasmoDB's latest \emph{P. falciparum} annotation (accessed November 5, 2021) pre-formatted using the \emph{curated} version of this function and ready for run.topGO.meta is included in this package (Pfal_geneID2GO).
#'
#' @seealso [formatGOdb.curated()]
#'
#' @export
formatGOdb <- function(url = "https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gaf/PlasmoDB-CURRENT_Pfalciparum3D7_GO.gaf.gzip",
           organism = "Pf",
           plasmoDB_version = "") {

  # old pre-v65 plasmoDB link: "https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gaf/PlasmoDB-56_Pfalciparum3D7_GO.gaf"
    # make connection to gaf file without downloading it, then read it in.
    con = gzcon(url(url))
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
    today = Sys.Date()
    today = format(today, format = "%b%d%Y")

    # create output file
    GOdb.file = paste("Routput/",
                      organism,
                      "_",
                      plasmoDB_version,
                      "_",
                      today,
                      "_GOdb.out",
                      sep = "")

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
    return(GO.db)
  }

## formatGOdb.curated----
#' @title
#' Generates new GO database from curated evidence-codes only for functional enrichment using run.topGO.meta.
#'
#' @description
#' Generates new GO database from curated evidence-codes only for functional enrichment using run.topGO.meta.
#'
#' @param url url or filepath to a .gaf or gz-compressed gaf file. Defaults to latest \emph{Pf} consortium .gaf.gzip annotation file from PlasmoDB. Note that zip archives won't work directly with this function. Zip files need to be downloaded and unzipped, then the 'url' argument should be updated to the local unzipped file path.
#'@param organism optional string to include in your output filename. Defaults to "Pf".
#'@param plasmoDB_version = optional string to include the PlasmoDB release version in your output filename (e.g. "v66").
#'
#' @details
#' Outputs a dated file named "*_GOdb.out" to the ./Routput folder. ./Routput will be created if it doesn't already exist.
#'
#' You will need to run `topGO::readmappings()` on the file generated with formatGOdb to use it as the geneID2GO parameter.
#'
#' Example: geneID2GO <- topGO::readmappings("./Routput/Pf_Mar022022_GOdb.out")
#'
#' @details # \strong{Notes on gaf.gz format}
#' The .gaf or .gaf.gz file should be in tabular format with 17 columns, one row for each GO term associated with a geneID. No formatting is necessary when using the provided url.
#'
#' A GOdb from Plasmo DB's latest \emph{P. falciparum} annotation (accessed Nov 5, 2021) pre-formatted using this function and ready for run.topGO.meta is included in this package (Pfal_geneID2GO).
#'
#' @seealso [formatGOdb()]
#' @export
formatGOdb.curated <-
  function(url = "https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gaf/PlasmoDB-CURRENT_Pfalciparum3D7_Curated_GO.gaf.gzip",
           organism = "Pf",
           plasmoDB_version = ""){

    # old plasmoDB gaf link structure (up until release 65): "https://plasmodb.org/common/downloads/release-56/Pfalciparum3D7/gaf/PlasmoDB-56_Pfalciparum3D7_GO.gaf"
      ## NOTE: seems plasmoDB is not following a consistent convention yet after release 64; v65 uses zip archives while v66 switches to gzip files (up until v65 files were not compressed). The method here will work to extract gzip and noncompressed files. Zip archives are more complicated and won't work with this function. Zip files need to be downloaded and unzipped, then update the 'url' argument to the local unzipped file path.

    # make connection to .gaf or .gaf.gz file without downloading it, then read it in.
    con = gzcon(url(url))

    input = readLines(con)
    x = utils::read.delim(textConnection(input),
                   sep = "\t",
                   comment.char = "!",
                   header = FALSE,
                   stringsAsFactors = TRUE)

    # weed out any electronically-inferred evidence-codes (only want curated GO terms)
#    x = x[x[,7]!="IEA",]
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
    today = Sys.Date()
    today = format(today, format = "%b%d%Y")

    # create output file
    GOdb.file = paste("Routput/",
                      organism,
                      "_",
                      plasmoDB_version,
                      "_curated",
                      today,
                      "_GOdb.out",
                      sep = "")

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
    return(GO.db)
    }


## get.GOdef -----
#' @title
#' Get GO terms and definitions
#'
#' @description
#' Extracts the term, term-definition, and ontology for all GOIDs (e.g. GO:0000027) in a geneID2GO object into a dataframe. Output also includes a column for all geneIDs mapping to each term as a comma-separated list. Not required presently for the GO enrichment pipeline, but provides useful context for results--very handy to use as a lookup-table for both GO terms and individual geneIDs.
#'
#'
#' @param geneID2GO  geneID2GO object. Defaults to Pfal_geneID2GO (if using the default, you must first load the data-object \link{Pfal_geneID2GO} ).
#'
#' @export
#'
#' @examples
#' # example code
#' data(Pfal_geneID2GO)  ## load Pfal_geneID2GO object
#' go2genesLookup.df <- get.GOdef(Pfal_geneID2GO)
#'
get.GOdef <- function(geneID2GO = "Pfal_geneID2GO"){
  # get go2genes
  go2genes = topGO::inverseList(geneID2GO)
  go2genes.df = data.frame(sapply(go2genes,stringr::str_flatten_comma))
  go2genes.df$GOID = rownames(go2genes.df)
  rownames(go2genes.df) = NULL
  colnames(go2genes.df) = c("geneIDsList","GOID")
  # grab all GO terms in GOdb
  terms = names(go2genes)
  # extract all relevant term definitions from GO.db package
  GOdef.df = AnnotationDbi::select(x = GO.db, keys = terms, columns = as.character(AnnotationDbi::columns(GO.db)), keytype = "GOID")
  # join GO term definitions with all genes mapping to each term (geneIDs in a comma-separated string for each GOterm)
  GOdef.df = dplyr::left_join(GOdef.df,go2genes.df)
  GOdef.df
}


## get.annot ----
#' @title
#' Extracts and formats annotations from a .gff file.
#'
#' @description
#' Extracts and formats all gene annotations (biotypes "protein_coding_gene", ncRNA_gene", and "pseudogene") from a .gff file. Not required presently for the GO enrichment pipeline, but provides useful context for results. Can be used as-is with a provided gff file as input, or is called by get_pfannot to get the gff file from PlasmoDB.
#'
#' @param x input gff file


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

## get.pfannot ----
#' @title
#' Extracts and formats annotations from a gff file from PlasmoDB.
#'
#' @description
#' Extracts and formats annotations from a gff file from PlasmoDB. Not required presently for the GO enrichment pipeline, but provides useful context for results. Opens a connection to the .gff file from PlasmoDB without downloading it, then calls get_annot() to extract and format the annotation.
#'
#' @param gff_url connection to .gff file. Defaults to \url{https://plasmodb.org/common/downloads/release-66/Pfalciparum3D7/gff/data/PlasmoDB-66_Pfalciparum3D7.gff}


#' @details # \strong{Notes on gff format}
#' The .gff file should be in tabular format with 9 columns, one for each annotated feature associated with a geneID. No formatting is necessary when using the provided url.
#'
#' an annotation created from PlasmoDB's latest \emph{P. falciparum} gff file (accessed December 16, 2021) pre-formatted using this function and ready for run.topGO.meta is included in this package (pf.annot).
#'
#' You \emph{will} need to update the gff url accordingly to the latest version when PlasmoDB is updated.
#'
#' @seealso [get.annot()]
#' @export
get.pfannot <-
  function(gff_url = "https://plasmodb.org/common/downloads/release-66/Pfalciparum3D7/gff/data/PlasmoDB-66_Pfalciparum3D7.gff") {
    # make connection to gff file without downloading it, then read it in.
    con = gzcon(url(gff_url))
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

