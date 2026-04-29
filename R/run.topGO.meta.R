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
#' You can use the included [formatGOdb.curated()] function to format a custom GO database from curated GeneDB/PlasmoDB annotations for several non-model organisms (or the [formatGOdb()] function to include all GO annotations, if you aren't picky about including automated electronic annotations). If you're studying a model organism, several annotations are already available and can be downloaded through the AnnotationDbi bioconductor package that loads with topGO.
#'
#' @seealso [topGO::topGO()]
#'
#' @examples
#'
#' data(exampleMydf)  ## load exampleMydf object (mutant classification-data from pooled 1k heat shock screen)
#' data(Pfal_geneID2GO_curated)  ## load Pfal_geneID2GO_curated object
#' run.topGO.meta(exampleMydf,Pfal_geneID2GO_curated)
#'
#' @export
run.topGO.meta <- function(mydf = "mydf", geneID2GO = "Pfal_geneID2GO_curated", pval = 0.05) {
#  require(topGO)
#  require(tidyverse)

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

  # print final *significant* GO results table to screen
  cat("\n==============================================================================\n")
  cat("Significant terms from final GO enrichment results table (saved to 'Routput/GO/all.combined.GO.results.tsv':\n")
  print(all.bin.combined.GO.output %>% dplyr::filter(topGO<=pval))

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
