#' @title All curated GO terms mapped to P. falciparum genes
#' @name Pfal_geneID2GO
#' @docType data
#' @description
#' A P. falciparum GO database containing all curated GO terms mapped to P. falciparum genes (accessed from GeneDB Dec 08, 2020).
#'
#' @usage data(Pfal_geneID2GO)
#' @format A list of 3381 named vectors--one vector for each Pf geneID to which GO terms are mapped. Each vector contains all curated GO-terms mapped to the geneID.
#' @source <ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/Pfalciparum/Pfalciparum.gaf.gz>
#' @description can be used as the geneID2GO input for run.topGO.meta.
#' @example run.topGO.meta(mydf = mydf, geneID2GO = Pfal_geneID2GO)
"Pfal_geneID2GO"


#' @title All P. falciparum gene-product annotations.
#' @name pf.annot
#' @docType data
#'
#' @description
#' A dataset containing all P. falciparum gene-product annotations.
#' @usage data(pf.annot)
#'
#' @format A data frame of 5545 rows and 5 columns.
#'
#' \describe{
#'   \item{organism_name}{P. falciparum)}
#'   \item{transcript_id}{transcript ID. Most genes in Pf are single-transcript (end in .1), but a few are multi-transcript.)}
#'   \item{gene_id}{gene ID. All enrichment-analyses in this package are based on gene ID, not transcript ID.)}
#'   \item{gene_name}{the gene symbol, or short name, if one exists)}
#'   \item{product}{functional annotation)}
#' }
#'
#' @source <ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/Pfalciparum/Pfalciparum.gaf.gz>
#'
#' @keywords dataset
#' @details
#' These data aren't explicitly required for running any enrichments with this package. They are included for reference to aid exploring your enrichment results.
#'
#' Note there are 5545 *transcripts*--not unique geneIDs. There are 5473 unique geneIDs in this dataset.
#'
#' Some redundant columns were filtered/removed from the original data source.
"pf.annot"

#' @title Curated P. falciparum gene-sets.
#' @name pf.genesets.mpmp
#' @docType data
#' @description
#' curated P. falciparum metabolic pathways, gene-sets, GO terms from MPMP
#' @usage data(pf.genesets.mpmp)
#'
#' @format A large data frame of 54090 rows and 5 columns.
#'
#' \describe{
#'   \item{Gene.ID}{P. falciparum gene ID)}
#'   \item{Map.Name}{name of the metabolic pathway map/other gene-set curated from publication; searchable at the MPMP website)}
#'   \item{Map_id}{the mpmp website address for the pathway map if applicable/the mpmp page for that gene-set. Follow these links to view pathways of interest.)}
#'   \item{Map.cat}{appears to be the gene symbol (short name), or gene-product where available)}
#'   \item{pla.id}{appears to be the GO term or enzyme-code where applicable/available)}
#' }
#'
#' @source the Malaria Parasite Metabolic Pathways database (MPMP). <https://mpmp.huji.ac.il>
#'
#' @details Not sure about date of access or exact origins of this dataset, but shared from Chengqi Wang in April 2020.
"pf.genesets.mpmp"
