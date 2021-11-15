#' @title All curated GO terms mapped to \emph{P. falciparum} genes
#' @name Pfal_geneID2GO
#' @docType data
#' @description
#' A \emph{P. falciparum} GO database containing all curated GO terms mapped to \emph{P. falciparum} genes (accessed from PlasmoDB Nov 5, 2021).
#'
#' @usage data(Pfal_geneID2GO)
#' @format A list of 3878 named vectors--one vector for each Pf geneID to which GO terms are mapped. Each vector contains all curated GO-terms mapped to the geneID.
#' @source <https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gaf/PlasmoDB-54_Pfalciparum3D7_GO.gaf>
#' @description can be used as the geneID2GO input for run.topGO.meta.
#' @examples
#' run.topGO.meta(mydf = exampleMydf, geneID2GO = Pfal_geneID2GO)
"Pfal_geneID2GO"




#' @title *piggyBac* mutant classification-data from pooled 1k heat shock screen (to use as example "mydf" parameter for run.topGO.meta()).
#' @name exampleMydf
#' @docType data
#' @description
#' A data frame with \emph{P. falciparum} geneIDs in column 1, and interest-category classifications in column 2. In pooled *piggyBac*-mutant screening, common interest-category classifications are "sensitive", "tolerant", and "neutral". All mutants detected above threshold in the pooled 1k-library heat shock screen are included in this table (n = 752). Derived from https://www.nature.com/articles/s41467-021-24814-1, Table S2.

#'
#' @usage data(exampleMydf)
#' @format A 2-column data frame with 752 \emph{P. falciparum} geneIDs in column 1, and interest-category classifications (here, phenotypes in pooled screening) in column 2.
#' @source <https://www.nature.com/articles/s41467-021-24814-1>
#' @description can be used as example mydf input for run.topGO.meta.
#' @examples
#' run.topGO.meta(mydf = exampleMydf, geneID2GO = Pfal_geneID2GO)
"exampleMydf"




#' @title All \emph{P. falciparum} gene-product annotations.
#' @name pf.annot
#' @docType data
#'
#' @description
#' A dataset containing all \emph{P. falciparum} gene-product annotations taken directly from gff file(reference strain 3D7). Only "gene" entries are kept to remove redundancy.
#'
#' @usage data(pf.annot)
#'
#' @format A data frame of 5562 rows and 9 columns.
#'
#' \describe{
#'   \item{seqid}{chromosome ID.}
#'   \item{source}{annotation source.}
#'   \item{type}{type of annotation. Here only protein-coding genes and ncRNA genes are kept.}
#'   \item{feature_start}{gene start coordinate.}
#'   \item{feature_end}{gene end coordinate.}
#'   \item{strand}{gene strand.}
#'   \item{geneID}{gene ID.}
#'   \item{description}{functional annotation}
#'   \item{geneName}{the gene symbol, or short name, if one exists}

#' }
#'
#' @source <https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gff/data/PlasmoDB-54_Pfalciparum3D7.gff>
#'
#' @keywords dataset
#' @details ## Keep in mind:
#' These data aren't explicitly required for running any enrichments with this package. They are included for reference to aid exploring your enrichment results.
#'
#'
#' Some redundant columns were filtered/removed from the original data source.
#' Updated versions can be generated using the get_pfannot function.
"pf.annot"



#' @title Curated \emph{P. falciparum} gene-sets.
#' @name pf.genesets.mpmp
#' @docType data
#' @description
#' curated \emph{P. falciparum} metabolic pathways, gene-sets, GO terms from MPMP
#' @usage data(pf.genesets.mpmp)
#'
#' @format A large data frame of 54090 rows and 5 columns.
#'
#' \describe{
#'   \item{Gene.ID}{*P. falciparum* gene ID}
#'   \item{Map.Name}{name of the metabolic pathway map/other gene-set curated from publication; searchable at the MPMP website}
#'   \item{Map_id}{the mpmp website address for the pathway map if applicable/the mpmp page for that gene-set. Follow these links to view pathways of interest.)}
#'   \item{Map.cat}{appears to be the gene symbol (short name), or gene-product where available}
#'   \item{pla.id}{appears to be the GO term or enzyme-code where applicable/available}
#' }
#'
#' @source the Malaria Parasite Metabolic Pathways database (MPMP). <https://mpmp.huji.ac.il>
#'
#' @details Not sure about date of access or exact origins of this dataset, but shared from Chengqi Wang in April 2021.
"pf.genesets.mpmp"


#' @title \emph{piggyBac} insertions reported in Zhang et al. 2018 (Science)
#' @name insertions
#' @docType data
#' @description
#' \emph{P. falciparum} insertion-data as published in Science 2018. Contains 3 of the published supplemental tables (Table S1, Pilot_library; Table S3, Saturation_library; Table S5, GenesbyMIS). Not need to run GO enrichment functions, but useful additional info to have for interpreting output.
#'
#' @usage data(insertions)
#' @format A list of 3 dataframes--Pilot_library, corresponding to published table S1; Saturation_library, corresponding to Table S3, and GenesbyMIS, corresponding to Table S5.
#' @source <https://www.science.org/doi/10.1126/science.aap7847>

"insertions"
