% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{Pfal_geneID2GO}
\alias{Pfal_geneID2GO}
\title{All GO terms mapped to \emph{P. falciparum} genes}
\format{
A list of 4477 named vectors--one vector for each \emph{Pf} geneID to which GO terms are mapped. Each vector contains all GO-terms mapped to the geneID.
}
\source{
\url{https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gaf/PlasmoDB-66_Pfalciparum3D7_GO.gaf.gzip}
}
\usage{
data(Pfal_geneID2GO)
}
\description{
A \emph{P. falciparum} GO database containing all GO terms mapped to \emph{P. falciparum} genes (from PlasmoDB release 66; accessed October 18, 2023). Includes automated annotations (evidence code "IEA").

Can be used as the geneID2GO input for run.topGO.meta.
}
\examples{
run.topGO.meta(mydf = exampleMydf, geneID2GO = Pfal_geneID2GO)
}
\keyword{datasets}
