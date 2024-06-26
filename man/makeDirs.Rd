% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{makeDirs}
\alias{makeDirs}
\alias{makeRoutput.dir}
\alias{makeGOoutput.dir}
\alias{makeGOhierarchy.dir}
\title{Creating directory-structure}
\usage{
makeRoutput.dir()

makeGOoutput.dir()

makeGOhierarchy.dir()
}
\arguments{
\item{...}{not required. Defaults to creating top-level Routput directory in current working directory.}
}
\description{
Creates the output-directory structure I use for my topGO pipeline. All the included makeDir functions evaluate to the newly created path, or to the existing path if it already exists.
}
\keyword{internal}
