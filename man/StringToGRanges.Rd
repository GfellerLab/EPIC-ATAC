% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EPIC_fun.R
\name{StringToGRanges}
\alias{StringToGRanges}
\title{Converting open chromatin regions (i.e. peaks) to GRanges objects.}
\usage{
StringToGRanges(regions, sep = c("-", "-"), ...)
}
\arguments{
\item{regions}{: vector of coordinates in the following format: "chr-start-end".}

\item{sep}{(optional, default is c("-", "-"): a vector of length 2 containing the separators used in the coordinates}
}
\value{
a GRanges object.
}
\description{
Converting open chromatin regions (i.e. peaks) to GRanges objects.
}
\keyword{internal}
