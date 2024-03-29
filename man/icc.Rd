% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_statistics.R
\name{icc}
\alias{icc}
\title{Intra class correlation}
\usage{
icc(partition, attribute)
}
\arguments{
\item{partition}{A partition}

\item{attribute}{A vector containing the values of the attribute}
}
\value{
A number corresponding to the ICC
}
\description{
This function computes the intra class correlation correlation
of attributes for 2 randomly drawn individuals in the same group.
}
\examples{
p <- c(1,2,2,3,3,4,4,4,5)
at <- c(3,5,23,2,1,0,3,9,2)
icc(p, at)
}
