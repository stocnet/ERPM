% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_exactcalculations.R
\name{compute_averagesize}
\alias{compute_averagesize}
\title{Compute the average size of a random partition}
\usage{
compute_averagesize(num.nodes)
}
\arguments{
\item{num.nodes}{number of nodes}
}
\value{
a numeric
}
\description{
Recursive function to compute the average size of a random partition for a given number of nodes
}
\examples{
n <- 6
compute_averagesize(n)

}
