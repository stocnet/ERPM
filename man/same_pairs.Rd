% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_statistics.R
\name{same_pairs}
\alias{same_pairs}
\title{Same pairs of individuals in a partition}
\usage{
same_pairs(partition, attribute, stat)
}
\arguments{
\item{partition}{A partition (vector)}

\item{attribute}{A vector containing the values of the attribute}

\item{stat}{The statistic to compute : 'avg_pergroup' for the average, 'sum_pergroup' for the sum,  'sum_perind' and 'avg_perind'  for the number of ties per individual
each individual has in its group.}
}
\value{
The statistic chosen in stat
}
\description{
This function computes the total number, the average number having the same value of a categorical variable
and the number of individuals a partition.
}
\examples{
p <- c(1,2,2,3,3,4,4,4,5)
at <- c(0,1,1,1,1,0,0,0,0)
same_pairs(p,at,'avg_pergroup')
}
