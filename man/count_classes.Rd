% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_utility.R
\name{count_classes}
\alias{count_classes}
\title{Function to count the number of partitions with a certain
group size structure, for all possible group size structure.
Function to use after calling the "find_all_partitions" function.}
\usage{
count_classes(allpartitions)
}
\arguments{
\item{allpartitions}{matrix containing all possible partitions for a nodeset}
}
\value{
integer(number of partitions with different group structures)
}
\description{
Function to count the number of partitions with a certain
group size structure, for all possible group size structure.
Function to use after calling the "find_all_partitions" function.
}
\examples{
#find partitions first
n <- 6
all_partitions <- find_all_partitions(n)
# count classes
counts_partition_classes <- count_classes(all_partitions)

}
