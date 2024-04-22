######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Functions for exact model calculations                           ##
## Author: Marion Hoffman                                           ##
######################################################################


## --------Calculations for simple models --------

#' Exact estimates number of groups
#'
#' This function finds the best estimate for a model only including the statistics of number of groups.
#' It does a grid search for a vector of potential parameters, for all numbers of groups.
#'
#' @param num.nodes number of nodes
#' @param pmin lowest parameter value
#' @param pmax highest parameter value
#' @param pinc increment between different parameter values
#' @param verbose logical: should intermediate results be printed or not? Defaults to FALSE.
#' @return a list
#' @export
exactestimates_numgroups <- function(num.nodes, pmin, pmax, pinc, verbose = FALSE) {

  allm <- 2:(num.nodes-1)
  allestimates <- rep(0,num.nodes-2)

  for(m in 1:length(allm)) {
    m.obs <- allm[m]
    if (verbose) cat(m.obs, "\n")

    allparameters <- seq(pmin,pmax,pinc)
    alllogLs <- rep(0,length(allparameters))

    for(i in 1:length(allparameters)) {
      alpha <- allparameters[i]
      numerator <- exp(alpha*m.obs)
      denominator <- compute_numgroups_denominator(num.nodes, alpha)
      # if (verbose) {
      #   cat(paste("numerator",numerator), "\n")
      #   cat(paste("denominator",denominator), "\n\n")
      # }
      alllogLs[i] <- log(numerator) - log(denominator)
    }

    allestimates[m] <- allparameters[which(alllogLs == max(alllogLs))[1]]
  }

  return(list(allm,allestimates))

}


#' Plot likelihood of number groups
#'
#' Function to plot the log-likelihood of the model with a single statistic (number of groups)
#' depending on the parameter value for this statistic
#'
#' @param m.obs observed number of groups
#' @param num.nodes number of nodes
#' @param pmin lowest parameter value
#' @param pmax highest parameter value
#' @param pinc increment between different parameter values
#' @return a vector
#' @export
plot_numgroups_likelihood <- function(m.obs, num.nodes, pmin, pmax, pinc) {

  allparameters <- seq(pmin,pmax,pinc)
  alllogLs <- rep(0,length(allparameters))

  for(i in 1:length(allparameters)) {
    alpha <- allparameters[i]
    numerator <- exp(alpha*m.obs)
    denominator <- compute_numgroups_denominator(num.nodes, alpha)
    #print(paste("numerator",numerator))
    #print(paste("denominator",denominator))
    alllogLs[i] <- numerator / denominator
  }

  plot(allparameters, alllogLs, main="log likelihood")

}


#' Compute denominator for model with number of groups
#' 
#' Recursive function to compute the value of the denominator for the model
#' with a single statistic which is the number of groups
#'
#'
#' @param num.nodes number of nodes
#' @param alpha parameter value
#' @return a numeric
#' @importFrom utils combn
#' @export
compute_numgroups_denominator <- function( num.nodes, alpha){

  # if no nodes, by convention we return 1
  if(num.nodes == 0) {
    return(1)

  } else {

    n <- num.nodes - 1
    sum <- 0
    for(k in 0:n){
      if( k==0 ) {
        sum <- sum + compute_numgroups_denominator(n, alpha)
      } else if (k == n){
        sum <- sum + compute_numgroups_denominator(0, alpha)
      } else {
        sum <- sum + dim(combn(n,k))[2]*compute_numgroups_denominator(n-k, alpha)
      }
    }

    return(exp(alpha)*sum)
  }

}
#' Plot average sizes
#'
#' Function to plot the average size of a random partition depending on the number of nodes
#'
#' @param nmin minimum number of nodes
#' @param nmax maximum number of nodes
#' @param ninc increment between the different number of nodes
#' @return a vector
#' @export
plot_averagesizes <- function(nmin, nmax, ninc) {

  allns <- seq(nmin,nmax,ninc)
  allsizes <- rep(0,length(allns))

  for(i in 1:length(allns)) {
    n <- allns[i]
    size <- compute_averagesize(n)
    print(paste("n",n))
    print(paste("average size",size))
    allsizes[i] <- size
  }

  plot(allns, allsizes, main="average size for a random partition depending on the number of nodes")

}


#' Compute the average size of a random partition 
#' 
#' Recursive function to compute the average size of a random partition for a given number of nodes
#'
#'
#' @param num.nodes number of nodes
#' @return a numeric
#' @importFrom numbers bell
#' @importFrom utils combn
#' @export 
#' @examples
#' n <- 6
#' compute_averagesize(n)
#' 
compute_averagesize <- function( num.nodes){

  # if no nodes, by convention we return 1
  if(num.nodes == 1) {
    return(1)

  } else {

    n <- num.nodes - 1
    sum <- 0

    for(k in 0:n){
      if( k==0 ) {
        sum <- sum + bell(n) *  (n+1) / (1/bell(n) + n/compute_averagesize(n) )
      } else if (k==n){
        sum <- sum + bell(0) * (n+1)
      } else {
        sum <- sum + dim(combn(n,k))[2] * bell(n-k) * (n+1) / (1/bell(n-k) + (n-k)/compute_averagesize(n-k) )
      }
    }

    return(sum/bell(num.nodes))
  }

}

## -------- Exact estimation simple models --------

#' Calculate Dirichlet probability
#'
#' Calculate the probability of observing a partition with a given number of groups
#' for a model with a single statistic for the number of groups and a given parameter value.
#' The set of possible partitions can be restricted to partitions with groups of a certain size.
#' 
#'
#' @param alpha parameter value
#' @param stat observed stat (number of groups)
#' @param n number of nodes
#' @param smin minimum size for a group
#' @param smax maximum size for a group
#' @return a numeric
#' @export
calculate_proba_Dirichlet_restricted <- function(alpha,stat,n,smin,smax){
  num <- exp(alpha*stat)
  denom <- calculate_denominator_Dirichlet_restricted(n,smin,smax,alpha,as.list(rep(-1,n+1)))
  return(num/denom$den)
}

#' Calculate Dirichlet denominator
#' 
#' Recursive function to calculate the denominator for the model with a single statistic 
#' for the number of groups and a given parameter value.
#' The set of possible partitions can be restricted to partitions with groups of a certain size.
#'
#' @param n number of nodes
#' @param smin minimum size for a group
#' @param smax maximum size for a group
#' @param alpha parameter value
#' @param results a list
#' @return a numeric
#' @export
calculate_denominator_Dirichlet_restricted <- function(n,smin,smax,alpha, results){
  if(n == 0){
    results[[1]] <- 1
    return(list(den=1,res=results))
  }
  if(n < smin){
    results[[n+1]] <- 0
    return(list(den=0,res=results))
  }
  if(n == smin){
    results[[smin+1]] <- exp(alpha)
    return(list(den=exp(alpha),res=results))
  }
  if(n > smin){
    if(results[[n+1]] > -1){
      return(list(den=results[[n+1]],res=results))
    } else {
      sum <- 0
      imin <- max(0,n-smax)
      imax <- min(n-1,n-smin)
      if(imax >= imin){
        for(i in imin:imax){
          if(results[[i+1]] > -1) {
            deni <- results[[i+1]]
          } else {
            calci <- calculate_denominator_Dirichlet_restricted(i,smin,smax,alpha,results)
            deni <- calci$den
            results <- calci$res
          }
          sum <- sum + choose(n-1,i)*exp(alpha)*deni
        }
      }
      results[[n+1]] <- sum
      return(list(den=sum,res=results))
    }

  }
}

#library(DEoptimR)
#optimize(f = function(x){ return(calculate_proba_restricted(x,14,58,3,5))},
#         interval = c(-17,10),
#         maximum = TRUE)

