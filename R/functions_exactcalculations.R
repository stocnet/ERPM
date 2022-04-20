######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Functions for exact model calculations                           ##
## Author: Marion Hoffman                                           ##
######################################################################


## --------Calculations for simple models --------

#' Exact estimates number of groups
#'
#'
#' @param num.nodes XXX
#' @param pmin XXX
#' @param pmax XXX
#' @param pmax XXX
#' @return XXX
#' @export
exactestimates_numgroups <- function(num.nodes, pmin, pmax, pinc) {

  allm <- 2:(num.nodes-1)
  allestimates <- rep(0,num.nodes-2)

  for(m in 1:length(allm)) {
    m.obs <- allm[m]
    print(m.obs)

    allparameters <- seq(pmin,pmax,pinc)
    alllogLs <- rep(0,length(allparameters))

    for(i in 1:length(allparameters)) {
      alpha <- allparameters[i]
      numerator <- exp(alpha*m.obs)
      denominator <- compute_numgroups_denominator(num.nodes, alpha)
      #print(paste("numerator",numerator))
      #print(paste("denominator",denominator))
      alllogLs[i] <- log(numerator) - log(denominator)
    }

    allestimates[m] <- allparameters[which(alllogLs == max(alllogLs))[1]]
  }

  return(list(allm,allestimates))
  #plot(allm, allestimates, main = "estimates for different number of groups observed (N=10)")

}


#' Plot likelihood of number groups
#'
#'
#' @param m.obs A partition (vector)
#' @param num.nodes XXX
#' @param pmin Node set (data frame)
#' @param pmax Effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param pinc Objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @return XXX
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


#' Recursive function to compute the value of the denominator
#'
#'
#' @param num.nodes XXX
#' @param alpha Node set (data frame)
#' @return XXX
#' @importFrom utils combn
#' @export
## recursive function to compute the value of the denominator
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
#' Plot the average size distribution
#'
#'
#' @param nmin A partition (vector)
#' @param nmax XXX
#' @return XXX
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


#' Recursive function to compute the value of the denominator
#'
#'
#' @param num.nodes XXX
#' @return XXX
#' @importFrom numbers bell
#' @importFrom utils combn
#' @export
## recursive function to compute the value of the denominator
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

#' Estimation basic model with size constraint
#'
#'
#' @param alpha XXX
#' @param stat XXX
#' @param n XXX
#' @param smin XXX
#' @param smax XXX
#' @return XXX
#' @export
# Estimation basic model with size constraint
calculate_proba_Dirichlet_restricted <- function(alpha,stat,n,smin,smax){
  num <- exp(alpha*stat)
  denom <- calculate_denominator_Dirichlet_restricted(n,smin,smax,alpha,as.list(rep(-1,n+1)))
  return(num/denom$den)
}

#' Estimation basic model with size constraint (same as compute_num_groups_denominator with restrictions)
#'
#'
#' @param n XXX
#' @param smin XXX
#' @param smax XXX
#' @param alpha XXX
#' @param results XXX
#' @return XXX
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
#         maximum=T)
