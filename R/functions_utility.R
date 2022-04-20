######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Utility functions                                                ##
## Author: Marion Hoffman                                           ##
######################################################################



## Functions used for counting partitions and descriptives ############

#' Function to calculate the number of partitions with groups
#' of sizes between smin and smax
#' 
#' @param n XXX
#' @param smin XXX
#' @param smax XXX
#' @return XXX
#' @export
Bell_constraints <- function(n,smin,smax){

  bell <- 0
  for(i in 1:n) {
    bell <- bell + Stirling2_constraints(n,i,smin,smax)
  }
  return(bell)
}


#' Function to calculate the number of partitions with k groups
#' of sizes between smin and smax
#' 
#' @param n XXX
#' @param k XXX
#' @param smin XXX
#' @param smax XXX
#' @return XXX
#' @export
Stirling2_constraints <- function(n,k,smin,smax){

  # base cases
  if(n < k*smin) return(0)
  if(n > k*smax) return(0)
  if(n == k*smin) return( factorial(n) / (factorial(k)*(factorial(smin)^k)) )

  # recurrence relation
  s2 <- 0
  imin <- max(n-smax,0)
  imax <- min(n-smin,n-1)
  for(i in imin:imax){
    fac <- factorial(n-1) / (factorial(i)*factorial(n-1-i))
    s2 <- s2 + fac* Stirling2_constraints(i,k-1,smin,smax)
  }
  return(s2)
}


# Function to calculate the number of partitions with k groups
Stirling <- function(n,k){

  # base cases
  if(n < k || k == 0) return(0)
  if(n == k) return(1)

  # recurrence relation
  s2 <- k * Stirling(n-1,k) + Stirling(n-1,k-1)
  return(s2)
}



## Functions used to create, order, and check partitions in the main code ############

# Find a good starting point for the simple estimation procedure
find_startingpoint_single <- function(nodes,
                                      sizes.allowed){

  num.nodes <- nrow(nodes)

  if(is.null(sizes.allowed)){
    first.partition <- 1 + rbinom(num.nodes, as.integer(num.nodes/2), 0.5)
  } else {
    smin <- min(sizes.allowed)
    cpt <- 0
    g <- 1
    first.partition <- rep(0,num.nodes)
    for(i in 1:num.nodes){
      if(cpt == smin) {
        g <- g + 1
        first.partition[i] <- g
        cpt <- 1
      } else {
        first.partition[i] <- g
        cpt <- cpt + 1
      }
    }
  }
  first.partition <- order_groupids(first.partition)

  return(first.partition)
}

# Find a good starting point for the multiple estimation procedure
find_startingpoint_multiple <- function(presence.tables,
                                        nodes,
                                        sizes.allowed){

  num.nodes <- nrow(nodes)
  num.obs <- ncol(presence.tables)

  first.partitions <- matrix(0, nrow=num.nodes, ncol = num.obs)

  for(o in 1:num.obs){
    nodes.o <- which(presence.tables[,o] == 1)
    num.nodes.o <- sum(presence.tables[,o])
    if(is.null(sizes.allowed)){
      first.partition <- 1 + rbinom(num.nodes.o, as.integer(num.nodes/2), 0.5)
    } else {
      smin <- min(sizes.allowed)
      cpt <- 0
      g <- 1
      first.partition <- rep(0,num.nodes.o)
      for(i in 1:num.nodes.o){
        if(cpt == smin) {
          g <- g + 1
          first.partition[i] <- g
          cpt <- 1
        } else {
          first.partition[i] <- g
          cpt <- cpt + 1
        }
      }
    }
    p.o <- rep(NA,num.nodes)
    p.o[nodes.o] <- order_groupids(first.partition)
    first.partitions[,o] <- p.o
  }

  return(first.partitions)
}


# Function to replace the ids of the group without forgetting an id
# and put in the first appearance order
# for example: [2 1 1 4 2] becomes [1 2 2 3 1]
order_groupids <- function(partition) {

  num.nodes <- length(partition)
  num.groups <- length(unique(partition))
  max.id <- max(partition)

  new.partition <- rep(0,num.nodes)
  cptg <- 1

  for(i in 1:num.nodes){
    if(i == 1) {
      new.partition[i] <- cptg
      cptg <- cptg + 1
      next
    }
    if(sum(partition[1:(i-1)] == partition[i]) == 0) {
      new.partition[i] <- cptg
      cptg <- cptg + 1
    } else {
      g <- new.partition[which(partition[1:(i-1)] == partition[i])[1]]
      new.partition[i] <- g
    }
  }

  return(new.partition)
}


# Function to determine whether a partition contains the allowed group sizes
check_sizes <- function(partition, sizes.allowed){
  allsizes <- unique(table(partition))
  check <- NA %in% match(allsizes,sizes.allowed)
  return(!check)
}




## Other useful functions for descriptives #####################################

min2 <- function(x){
  if(sum(x>0)) {
    return(min(x[x>0]))
  } else {
    return(0)
  }
}

calculate_confidence_interval <- function(vector, conf) {

  average <- mean(vector)
  sd <- sd(vector)
  n <- length(vector)
  error <- qt((conf + 1)/2, df = length(vector) - 1) * sd / sqrt(length(vector))
  result <- c("lower" = average - error, "upper" = average + error)
  return(result)
}
