######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Utility functions                                                ##
## Author: Marion Hoffman                                           ##
######################################################################



## Functions used for counting partitions and descriptives ############

#' Function to calculate the number of partitions with groups
#' of sizes between smin and smax
#' 
#' @param n number of nodes
#' @param smin minimum group size possible in the partition
#' @param smax minimum group size possible in the partition
#' @return a numeric
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
#' @param n number of nodes
#' @param k number of groups
#' @param smin minimum group size possible in the partition
#' @param smax maximum group size possible in the partition
#' @return a numeric
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


## Functions to enumerate partitions #########


#' Function to enumerate all partitions for a given n
#' 
#' @param n XXX
#' @return XXX
#' @export
find_all_partitions <- function(n){
  
  if(n == 0) {
    
    # if empty set, nothing
    return(c())
    
  } else {
    
    # find the possibilities for 1 to n-1
    setsn_1 <- find_all_partitions(n-1)
    
    # if n = 1 just return 1
    if(is.null(setsn_1)) 
      return(as.matrix(1))
    
    nsets <- dim(setsn_1)[1]
    res <- c()
    
    # for all solutions, add the singleton n
    for(s in 1:nsets){
      partition <- setsn_1[s,] 
      new <- c(partition,max(partition)+1)
      res <- rbind(res, new)
    }
    
    # for all solutions, add n in all possible groups
    for(s in 1:nsets){
      partition <- setsn_1[s,] 
      ngroups <- max(partition)
      for(g in 1:ngroups) {
        new <- c(partition,g)
        res <- rbind(res, new)
      }
    }
    
    return(res)
  }
  
}


#' Function to count the number of partitions with a certain
#' group size structure, for all possible group size structure
#' 
#' @param allpartitions XXX
#' @return XXX
#' @export
count_classes <- function(allpartitions){
  
  np <- dim(allpartitions)[1]
  n <- dim(allpartitions)[2]
  
  cptclass <- 0
  classes <- matrix()
  classes_count <- c()
  
  for(i in 1:np) {
    partition <- allpartitions[i,]
    
    # find distirbution of blocks
    ngs <- rep(0,n)
    for(g in 1:max(partition)) {
      ngs[sum(partition == g)] <- ngs[sum(partition == g)] + 1 
    }
    
    # find whether there is already classes
    if(cptclass == 0){
      
      cptclass <- 1
      classes <- t(as.matrix(ngs))
      classes_count <- 1
      
    } else {
      
      # check previous classes
      found <- FALSE
      for(c in 1:dim(classes)[1]){
        if(all(classes[c,] == ngs)) {
          found <- TRUE
          classes_count[c] <- classes_count[c] + 1
        }
      }
      
      # if new class
      if(!found){
        cptclass <- cptclass + 1
        classes <- rbind(classes,ngs)
        classes_count <- rbind(classes_count,1)
      }
      
    }
  }
  
  return(list(classes = classes,
              counts = classes_count))
  
}

## Functions used to create, order, and check partitions in the main code ############

#' Function to replace the ids of the group without forgetting an id
#' and put in the first appearance order
#' for example: `[2 1 1 4 2]` becomes `[1 2 2 3 1]`
#' 
#' @param partition observed partition
#' @return a vector (partition)
#' @export
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

#' Function to determine whether a partition contains the allowed group sizes
#' 
#' @param partition observed partition
#' @param sizes.allowed vector containing possible group sizes in the partition
#' @param numgroups.allowed vector containing possible number of groups in the partition
#' @return boolean
#' @export
check_sizes <- function(partition, sizes.allowed, numgroups.allowed){
  allsizes <- unique(table(partition))
  check1 <- (!NA %in% match(allsizes,sizes.allowed))
  numgroup <- length(table(partition))
  check2 <- numgroup %in% numgroups.allowed
  return(check1 & check2)
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




## DEPRECATED ##

## Functions to find starting points for phase 2.
## Deprecated now because we always use the observed partition(s)
## Also, not compatible with the restriction on number of groups

# # Find a good starting point for the simple estimation procedure
# find_startingpoint_single <- function(nodes,
#                                       numgroups.allowed,
#                                       sizes.allowed){
#   
#   num.nodes <- nrow(nodes)
#   
#   if(is.null(sizes.allowed)){
#     first.partition <- 1 + rbinom(num.nodes, as.integer(num.nodes/2), 0.5)
#   } else {
#     smin <- min(sizes.allowed)
#     cpt <- 0
#     g <- 1
#     first.partition <- rep(0,num.nodes)
#     for(i in 1:num.nodes){
#       if(cpt == smin) {
#         g <- g + 1
#         first.partition[i] <- g
#         cpt <- 1
#       } else {
#         first.partition[i] <- g
#         cpt <- cpt + 1
#       }
#     }
#   }
#   first.partition <- order_groupids(first.partition)
#   
#   return(first.partition)
# }
# 
# # Find a good starting point for the multiple estimation procedure
# find_startingpoint_multiple <- function(presence.tables,
#                                         nodes,
#                                         numgroups.allowed,
#                                         sizes.allowed){
#   
#   num.nodes <- nrow(nodes)
#   num.obs <- ncol(presence.tables)
#   
#   first.partitions <- matrix(0, nrow=num.nodes, ncol = num.obs)
#   
#   for(o in 1:num.obs){
#     nodes.o <- which(presence.tables[,o] == 1)
#     num.nodes.o <- sum(presence.tables[,o])
#     if(is.null(sizes.allowed)){
#       first.partition <- 1 + rbinom(num.nodes.o, as.integer(num.nodes/2), 0.5)
#     } else {
#       smin <- min(sizes.allowed)
#       cpt <- 0
#       g <- 1
#       first.partition <- rep(0,num.nodes.o)
#       for(i in 1:num.nodes.o){
#         if(cpt == smin) {
#           g <- g + 1
#           first.partition[i] <- g
#           cpt <- 1
#         } else {
#           first.partition[i] <- g
#           cpt <- cpt + 1
#         }
#       }
#     }
#     p.o <- rep(NA,num.nodes)
#     p.o[nodes.o] <- order_groupids(first.partition)
#     first.partitions[,o] <- p.o
#   }
#   
#   return(first.partitions)
# }
