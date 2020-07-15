######################################################################
## Simulation and estimation of Exponential Random Partition Models ## 
## Utility functions                                                ##
## Author: Marion Hoffman                                           ##
######################################################################



## Functions used for counting partitions and descriptives ############

# Function to calculate the number of partitions with groups
# of sizes between smin and smax
Bell_constraints <- function(n,smin,smax){

  bell <- 0
  for(i in 1:n) {
    bell <- bell + Stirling2_constraints(n,i,smin,smax)
  }
  return(bell)
}


# Function to calculate the number of partitions with k groups
# of sizes between smin and smax
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


# Function to calculate the average number of groups
# in a random partition for a given n
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
        sum <- sum + choose(n,k) * bell(n-k) * (n+1) / (1/bell(n-k) + (n-k)/compute_averagesize(n-k) )
      }
    }
    
    return(sum/bell(num.nodes))
  }
  
}


# Function to enumerate all partitions for a given n
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


# Function to count the number of partitions with a certain
# group size structure, for all possible group size structure
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
  # while(max.id != num.groups) {
  #   
  #   # find the first missing id
  #   all.ids <- 1:max.id
  #   first <- all.ids[!(all.ids %in% partition)][1]
  #   
  #   # remove 1 from the ids following the missing id
  #   ids.more <- which(partition > first)
  #   partition[ids.more] <- partition[ids.more] - 1
  #   
  #   # go to the next
  #   max.id <- max.id - 1
  #   
  # }
  # 
  # return(partition)
}


# Function to determine whether a partition contains the allowed group sizes
check_sizes <- function(partition, sizes.allowed){
  allsizes <- unique(table(partition))
  check <- NA %in% match(allsizes,sizes.allowed)
  return(!check)
}






## Other useful functions for descriptives #####################################

library(clue)
calculate_min_transferdistances <- function(allpartitions, allclasses){
  
  n <- dim(allpartitions)[2]
  np <- dim(allpartitions)[1]
  nc <- dim(allclasses)[1]
  distmatrix <- matrix(1,nc,nc)
  
  for(p1 in 1:(np-1)){

    #define the first partition
    part1 <- allpartitions[p1,]
    ngs <- rep(0,n)
    for(g in 1:max(part1)) {
      ngs[sum(part1 == g)] <- ngs[sum(part1 == g)] + 1
    }
    for(c in 1:dim(allclasses)[1]){
      if(all(allclasses[c,] == ngs)) {
        class1 <- c
      }
    }

    for(p2 in (p1+1):np) {

      # define the second
      part2 <- allpartitions[p2,]
      ngs <- rep(0,n)
      for(g in 1:max(part2)) {
        ngs[sum(part2 == g)] <- ngs[sum(part2 == g)] + 1
      }
      for(c in 1:dim(allclasses)[1]){
        if(all(allclasses[c,] == ngs)) {
          class2 <- c
        }
      }

      # calculate the similarity matrix
      m <- max(max(part1),max(part2))
      sim <- matrix(0,m,m)
      for(node in 1:n){
        sim[part1[node],part2[node]] <- sim[part1[node],part2[node]] + 1
      }

      #calculate the maximal assignment with hungarian method and similarity
      maxass <- solve_LSAP(sim, maximum=T)
      s <- 0
      for(g in 1:m){
        s <- s + sim[g,maxass[g]]
      }

      # calculate the transfer distance
      t <- n - s
      distmatrix[c1,c2] <- t
      distmatrix[c2,c1] <- t
          
      if(t < distmatrix[class1,class2]){
        distmatrix[class1,class2] <- t
        distmatrix[class2,class1] <- t
      }

    }

  }

  diag(distmatrix) <- 1
  
  return(distmatrix)
}

calculate_min_Randdistances <- function(allpartitions, allclasses){
  
  n <- dim(allpartitions)[2]
  np <- dim(allpartitions)[1]
  nc <- dim(allclasses)[1]
  distmatrix <- matrix(1,nc,nc)
  
  for(c1 in 1:(nc-1)){
    ngs1 <- allclasses[c1,]
    
    for(c2 in (c1+1):nc){
      ngs2 <- allclasses[c2,]
      
      # find first ordered partition
      cpt1 <- 1
      cpt2 <- n
      g <- 1
      part1 <- rep(0,n)
      tempng <- ngs1
      
      while(cpt1 <= n){
        foundcpt2 <- (tempng[cpt2] > 0)
        while(!foundcpt2){
          cpt2 <- cpt2 - 1
          foundcpt2 <- (tempng[cpt2] > 0)
        }
        
        part1[cpt1:(cpt1+cpt2-1)] <- g
        g <- g+1
        cpt1 <- cpt1+cpt2
        tempng[cpt2] <- tempng[cpt2]-1
      }
      
      # find second ordered partition
      cpt1 <- 1
      cpt2 <- n
      g <- 1
      part2 <- rep(0,n)
      tempng <- ngs2
      
      while(cpt1 <= n){
        foundcpt2 <- (tempng[cpt2] > 0)
        while(!foundcpt2){
          cpt2 <- cpt2 - 1
          foundcpt2 <- (tempng[cpt2] > 0)
        }
        
        part2[cpt1:(cpt1+cpt2-1)] <- g
        g <- g+1
        cpt1 <- cpt1+cpt2
        tempng[cpt2] <- tempng[cpt2]-1
      }
      
      # calculate r
      r <- adjustedRand(part1, part2, randMethod = "Rand")
      distmatrix[c1,c2] <- r
      distmatrix[c2,c1] <- r
    }
    
  }
  
  
  # for(p1 in 1:(np-1)){
  #   
  #   #define the first partition
  #   part1 <- allpartitions[p1,]
  #   ngs <- rep(0,n)
  #   for(g in 1:max(part1)) {
  #     ngs[sum(part1 == g)] <- ngs[sum(part1 == g)] + 1 
  #   }
  #   for(c in 1:dim(allclasses)[1]){
  #     if(all(allclasses[c,] == ngs)) {
  #       class1 <- c
  #     }
  #   }
  #   
  #   for(p2 in (p1+1):np) {
  #     
  #     # define the second
  #     part2 <- allpartitions[p2,]
  #     ngs <- rep(0,n)
  #     for(g in 1:max(part2)) {
  #       ngs[sum(part2 == g)] <- ngs[sum(part2 == g)] + 1 
  #     }
  #     for(c in 1:dim(allclasses)[1]){
  #       if(all(allclasses[c,] == ngs)) {
  #         class2 <- c
  #       }
  #     }
  #     
  #     r <- adjustedRand(part1,part2,randMethod="Rand")
  #     
  #     if(r > distmatrix[class1,class2]){
  #       distmatrix[class1,class2] <- r
  #       distmatrix[class2,class1] <- r
  #     }
  #     
  #   }
  #   
  # }
  # 
  # diag(distmatrix) <- 1
  
  return(distmatrix)
}

calculate_average_Hammingdistances <- function(allclasses) {
  
  nc <- dim(allclasses)[1]
  n <- dim(allclasses)[2]
  distmatrix <- matrix(1,nc,nc)
  
  for(c1 in 1:(nc-1)){
    ngs1 <- allclasses[c1,]
    
    # find number of elements off diagonal on average
    noffdiag <- 0
    for(g in 1:n){
      noffdiag <- noffdiag + ngs1[g]*(g*(g-1))
    }
    av_partition1 <- matrix(noffdiag/(n*(n-1)),n,n)
    diag(av_partition1) <- 1
    
    for(c2 in (c1+1):nc){
      ngs2 <- allclasses[c2,]
      
      # find number of elements off diagonal on average
      noffdiag <- 0
      for(g in 1:n){
        noffdiag <- noffdiag + ngs2[g]*(g*(g-1))
      }
      av_partition2 <- matrix(noffdiag/(n*(n-1)),n,n)
      diag(av_partition2) <- 1
      
      # calculate r
      r <- adjustedRand(av_partition1, av_partition2, randMethod = "Rand")
      distmatrix[c1,c2] <- r
      distmatrix[c2,c1] <- r
    }
    
  }
  
}

# ICC for binary/continuous variables
computeicc <- function(attribute,teams) {
  
  n <- length(attribute)
  m <- max(teams,na.rm=T)
  totalmean <- mean(attribute,na.rm=T)
  
  s2b <- 0
  s2w <- 0
  
  for(h in 1:m){
    #if(h==6) h <- 15
    members <- which(teams==h)
    groupmean <- mean(attribute[members])
    tempsumb <- 0
    tempsumw <- 0
    for(i in 1:length(members)){
      tempsumb <- tempsumb + attribute[i] - totalmean
      tempsumw <- tempsumw + (attribute[i] - groupmean)^2
    }
    tempsumb <- tempsumb/length(members)
    s2b <- s2b + tempsumb^2 
    s2w <- s2w + tempsumw
  }
  s2b <- s2b/(m-1)
  s2w <- s2w/(n-m)
  
  icc <- s2b / (s2b + s2w)
  return(icc)
}

# Density for network/similarity dyadic variables (binary, symetric)
computedensity <- function(matrix,teams){
  
  m <- max(teams,na.rm=T)
  avd <- 0
  std <- 0
  
  for(h in 1:m){
    members <- which(teams==h)
    
    if(length(members)>1){
      beforeties <- 0
      allties <- length(members)*(length(members)-1)/2
    
      for(i in 1:(length(members)-1)){
        for(j in (i+1):length(members)){
          beforeties <- beforeties + matrix[members[i],members[j]]
        }
      }
    
      avd <- avd + beforeties/allties
      std <- std + (beforeties/allties)^2
    }
  }
  
  avd <- avd/m
  std <- sqrt(1/m*std - avd^2)
  return(list(average=avd,
              standard.deviation=std))
}


# Correlation between variable and size
computecorrelationwithsize <- function(attribute,teams){
  
  n <- length(attribute)
  m <- max(teams)
  
  sizes <- 0
  
  for(i in 1:n){
    t <- teams[i]
    sizes[i] <- length(which(teams==t))
  }
  
  c <- cor(attribute,sizes)
  return(c)
}


# convert networks into partitions
fromnetworktopartition <- function(network){
  
  n <- nrow(network)
  partition <- rep(0,n)
  cptg <- 1
  
  for(i in 1:n){
    if(i == 1){
      partition[1] <- cptg
      cptg <- cptg+1
      next
    }
    if(sum(network[i,1:(i-1)])==0) {
      partition[i] <- cptg
      cptg <- cptg + 1 
    } else {
      g <- partition[which(network[i,1:(i-1)] > 0)[1]]
      partition[i] <- g
    }
  }
  
  return(partition)
}
