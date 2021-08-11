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


# Function to calculate the number of partitions with k groups
Stirling <- function(n,k){
  
  # base cases
  if(n < k || k == 0) return(0)
  if(n == k) return(1)
  
  # recurrence relation
  s2 <- k * Stirling(n-1,k) + Stirling(n-1,k-1)
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

#library(clue)
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


# Proportion of majority attribute in groups (for binary or two categories)
computeproportion <- function(attribute,teams){
  
  attribute <- as.numeric(factor(attribute))
  maxatt <- max(attribute)
  minatt <- min(attribute)
  
  sum <- 0
  for(g in 1:max(teams, na.rm = T)){
    members <- which(teams == g)
    nmax <- sum(attribute[members] == maxatt)
    nmin <- sum(attribute[members] == minatt)
    sum <- sum + (max(nmax,nmin) - min(nmax,nmin)) / length(members)
  }
  
  return(list(sum = sum,
              average = sum/max(teams, na.rm = T)))
}

# Range of attribute in groups (for non binary)
computerange <- function(attribute,teams){

  sum <- 0
  for(g in 1:max(teams, na.rm = T)){
    members <- which(teams == g)
    sum <- sum + max(attribute[members], na.rm = T) - min(attribute[members], na.rm = T)
  }
  
  return(list(sum = sum,
              average = sum/max(teams, na.rm = T)))
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

## CUG functions ############################♠#####################################

# CUG test for network: binary attribute (also works for categorical)
CUG_network_binary <- function(network_sample, network_obs, attribute) {
  
  dist <- as.matrix(dist(as.numeric(factor(attribute)), diag = T, upper = T))
  distances <- dist
  distances[dist == 0] <- 1
  distances[dist != 0] <- 0
  diag(distances) <- 0
  
  sameties_obs <- sum(distances * network_obs)
  avsameties_obs <- sum(distances * network_obs) / sum(network_obs)
  
  allsameties <- rep(0,length(network_sample))
  allavsameties <- rep(0,length(network_sample))
  for(s in 1:length(network_sample)){
    allsameties[s] <- sum(distances * network_sample[[s]])
    allavsameties[s] <- sum(distances * network_sample[[s]]) / sum(network_sample[[s]])
  }
  
  minsame <- quantile(allsameties, probs = 0.05)
  maxsame <- quantile(allsameties, probs = 0.95)
  minavsame <- quantile(allavsameties, probs = 0.05)
  maxavsame <- quantile(allavsameties, probs = 0.95)
  
  p1 <- sum(allsameties <= sameties_obs) / length(network_sample)
  p2 <- sum(allsameties >= sameties_obs) / length(network_sample)
  if(p1 == 0) {
    p_same <- p2
  } else if(p2 == 0) {
    p_same <- p1
  } else {
    p_same <- 2*min(p1,p2)
  }
  
  p1 <- sum(allavsameties <= avsameties_obs) / length(network_sample)
  p2 <- sum(allavsameties >= avsameties_obs) / length(network_sample)
  if(p1 == 0) {
    p_avsame <- p2
  } else if(p2 == 0) {
    p_avsame <- p1
  } else {
    p_avsame <- 2*min(p1,p2)
  }
  
  res <- list(number_same_ties = list(observed = sameties_obs,
                                      mean_sample = mean(allsameties),
                                      sd_sample = sd(allsameties),
                                      min_CI_95 = minsame,
                                      max_CI_95 = maxsame,
                                      p = p_same),
              average_number_same_ties = list(observed = avsameties_obs,
                                              mean_sample = mean(allavsameties),
                                              sd_sample = sd(allavsameties),
                                              min_CI_95 = minavsame,
                                              max_CI_95 = maxavsame,
                                              p = p_avsame))
    
  return(res)
}

# CUG test for network: dyadic attribute
CUG_network_dyadic <- function(network_sample, network_obs, net_attribute) {
  
  numties_obs <- sum(net_attribute * network_obs)
  avnumties_obs <- sum(net_attribute * network_obs) / sum(network_obs)
  
  allnumties <- rep(0,length(network_sample))
  allavnumties <- rep(0,length(network_sample))
  for(s in 1:length(network_sample)){
    allnumties[s] <- sum(net_attribute * network_sample[[s]])
    allavnumties[s] <- sum(net_attribute * network_sample[[s]]) / sum(network_sample[[s]])
  }
  
  minnum <- quantile(allnumties, probs = 0.05)
  maxnum <- quantile(allnumties, probs = 0.95)
  minavnum <- quantile(allavnumties, probs = 0.05)
  maxavnum <- quantile(allavnumties, probs = 0.95)
  
  p1 <- sum(allnumties <= numties_obs) / length(network_sample)
  p2 <- sum(allnumties >= numties_obs) / length(network_sample)
  if(p1 == 0) {
    p_num <- p2
  } else if(p2 == 0) {
    p_num <- p1
  } else {
    p_num <- 2*min(p1,p2)
  }
  
  p1 <- sum(allavnumties <= avnumties_obs) / length(network_sample)
  p2 <- sum(allavnumties >= avnumties_obs) / length(network_sample)
  if(p1 == 0) {
    p_avnum <- p2
  } else if(p2 == 0) {
    p_avnum <- p1
  } else {
    p_avnum <- 2*min(p1,p2)
  }
  
  res <- list(number_ties = list(observed = numties_obs,
                                 mean_sample = mean(allnumties),
                                 sd_sample = sd(allnumties),
                                 min_CI_95 = minnum,
                                 max_CI_95 = maxnum,
                                 p = p_num),
              average_number_ties = list(observed = avnumties_obs,
                                         mean_sample = mean(allavnumties),
                                         sd_sample = sd(allavnumties),
                                         min_CI_95 = minavnum,
                                         max_CI_95 = maxavnum,
                                         p = p_avnum))
  
  return(res)
}

# CUG test for network: continuous attribute (or almost continuous)
CUG_network_continous <- function(network_sample, network_obs, attribute) {
  
  distances <- as.matrix(dist(attribute, diag = T, upper = T))
  
  diffs_obs <- sum(distances * network_obs)
  avdiffs_obs <- sum(distances * network_obs) / sum(network_obs)
  
  alldiffs <- rep(0,length(network_sample))
  allavdiffs <- rep(0,length(network_sample))
  for(s in 1:length(network_sample)){
    alldiffs[s] <- sum(distances * network_sample[[s]])
    alldiffs[s] <- sum(distances * network_sample[[s]]) / sum(network_sample[[s]])
  }
  
  mindiffs <- quantile(alldiffs, probs = 0.05)
  maxdiffs <- quantile(alldiffs, probs = 0.95)
  minavdiffs <- quantile(allavdiffs, probs = 0.05)
  maxavdiffs <- quantile(allavdiffs, probs = 0.95)
  
  p1 <- sum(alldiffs <= diffs_obs) / length(network_sample)
  p2 <- sum(alldiffs >= diffs_obs) / length(network_sample)
  if(p1 == 0) {
    p_diffs <- p2
  } else if(p2 == 0) {
    p_diffs <- p1
  } else {
    p_diffs <- 2*min(p1,p2)
  }
  
  p1 <- sum(allavdiffs <= avdiffs_obs) / length(network_sample)
  p2 <- sum(allavdiffs >= avdiffs_obs) / length(network_sample)
  if(p1 == 0) {
    p_avdiffs <- p2
  } else if(p2 == 0) {
    p_avdiffs <- p1
  } else {
    p_avdiffs <- 2*min(p1,p2)
  }
  
  res <- list(sum_absolute_differences = list(observed = diffs_obs,
                                              mean_sample = mean(alldiffs),
                                              sd_sample = sd(alldiffs),
                                              min_CI_95 = mindiffs,
                                              max_CI_95 = maxdiffs,
                                              p = p_diffs),
              average_absolute_difference = list(observed = avdiffs_obs,
                                                 mean_sample = mean(allavdiffs),
                                                 sd_sample = sd(allavdiffs),
                                                 min_CI_95 = minavdiffs,
                                                 max_CI_95 = maxavdiffs,
                                                 p = p_avdiffs))
  
  return(res)
}

# CUG test for partition: ONE attribute and a given tye for the partition groups
CUG_partition_uniquebinary_withtypes <- function(partition_sample, partition_obs, types_obs, attribute, type) {
  
  partition_obs2 <- partition_obs
  out_groups <- which(types_obs != type)
  partition_obs2[partition_obs2 %in% out_groups] <- NA
  partition_obs2[!is.na(partition_obs2)] <- order_groupids(partition_obs2[!is.na(partition_obs2)])
  
  partition_sample2 <- list()
  for(s in 1:length(partition_sample)) {
    p2 <- partition_sample[[s]]$partition
    out_groups <- which(partition_sample[[s]]$types != type)
    p2[p2 %in% out_groups] <- NA
    p2[!is.na(p2)] <- order_groupids(p2[!is.na(p2)])
    partition_sample2[[s]] <- p2
  }
  
  numinds_obs <- sum(attribute * !is.na(partition_obs2)) 
  avgroupnuminds_obs <- sum(attribute * !is.na(partition_obs2))  / (max(partition_obs2, na.rm = T))

  allnuminds <- rep(0,length(partition_sample))
  allavgroupnuminds <- rep(0,length(partition_sample))
  for(s in 1:length(partition_sample)){
    allnuminds[s] <- sum(attribute * !is.na(partition_sample2[[s]])) 
    allavgroupnuminds[s] <- sum(attribute * !is.na(partition_sample2[[s]])) / (max(partition_sample2[[s]], na.rm=T))
  }
  
  minnum <- quantile(allnuminds, probs = 0.05)
  maxnum <- quantile(allnuminds, probs = 0.95)
  minavgroupnum <- quantile(allavgroupnuminds, probs = 0.05)
  maxavgroupnum <- quantile(allavgroupnuminds, probs = 0.95)
  
  p1 <- sum(allnuminds <= numinds_obs) / length(partition_sample)
  p2 <- sum(allnuminds >= numinds_obs) / length(partition_sample)
  if(p1 == 0) {
    p_num <- p2
  } else if(p2 == 0) {
    p_num <- p1
  } else {
    p_num <- 2*min(p1,p2)
  }
  
  p1 <- sum(allavgroupnuminds <= avgroupnuminds_obs) / length(partition_sample)
  p2 <- sum(allavgroupnuminds >= avgroupnuminds_obs) / length(partition_sample)
  if(p1 == 0) {
    p_avgroupnum <- p2
  } else if(p2 == 0) {
    p_avgroupnum <- p1
  } else {
    p_avgroupnum <- 2*min(p1,p2)
  }
  
  res <- list(number_individuals = list(observed = numinds_obs,
                                        mean_sample = mean(allnuminds),
                                        sd_sample = sd(allnuminds),
                                        min_CI_95 = minnum,
                                        max_CI_95 = maxnum,
                                        p = p_num),
              averagepergroup_number_individuals = list(observed = avgroupnuminds_obs,
                                                        mean_sample = mean(allavgroupnuminds),
                                                        sd_sample = sd(allavgroupnuminds),
                                                        min_CI_95 = minavgroupnum,
                                                        max_CI_95 = maxavgroupnum,
                                                        p = p_avgroupnum))
  
  return(res)
}

# CUG test for partition: ONE continuous attribute and a given tye for the partition groups
CUG_partition_uniquecontinuous_withtypes <- function(partition_sample, partition_obs, types_obs, attribute, type) {
  
  partition_obs2 <- partition_obs
  out_groups <- which(types_obs != type)
  partition_obs2[partition_obs2 %in% out_groups] <- NA
  partition_obs2[!is.na(partition_obs2)] <- order_groupids(partition_obs2[!is.na(partition_obs2)])
  
  partition_sample2 <- list()
  for(s in 1:length(partition_sample)) {
    p2 <- partition_sample[[s]]$partition
    out_groups <- which(partition_sample[[s]]$types != type)
    p2[p2 %in% out_groups] <- NA
    p2[!is.na(p2)] <- order_groupids(p2[!is.na(p2)])
    partition_sample2[[s]] <- p2
  }
  
  avatt_obs <- mean(attribute * !is.na(partition_obs2)) 
  avgroupavatt_obs <- 0
  for(g in 1:max(partition_obs2, na.rm = T)){
    avgroupavatt_obs <- avgroupavatt_obs + mean(attribute[which(partition_obs2 == g)])
  }
  avgroupavatt_obs <- avgroupavatt_obs / (max(partition_obs2, na.rm = T))
  
  allavatt <- rep(0,length(partition_sample))
  allavgroupavatt <- rep(0,length(partition_sample))
  for(s in 1:length(partition_sample)){
    allavatt[s] <- mean(attribute * !is.na(partition_sample2[[s]])) 
    allavgroupavatt[s] <- 0
    for(g in 1:max(partition_sample2[[s]], na.rm=T)){
      allavgroupavatt[s] <- allavgroupavatt[s] + mean(attribute[which(partition_sample2[[s]] == g)])
    }
    allavgroupavatt[s] <- allavgroupavatt[s] / max(partition_sample2[[s]], na.rm = T)
  }
  
  minavatt <- quantile(allavatt, probs = 0.05)
  maxavatt <- quantile(allavatt, probs = 0.95)
  minavgroupavatt <- quantile(allavgroupavatt, probs = 0.05)
  maxavgroupavatt <- quantile(allavgroupavatt, probs = 0.95)
  
  p1 <- sum(allavatt <= avatt_obs) / length(partition_sample)
  p2 <- sum(allavatt >= avatt_obs) / length(partition_sample)
  if(p1 == 0) {
    p_avatt <- p2
  } else if(p2 == 0) {
    p_avatt <- p1
  } else {
    p_avatt <- 2*min(p1,p2)
  }
  
  p1 <- sum(allavgroupavatt <= avgroupavatt_obs) / length(partition_sample)
  p2 <- sum(allavgroupavatt >= avgroupavatt_obs) / length(partition_sample)
  if(p1 == 0) {
    p_avgroupavatt <- p2
  } else if(p2 == 0) {
    p_avgroupavatt <- p1
  } else {
    p_avgroupavatt <- 2*min(p1,p2)
  }
  
  res <- list(mean_attribute = list(observed = avatt_obs,
                                    mean_sample = mean(allavatt),
                                    sd_sample = sd(allavatt),
                                    min_CI_95 = minavatt,
                                    max_CI_95 = maxavatt,
                                    p = p_avatt),
              averagepergroup_mean_attribute = list(observed = avgroupavatt_obs,
                                                    mean_sample = mean(allavgroupavatt),
                                                    sd_sample = sd(allavgroupavatt),
                                                    min_CI_95 = minavgroupavatt,
                                                    max_CI_95 = maxavgroupavatt,
                                                    p = p_avgroupavatt))
  
  return(res)
}


# CUG test for partition: binary attribute (also works for categorical)
CUG_partition_binary <- function(partition_sample, partition_obs, attribute) {

  
  dist<- as.matrix(dist(partition_obs, diag = T, upper = T))
  net_partition_obs <- dist
  net_partition_obs[dist == 0] <- 1
  net_partition_obs[dist != 0] <- 0
  net_partition_obs[is.na(dist)] <- 0
  diag(net_partition_obs) <- 0
  
  net_partition_sample <- list()
  for(s in 1:length(partition_sample)) {
    dist<- as.matrix(dist(partition_sample[[s]], diag = T, upper = T))
    net <- dist
    net[dist == 0] <- 1
    net[dist != 0] <- 0
    net[is.na(dist)] <- 0
    diag(net) <- 0
    net_partition_sample[[s]] <- net
  }
  
  dist <- as.matrix(dist(as.numeric(factor(attribute)), diag = T, upper = T))
  distances <- dist
  distances[dist == 0] <- 1
  distances[dist != 0] <- 0
  diag(distances) <- 0
  
  sameties_obs <- sum(distances * net_partition_obs) / 2
  avgroupsameties_obs <- sum(distances * net_partition_obs) / (2*max(partition_obs, na.rm = T))
  densities_obs <- computedensity(distances, partition_obs)$average
  indsameties_obs <- sum( rowSums(distances * net_partition_obs) > 0 ) 
  proportions_obs <- computeproportion(attribute, partition_obs)$sum
  avproportions_obs <- computeproportion(attribute, partition_obs)$average
  
  allsameties <- rep(0,length(net_partition_sample))
  allavgroupsameties <- rep(0,length(net_partition_sample))
  alldensities  <- rep(0,length(net_partition_sample))
  allindsameties <- rep(0,length(net_partition_sample))
  allproportions <- rep(0,length(net_partition_sample))
  allavproportions <- rep(0,length(net_partition_sample))
  for(s in 1:length(net_partition_sample)){
    allsameties[s] <- sum(distances * net_partition_sample[[s]]) / 2
    allavgroupsameties[s] <- sum(distances * net_partition_sample[[s]]) / (2*max(partition_sample[[s]], na.rm=T))
    alldensities[s] <- computedensity(distances, partition_sample[[s]])$average
    allindsameties[s] <- sum( rowSums(distances * net_partition_sample[[s]]) > 0 )
    allproportions[s] <- computeproportion(attribute, partition_sample[[s]])$sum
    allavproportions[s] <- computeproportion(attribute, partition_sample[[s]])$average
  }
  
  minsame <- quantile(allsameties, probs = 0.05)
  maxsame <- quantile(allsameties, probs = 0.95)
  minavgroupsame <- quantile(allavgroupsameties, probs = 0.05)
  maxavgroupsame <- quantile(allavgroupsameties, probs = 0.95)
  mindensities <- quantile(alldensities, probs = 0.05)
  maxdensities <- quantile(alldensities, probs = 0.95)
  minindsame <- quantile(allindsameties, probs = 0.05)
  maxindsame <- quantile(allindsameties, probs = 0.95)
  minproportions <- quantile(allproportions, probs = 0.05)
  maxproportions <- quantile(allproportions, probs = 0.95)
  minavproportions <- quantile(allavproportions, probs = 0.05)
  maxavproportions <- quantile(allavproportions, probs = 0.95)
  
  
  p1 <- sum(allsameties <= sameties_obs) / length(net_partition_sample)
  p2 <- sum(allsameties >= sameties_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_same <- p2
  } else if(p2 == 0) {
    p_same <- p1
  } else {
    p_same <- 2*min(p1,p2)
  }
  
  p1 <- sum(allavgroupsameties <= avgroupsameties_obs) / length(net_partition_sample)
  p2 <- sum(allavgroupsameties >= avgroupsameties_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_avgroupsame <- p2
  } else if(p2 == 0) {
    p_avgroupsame <- p1
  } else {
    p_avgroupsame <- 2*min(p1,p2)
  }
  
  p1 <- sum(alldensities <= densities_obs) / length(net_partition_sample)
  p2 <- sum(alldensities >= densities_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_densities <- p2
  } else if(p2 == 0) {
    p_densities <- p1
  } else {
    p_densities <- 2*min(p1,p2)
  }
  
  p1 <- sum(allindsameties <= indsameties_obs) / length(net_partition_sample)
  p2 <- sum(allindsameties >= indsameties_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_indsame <- p2
  } else if(p2 == 0) {
    p_indsame <- p1
  } else {
    p_indsame <- 2*min(p1,p2)
  }
  
  p1 <- sum(allproportions <= proportions_obs) / length(partition_sample)
  p2 <- sum(allproportions >= proportions_obs) / length(partition_sample)
  if(p1 == 0) {
    p_proportions <- p2
  } else if(p2 == 0) {
    p_proportions <- p1
  } else {
    p_proportions <- 2*min(p1,p2)
  }
  
  p1 <- sum(allavproportions <= avproportions_obs) / length(partition_sample)
  p2 <- sum(allavproportions >= avproportions_obs) / length(partition_sample)
  if(p1 == 0) {
    p_avproportions <- p2
  } else if(p2 == 0) {
    p_avproportions <- p1
  } else {
    p_avproportions <- 2*min(p1,p2)
  }
  
  res <- list(number_same_ties = list(observed = sameties_obs,
                                      mean_sample = mean(allsameties),
                                      sd_sample = sd(allsameties),
                                      min_CI_95 = minsame,
                                      max_CI_95 = maxsame,
                                      p = p_same),
              averagepergroup_number_same_ties = list(observed = avgroupsameties_obs,
                                                      mean_sample = mean(allavgroupsameties),
                                                      sd_sample = sd(allavgroupsameties),
                                                      min_CI_95 = minavgroupsame,
                                                      max_CI_95 = maxavgroupsame,
                                                      p = p_avgroupsame),
              intratie_densities = list(observed = densities_obs,
                                        mean_sample = mean(alldensities),
                                        sd_sample = sd(alldensities),
                                        min_CI_95 = mindensities,
                                        max_CI_95 = maxdensities,
                                        p = p_densities),
              individual_number_same_ties = list(observed = indsameties_obs,
                                                 mean_sample = mean(allindsameties),
                                                 sd_sample = sd(allindsameties),
                                                 min_CI_95 = minindsame,
                                                 max_CI_95 = maxindsame,
                                                 p = p_indsame),
              sum_proportions = list(observed = proportions_obs,
                                     mean_sample = mean(allproportions),
                                     sd_sample = sd(allproportions),
                                     min_CI_95 = minproportions,
                                     max_CI_95 = maxproportions,
                                     p = p_proportions),
              averagepergroup_proportions = list(observed = avproportions_obs,
                                                 mean_sample = mean(allavproportions),
                                                 sd_sample = sd(allavproportions),
                                                 min_CI_95 = minavproportions,
                                                 max_CI_95 = maxavproportions,
                                                 p = p_avproportions))
  
  return(res)
}

# CUG test for partition: dyadic attribute
CUG_partition_dyadic <- function(partition_sample, partition_obs, net_attribute) {
  
  dist<- as.matrix(dist(partition_obs, diag = T, upper = T))
  net_partition_obs <- dist
  net_partition_obs[dist == 0] <- 1
  net_partition_obs[dist != 0] <- 0
  net_partition_obs[is.na(dist)] <- 0
  diag(net_partition_obs) <- 0
  
  net_partition_sample <- list()
  for(s in 1:length(partition_sample)) {
    dist<- as.matrix(dist(partition_sample[[s]], diag = T, upper = T))
    net <- dist
    net[dist == 0] <- 1
    net[dist != 0] <- 0
    net[is.na(dist)] <- 0
    diag(net) <- 0
    net_partition_sample[[s]] <- net
  }
  
  numties_obs <- sum(net_attribute * net_partition_obs) / 2
  avgroupnumties_obs <- sum(net_attribute * net_partition_obs) / (2*max(partition_obs, na.rm = T))
  densities_obs <- computedensity(net_attribute, partition_obs)$average
  indnumties_obs <- sum( rowSums(net_attribute * net_partition_obs) > 0 ) 
  
  allnumties <- rep(0,length(net_partition_sample))
  allavgroupnumties <- rep(0,length(net_partition_sample))
  alldensities  <- rep(0,length(net_partition_sample))
  allindnumties <- rep(0,length(net_partition_sample))
  for(s in 1:length(net_partition_sample)){
    allnumties[s] <- sum(net_attribute * net_partition_sample[[s]]) / 2
    allavgroupnumties[s] <- sum(net_attribute * net_partition_sample[[s]]) / (2*max(partition_sample[[s]], na.rm=T))
    alldensities[s] <- computedensity(net_attribute, partition_sample[[s]])$average
    allindnumties[s] <- sum( rowSums(net_attribute * net_partition_sample[[s]]) > 0 )
  }
  
  minnum <- quantile(allnumties, probs = 0.05)
  maxnum <- quantile(allnumties, probs = 0.95)
  minavgroupnum <- quantile(allavgroupnumties, probs = 0.05)
  maxavgroupnum <- quantile(allavgroupnumties, probs = 0.95)
  mindensities <- quantile(alldensities, probs = 0.05)
  maxdensities <- quantile(alldensities, probs = 0.95)
  minindnum <- quantile(allindnumties, probs = 0.05)
  maxindnum <- quantile(allindnumties, probs = 0.95)

  
  p1 <- sum(allnumties <= numties_obs) / length(net_partition_sample)
  p2 <- sum(allnumties >= numties_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_num <- p2
  } else if(p2 == 0) {
    p_num <- p1
  } else {
    p_num <- 2*min(p1,p2)
  }
  
  p1 <- sum(allavgroupnumties <= avgroupnumties_obs) / length(net_partition_sample)
  p2 <- sum(allavgroupnumties >= avgroupnumties_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_avgroupnum <- p2
  } else if(p2 == 0) {
    p_avgroupnum <- p1
  } else {
    p_avgroupnum <- 2*min(p1,p2)
  }
  
  p1 <- sum(alldensities <= densities_obs) / length(net_partition_sample)
  p2 <- sum(alldensities >= densities_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_densities <- p2
  } else if(p2 == 0) {
    p_densities <- p1
  } else {
    p_densities <- 2*min(p1,p2)
  }
  
  p1 <- sum(allindnumties <= indnumties_obs) / length(net_partition_sample)
  p2 <- sum(allindnumties >= indnumties_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_indnum <- p2
  } else if(p2 == 0) {
    p_indnum <- p1
  } else {
    p_indnum <- 2*min(p1,p2)
  }
  
  res <- list(number_ties = list(observed = numties_obs,
                                 mean_sample = mean(allnumties),
                                 sd_sample = sd(allnumties),
                                 min_CI_95 = minnum,
                                 max_CI_95 = maxnum,
                                 p = p_num),
              averagepergroup_number_ties = list(observed = avgroupnumties_obs,
                                                 mean_sample = mean(allavgroupnumties),
                                                 sd_sample = sd(allavgroupnumties),
                                                 min_CI_95 = minavgroupnum,
                                                 max_CI_95 = maxavgroupnum,
                                                 p = p_avgroupnum),
              intratie_densities = list(observed = densities_obs,
                                        mean_sample = mean(alldensities),
                                        sd_sample = sd(alldensities),
                                        min_CI_95 = mindensities,
                                        max_CI_95 = maxdensities,
                                        p = p_densities),
              individual_number_ties = list(observed = indnumties_obs,
                                            mean_sample = mean(allindnumties),
                                            sd_sample = sd(allindnumties),
                                            min_CI_95 = minindnum,
                                            max_CI_95 = maxindnum,
                                            p = p_indnum))
  
  return(res)
}

# CUG test for partition: dyadic attribute and a given tye for the partition groups
CUG_partition_dyadic_withtypes <- function(partition_sample, partition_obs, types_obs, net_attribute, type) {
  
  partition_obs2 <- partition_obs
  out_groups <- which(types_obs != type)
  partition_obs2[partition_obs2 %in% out_groups] <- NA
  partition_obs2[!is.na(partition_obs2)] <- order_groupids(partition_obs2[!is.na(partition_obs2)])
  
  dist<- as.matrix(dist(partition_obs2, diag = T, upper = T))
  net_partition_obs2 <- dist
  net_partition_obs2[dist == 0] <- 1
  net_partition_obs2[dist != 0] <- 0
  net_partition_obs2[is.na(dist)] <- 0
  diag(net_partition_obs2) <- 0
  
  partition_sample2 <- list()
  net_partition_sample2 <- list()
  for(s in 1:length(partition_sample)) {
    p2 <- partition_sample[[s]]$partition
    out_groups <- which(partition_sample[[s]]$types != type)
    p2[p2 %in% out_groups] <- NA
    p2[!is.na(p2)] <- order_groupids(p2[!is.na(p2)])
    partition_sample2[[s]] <- p2
    
    dist<- as.matrix(dist(p2, diag = T, upper = T))
    net2 <- dist
    net2[dist == 0] <- 1
    net2[dist != 0] <- 0
    net2[is.na(dist)] <- 0
    diag(net2) <- 0
    net_partition_sample2[[s]] <- net2
  }
  
  numties_obs <- sum(net_attribute * net_partition_obs2) / 2
  avgroupnumties_obs <- sum(net_attribute * net_partition_obs2) / (2*max(partition_obs2, na.rm = T))
  densities_obs <- computedensity(net_attribute, partition_obs2)$average
  indnumties_obs <- sum( rowSums(net_attribute * net_partition_obs2) > 0 ) 
  
  allnumties <- rep(0,length(net_partition_sample))
  allavgroupnumties <- rep(0,length(net_partition_sample))
  alldensities  <- rep(0,length(net_partition_sample))
  allindnumties <- rep(0,length(net_partition_sample))
  for(s in 1:length(net_partition_sample)){
    allnumties[s] <- sum(net_attribute * net_partition_sample2[[s]]) / 2
    allavgroupnumties[s] <- sum(net_attribute * net_partition_sample2[[s]]) / (2*max(partition_sample2[[s]], na.rm=T))
    alldensities[s] <- computedensity(net_attribute, partition_sample2[[s]])$average
    allindnumties[s] <- sum( rowSums(net_attribute * net_partition_sample2[[s]]) > 0 )
  }
  
  minnum <- quantile(allnumties, probs = 0.05)
  maxnum <- quantile(allnumties, probs = 0.95)
  minavgroupnum <- quantile(allavgroupnumties, probs = 0.05)
  maxavgroupnum <- quantile(allavgroupnumties, probs = 0.95)
  mindensities <- quantile(alldensities, probs = 0.05)
  maxdensities <- quantile(alldensities, probs = 0.95)
  minindnum <- quantile(allindnumties, probs = 0.05)
  maxindnum <- quantile(allindnumties, probs = 0.95)
  
  
  p1 <- sum(allnumties <= numties_obs) / length(net_partition_sample)
  p2 <- sum(allnumties >= numties_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_num <- p2
  } else if(p2 == 0) {
    p_num <- p1
  } else {
    p_num <- 2*min(p1,p2)
  }
  
  p1 <- sum(allavgroupnumties <= avgroupnumties_obs) / length(net_partition_sample)
  p2 <- sum(allavgroupnumties >= avgroupnumties_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_avgroupnum <- p2
  } else if(p2 == 0) {
    p_avgroupnum <- p1
  } else {
    p_avgroupnum <- 2*min(p1,p2)
  }
  
  p1 <- sum(alldensities <= densities_obs) / length(net_partition_sample)
  p2 <- sum(alldensities >= densities_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_densities <- p2
  } else if(p2 == 0) {
    p_densities <- p1
  } else {
    p_densities <- 2*min(p1,p2)
  }
  
  p1 <- sum(allindnumties <= indnumties_obs) / length(net_partition_sample)
  p2 <- sum(allindnumties >= indnumties_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_indnum <- p2
  } else if(p2 == 0) {
    p_indnum <- p1
  } else {
    p_indnum <- 2*min(p1,p2)
  }
  
  res <- list(number_ties = list(observed = numties_obs,
                                 mean_sample = mean(allnumties),
                                 sd_sample = sd(allnumties),
                                 min_CI_95 = minnum,
                                 max_CI_95 = maxnum,
                                 p = p_num),
              averagepergroup_number_ties = list(observed = avgroupnumties_obs,
                                                 mean_sample = mean(allavgroupnumties),
                                                 sd_sample = sd(allavgroupnumties),
                                                 min_CI_95 = minavgroupnum,
                                                 max_CI_95 = maxavgroupnum,
                                                 p = p_avgroupnum),
              intratie_densities = list(observed = densities_obs,
                                        mean_sample = mean(alldensities),
                                        sd_sample = sd(alldensities),
                                        min_CI_95 = mindensities,
                                        max_CI_95 = maxdensities,
                                        p = p_densities),
              individual_number_ties = list(observed = indnumties_obs,
                                            mean_sample = mean(allindnumties),
                                            sd_sample = sd(allindnumties),
                                            min_CI_95 = minindnum,
                                            max_CI_95 = maxindnum,
                                            p = p_indnum))
  
  return(res)
}

# CUG test for partition: continuous attribute
CUG_partition_continuous <- function(partition_sample, partition_obs, attribute) {
  
  dist<- as.matrix(dist(partition_obs, diag = T, upper = T))
  net_partition_obs <- dist
  net_partition_obs[dist == 0] <- 1
  net_partition_obs[dist != 0] <- 0
  net_partition_obs[is.na(dist)] <- 0
  diag(net_partition_obs) <- 0
  
  net_partition_sample <- list()
  for(s in 1:length(partition_sample)) {
    dist<- as.matrix(dist(partition_sample[[s]], diag = T, upper = T))
    net <- dist
    net[dist == 0] <- 1
    net[dist != 0] <- 0
    net[is.na(dist)] <- 0
    diag(net) <- 0
    net_partition_sample[[s]] <- net
  }
  
  distances <- as.matrix(dist(attribute, diag = T, upper = T))
  
  diffs_obs <- sum(distances * net_partition_obs)
  avdiffs_obs <- sum(distances * net_partition_obs) / max(partition_obs, na.rm = T)
  icc_obs <- computeicc(attribute, partition_obs)
  inddiffs_obs <- sum( apply(distances * net_partition_obs, 1, FUN=min2) )
  ranges_obs <- computerange(attribute, partition_obs)$sum
  avranges_obs <- computerange(attribute, partition_obs)$average
  
  alldiffs <- rep(0,length(net_partition_sample))
  allavdiffs <- rep(0,length(net_partition_sample))
  allicc <- rep(0,length(net_partition_sample))
  allinddiffs <- rep(0,length(net_partition_sample))
  allranges <- rep(0,length(net_partition_sample))
  allavranges <- rep(0,length(net_partition_sample))
  for(s in 1:length(net_partition_sample)){
    alldiffs[s] <- sum(distances * net_partition_sample[[s]])
    allavdiffs[s] <- sum(distances * net_partition_sample[[s]]) / max(partition_sample[[s]], na.rm = T)
    allicc[s] <- computeicc(attribute, partition_sample[[s]])
    allinddiffs[s] <- sum( apply(distances * net_partition_sample[[s]], 1, FUN=min2) )
    allranges[s] <- computerange(attribute, partition_sample[[s]])$sum
    allavranges[s] <- computerange(attribute, partition_sample[[s]])$average
  }
  
  mindiffs <- quantile(alldiffs, probs = 0.05)
  maxdiffs <- quantile(alldiffs, probs = 0.95)
  minavdiffs <- quantile(allavdiffs, probs = 0.05)
  maxavdiffs <- quantile(allavdiffs, probs = 0.95)
  minicc <- quantile(allicc, probs = 0.05)
  maxicc <- quantile(allicc, probs = 0.95)
  mininddiffs <- quantile(allinddiffs, probs = 0.05)
  maxinddiffs <- quantile(allinddiffs, probs = 0.95)
  minranges <- quantile(allranges, probs = 0.05)
  maxranges <- quantile(allranges, probs = 0.95)
  minavranges <- quantile(allavranges, probs = 0.05)
  maxavranges <- quantile(allavranges, probs = 0.95)
  
  
  p1 <- sum(alldiffs <= diffs_obs) / length(net_partition_sample)
  p2 <- sum(alldiffs >= diffs_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_diffs <- p2
  } else if(p2 == 0) {
    p_diffs <- p1
  } else {
    p_diffs <- 2*min(p1,p2)
  }
  
  p1 <- sum(allavdiffs <= avdiffs_obs) / length(net_partition_sample)
  p2 <- sum(allavdiffs >= avdiffs_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_avdiffs <- p2
  } else if(p2 == 0) {
    p_avdiffs <- p1
  } else {
    p_avdiffs <- 2*min(p1,p2)
  }
  
  p1 <- sum(allicc <= icc_obs) / length(partition_sample)
  p2 <- sum(allicc >= icc_obs) / length(partition_sample)
  if(p1 == 0) {
    p_icc <- p2
  } else if(p2 == 0) {
    p_icc <- p1
  } else {
    p_icc <- 2*min(p1,p2)
  }
  
  p1 <- sum(allinddiffs <= inddiffs_obs) / length(net_partition_sample)
  p2 <- sum(allinddiffs >= inddiffs_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_inddiffs <- p2
  } else if(p2 == 0) {
    p_inddiffs <- p1
  } else {
    p_inddiffs <- 2*min(p1,p2)
  }
  
  p1 <- sum(allranges <= ranges_obs) / length(partition_sample)
  p2 <- sum(allranges >= ranges_obs) / length(partition_sample)
  if(p1 == 0) {
    p_ranges <- p2
  } else if(p2 == 0) {
    p_ranges <- p1
  } else {
    p_ranges <- 2*min(p1,p2)
  }
  
  p1 <- sum(allavranges <= avranges_obs) / length(partition_sample)
  p2 <- sum(allavranges >= avranges_obs) / length(partition_sample)
  if(p1 == 0) {
    p_avranges <- p2
  } else if(p2 == 0) {
    p_avranges <- p1
  } else {
    p_avranges <- 2*min(p1,p2)
  }
  
  res <- list(sum_absolute_differences = list(observed = diffs_obs,
                                              mean_sample = mean(alldiffs),
                                              sd_sample = sd(alldiffs),
                                              min_CI_95 = mindiffs,
                                              max_CI_95 = maxdiffs,
                                              p = p_diffs),
              averagepergroup_absolute_differences = list(observed = avdiffs_obs,
                                                          mean_sample = mean(allavdiffs),
                                                          sd_sample = sd(allavdiffs),
                                                          min_CI_95 = minavdiffs,
                                                          max_CI_95 = maxavdiffs,
                                                          p = p_avdiffs),
              icc = list(observed = icc_obs,
                         mean_sample = mean(allicc),
                         sd_sample = sd(allicc),
                         min_CI_95 = minicc,
                         max_CI_95 = maxicc,
                         p = p_icc),
              individual_absolute_differences = list(observed = inddiffs_obs,
                                                     mean_sample = mean(allinddiffs),
                                                     sd_sample = sd(allinddiffs),
                                                     min_CI_95 = mininddiffs,
                                                     max_CI_95 = maxinddiffs,
                                                     p = p_inddiffs),
              sum_ranges = list(observed = ranges_obs,
                                mean_sample = mean(allranges),
                                sd_sample = sd(allranges),
                                min_CI_95 = minranges,
                                max_CI_95 = maxranges,
                                p = p_ranges),
              averagepergroup_ranges = list(observed = avranges_obs,
                                            mean_sample = mean(allavranges),
                                            sd_sample = sd(allavranges),
                                            min_CI_95 = minavranges,
                                            max_CI_95 = maxavranges,
                                            p = p_avranges))
  
  return(res)
}


# CUG test for partition: continuous attributes with threshold for the difference
CUG_partition_continuous_threshold <- function(partition_sample, partition_obs, attribute, threshold) {
  
  
  dist<- as.matrix(dist(partition_obs, diag = T, upper = T))
  net_partition_obs <- dist
  net_partition_obs[dist == 0] <- 1
  net_partition_obs[dist != 0] <- 0
  net_partition_obs[is.na(dist)] <- 0
  diag(net_partition_obs) <- 0
  
  net_partition_sample <- list()
  for(s in 1:length(partition_sample)) {
    dist<- as.matrix(dist(partition_sample[[s]], diag = T, upper = T))
    net <- dist
    net[dist == 0] <- 1
    net[dist != 0] <- 0
    net[is.na(dist)] <- 0
    diag(net) <- 0
    net_partition_sample[[s]] <- net
  }
  
  dist <- as.matrix(dist(attribute, diag = T, upper = T))
  distances <- dist
  distances[dist <= threshold] <- 1
  distances[dist > threshold] <- 0
  diag(distances) <- 0
  
  sameties_obs <- sum(distances * net_partition_obs) / 2
  avgroupsameties_obs <- sum(distances * net_partition_obs) / (2*max(partition_obs, na.rm = T))
  densities_obs <- computedensity(distances, partition_obs)$average
  indsameties_obs <- sum( rowSums(distances * net_partition_obs) > 0 ) 

  allsameties <- rep(0,length(net_partition_sample))
  allavgroupsameties <- rep(0,length(net_partition_sample))
  alldensities  <- rep(0,length(net_partition_sample))
  allindsameties <- rep(0,length(net_partition_sample))
  for(s in 1:length(net_partition_sample)){
    allsameties[s] <- sum(distances * net_partition_sample[[s]]) / 2
    allavgroupsameties[s] <- sum(distances * net_partition_sample[[s]]) / (2*max(partition_sample[[s]], na.rm = T))
    alldensities[s] <- computedensity(distances, partition_sample[[s]])$average
    allindsameties[s] <- sum( rowSums(distances * net_partition_sample[[s]]) > 0 )
  }
  
  minsame <- quantile(allsameties, probs = 0.05)
  maxsame <- quantile(allsameties, probs = 0.95)
  minavgroupsame <- quantile(allavgroupsameties, probs = 0.05)
  maxavgroupsame <- quantile(allavgroupsameties, probs = 0.95)
  mindensities <- quantile(alldensities, probs = 0.05)
  maxdensities <- quantile(alldensities, probs = 0.95)
  minindsame <- quantile(allindsameties, probs = 0.05)
  maxindsame <- quantile(allindsameties, probs = 0.95)
  
  
  p1 <- sum(allsameties <= sameties_obs) / length(net_partition_sample)
  p2 <- sum(allsameties >= sameties_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_same <- p2
  } else if(p2 == 0) {
    p_same <- p1
  } else {
    p_same <- 2*min(p1,p2)
  }
  
  p1 <- sum(allavgroupsameties <= avgroupsameties_obs) / length(net_partition_sample)
  p2 <- sum(allavgroupsameties >= avgroupsameties_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_avgroupsame <- p2
  } else if(p2 == 0) {
    p_avgroupsame <- p1
  } else {
    p_avgroupsame <- 2*min(p1,p2)
  }
  
  p1 <- sum(alldensities <= densities_obs) / length(partition_sample)
  p2 <- sum(alldensities >= densities_obs) / length(partition_sample)
  if(p1 == 0) {
    p_densities <- p2
  } else if(p2 == 0) {
    p_densities <- p1
  } else {
    p_densities <- 2*min(p1,p2)
  }
  
  p1 <- sum(allindsameties <= indsameties_obs) / length(net_partition_sample)
  p2 <- sum(allindsameties >= indsameties_obs) / length(net_partition_sample)
  if(p1 == 0) {
    p_indsame <- p2
  } else if(p2 == 0) {
    p_indsame <- p1
  } else {
    p_indsame <- 2*min(p1,p2)
  }

  
  res <- list(number_same_ties = list(observed = sameties_obs,
                                      mean_sample = mean(allsameties),
                                      sd_sample = sd(allsameties),
                                      min_CI_95 = minsame,
                                      max_CI_95 = maxsame,
                                      p = p_same),
              averagepergroup_number_same_ties = list(observed = avgroupsameties_obs,
                                                      mean_sample = mean(allavgroupsameties),
                                                      sd_sample = sd(allavgroupsameties),
                                                      min_CI_95 = minavgroupsame,
                                                      max_CI_95 = maxavgroupsame,
                                                      p = p_avgroupsame),
              intratie_densities = list(observed = densities_obs,
                                        mean_sample = mean(alldensities),
                                        sd_sample = sd(alldensities),
                                        min_CI_95 = mindensities,
                                        max_CI_95 = maxdensities,
                                        p = p_densities),
              individual_number_same_ties = list(observed = indsameties_obs,
                                                 mean_sample = mean(allindsameties),
                                                 sd_sample = sd(allindsameties),
                                                 min_CI_95 = minindsame,
                                                 max_CI_95 = maxindsame,
                                                 p = p_indsame))
  
  return(res)
}



# CUG test for multiple partitions: binary attribute (also works for categorical)
CUG_partitions_binary <- function(partitions_sample, partitions_obs, attribute) {
  
  n_obs <- ncol(partition_obs)
  
  nets_partitions_obs <- list()
  for(o in 1:n_obs){
    dist<- as.matrix(dist(partitions_obs[,o], diag = T, upper = T))
    net <- dist
    net[dist == 0] <- 1
    net[dist != 0] <- 0
    net[is.na(dist)] <- 0
    diag(net) <- 0
    nets_partitions_obs[[o]] <- net
  }
  
  nets_partitions_sample <- list()
  for(s in 1:length(partitions_sample)) {
    nets_partitions_sample[[s]] <- list()
    for(o in 1:n_obs){
      dist<- as.matrix(dist(partitions_sample[[s]][,o], diag = T, upper = T))
      net <- dist
      net[dist == 0] <- 1
      net[dist != 0] <- 0
      net[is.na(dist)] <- 0
      diag(net) <- 0
      nets_partitions_sample[[s]][[o]] <- net
    }
  }
  
  dist <- as.matrix(dist(as.numeric(factor(attribute)), diag = T, upper = T))
  distances <- dist
  distances[dist == 0] <- 1
  distances[dist != 0] <- 0
  diag(distances) <- 0
  
  sameties_obs <- rep(0,n_obs)
  avgroupsameties_obs <- rep(0,n_obs)
  densities_obs <- rep(0,n_obs)
  indsameties_obs <- rep(0,n_obs)
  proportions_obs <- rep(0,n_obs)
  avproportions_obs <- rep(0,n_obs)
  for(o in 1:n_obs){
    sameties_obs[o] <- sum(distances * nets_partitions_obs[[o]]) / 2
    avgroupsameties_obs[o] <- sum(distances * nets_partitions_obs[[o]]) / (2*max(partitions_obs[,o], na.rm = T))
    densities_obs[o] <- computedensity(distances, partitions_obs[,o])$average
    indsameties_obs[o] <- sum( rowSums(distances * nets_partitions_obs[o]) > 0 ) 
    proportions_obs[o] <- computeproportion(attribute, partitions_obs[,o])$sum
    avproportions_obs[o] <- computeproportion(attribute, partitions_obs[,o])$average
  }
  
  allsameties <- matrix(0,length(nets_partitions_sample),n_obs)
  allavgroupsameties <- matrix(0,length(nets_partitions_sample),n_obs)
  alldensities  <- matrix(0,length(nets_partitions_sample),n_obs)
  allindsameties <- matrix(0,length(nets_partitions_sample),n_obs)
  allproportions <- matrix(0,length(nets_partitions_sample),n_obs)
  allavproportions <- matrix(0,length(nets_partitions_sample),n_obs)
  for(s in 1:length(nets_partitions_sample)){
    for(o in 1:n_obs){
      allsameties[s,o] <- sum(distances * nets_partitions_sample[[s]][[o]]) / 2
      allavgroupsameties[s,o] <- sum(distances * nets_partitions_sample[[s]][[o]]) / (2*max(partitions_sample[[s]][,o], na.rm = T))
      alldensities[s,o] <- computedensity(distances, partitions_sample[[s]][,o])$average
      allindsameties[s,o] <- sum( rowSums(distances * nets_partitions_sample[[s]][[o]]) > 0 )
      allproportions[s,o] <- computeproportion(attribute, partitions_sample[[s]][,o])$sum
      allavproportions[s,o] <- computeproportion(attribute, partitions_sample[[s]][,o])$average
    }
  }
  
  minsame <- apply(allsameties,2, FUN = function(x){quantile(x, probs = 0.05)})
  maxsame <- apply(allsameties,2, FUN = function(x){quantile(x, probs = 0.95)})
  minavgroupsame <- apply(allavgroupsameties,2, FUN = function(x){quantile(x, probs = 0.05)})
  maxavgroupsame <- apply(allavgroupsameties,2, FUN = function(x){quantile(x, probs = 0.95)})
  mindensities <- apply(alldensities,2, FUN = function(x){quantile(x, probs = 0.05)})
  maxdensities <- apply(alldensities,2, FUN = function(x){quantile(x, probs = 0.95)})
  minindsame <- apply(allindsameties,2, FUN = function(x){quantile(x, probs = 0.05)})
  maxindsame <- apply(allindsameties,2, FUN = function(x){quantile(x, probs = 0.95)})
  minproportions <- apply(allproportions,2, FUN = function(x){quantile(x, probs = 0.05)})
  maxproportions <- apply(allproportions,2, FUN = function(x){quantile(x, probs = 0.95)})
  minavproportions <- apply(allavproportions,2, FUN = function(x){quantile(x, probs = 0.05)})
  maxavproportions <- apply(allavproportions,2, FUN = function(x){quantile(x, probs = 0.95)})
  
  p_same <- rep(0,n_obs)
  p_avgroupsame <- rep(0,n_obs)
  p_densities <- rep(0,n_obs)
  p_indsame <- rep(0,n_obs)
  p_proportions <- rep(0,n_obs)
  p_avproportions <- rep(0,n_obs)
  
  for(o in 1:n_obs){
    p1 <- sum(allsameties[,o] <= sameties_obs[o]) / length(nets_partitions_sample)
    p2 <- sum(allsameties[,o] >= sameties_obs[o]) / length(nets_partitions_sample)
    if(p1 == 0) {
      p_same[o] <- p2
    } else if(p2 == 0) {
      p_same[o] <- p1
    } else {
      p_same[o] <- 2*min(p1,p2)
    }
  
  p1 <- sum(allavgroupsameties[,o] <= avgroupsameties_obs[o]) / length(nets_partitions_sample)
  p2 <- sum(allavgroupsameties[,o] >= avgroupsameties_obs[o]) / length(nets_partitions_sample)
  if(p1 == 0) {
    p_avgroupsame[o] <- p2
  } else if(p2 == 0) {
    p_avgroupsame[o] <- p1
  } else {
    p_avgroupsame[o] <- 2*min(p1,p2)
  }
  
  p1 <- sum(alldensities[,o] <= densities_obs[o]) / length(nets_partitions_sample)
  p2 <- sum(alldensities[,o] >= densities_obs[o]) / length(nets_partitions_sample)
  if(p1 == 0) {
    p_densities[o] <- p2
  } else if(p2 == 0) {
    p_densities[o] <- p1
  } else {
    p_densities[o] <- 2*min(p1,p2)
  }
  
  p1 <- sum(allindsameties[,o] <= indsameties_obs[o]) / length(nets_partitions_sample)
  p2 <- sum(allindsameties[,o] >= indsameties_obs[o]) / length(nets_partitions_sample)
  if(p1 == 0) {
    p_indsame[o] <- p2
  } else if(p2 == 0) {
    p_indsame[o] <- p1
  } else {
    p_indsame[o] <- 2*min(p1,p2)
  }
  
  p1 <- sum(allproportions[,o] <= proportions_obs[o]) / length(partitions_sample)
  p2 <- sum(allproportions[,o] >= proportions_obs[o]) / length(partitions_sample)
  if(p1 == 0) {
    p_proportions[o] <- p2
  } else if(p2 == 0) {
    p_proportions[o] <- p1
  } else {
    p_proportions[o] <- 2*min(p1,p2)
  }
  
  p1 <- sum(allavproportions[,o] <= avproportions_obs[o]) / length(partitions_sample)
  p2 <- sum(allavproportions[,o] >= avproportions_obs[o]) / length(partitions_sample)
  if(p1 == 0) {
    p_avproportions[o] <- p2
  } else if(p2 == 0) {
    p_avproportions[o] <- p1
  } else {
    p_avproportions[o] <- 2*min(p1,p2)
  }
  }
  
  
  res <- list(number_same_ties = list(observed = sameties_obs,
                                      mean_sample = apply(allsameties,2,mean),
                                      sd_sample = apply(allsameties,2,sd),
                                      min_CI_95 = minsame,
                                      max_CI_95 = maxsame,
                                      p = p_same),
              averagepergroup_number_same_ties = list(observed = avgroupsameties_obs,
                                                      mean_sample = apply(allavgroupsameties,2,mean),
                                                      sd_sample = apply(allavgroupsameties,2,sd),
                                                      min_CI_95 = minavgroupsame,
                                                      max_CI_95 = maxavgroupsame,
                                                      p = p_avgroupsame),
              intratie_densities = list(observed = densities_obs,
                                        mean_sample = apply(alldensities,2,mean),
                                        sd_sample = apply(alldensities,2,sd),
                                        min_CI_95 = mindensities,
                                        max_CI_95 = maxdensities,
                                        p = p_densities),
              individual_number_same_ties = list(observed = indsameties_obs,
                                                 mean_sample = apply(allindsameties,2,mean),
                                                 sd_sample = apply(allindsameties,2,sd),
                                                 min_CI_95 = minindsame,
                                                 max_CI_95 = maxindsame,
                                                 p = p_indsame),
              sum_proportions = list(observed = proportions_obs,
                                     mean_sample = apply(allproportions,2,mean),
                                     sd_sample = apply(allproportions,2,sd),
                                     min_CI_95 = minproportions,
                                     max_CI_95 = maxproportions,
                                     p = p_proportions),
              averagepergroup_proportions = list(observed = avproportions_obs,
                                                 mean_sample = apply(allavproportions,2,mean),
                                                 sd_sample = apply(allavproportions,2,sd),
                                                 min_CI_95 = minavproportions,
                                                 max_CI_95 = maxavproportions,
                                                 p = p_avproportions))
  
  return(res)
}


# CUG test for multiple partitions: dyadic attribute
CUG_partitions_dyadic <- function(partitions_sample, partitions_obs, net_attribute) {
  
  n_obs <- ncol(partition_obs)
  
  nets_partitions_obs <- list()
  for(o in 1:n_obs){
    dist<- as.matrix(dist(partitions_obs[,o], diag = T, upper = T))
    net <- dist
    net[dist == 0] <- 1
    net[dist != 0] <- 0
    net[is.na(dist)] <- 0
    diag(net) <- 0
    nets_partitions_obs[[o]] <- net
  }
  
  nets_partitions_sample <- list()
  for(s in 1:length(partitions_sample)) {
    nets_partitions_sample[[s]] <- list()
    for(o in 1:n_obs){
      dist<- as.matrix(dist(partitions_sample[[s]][,o], diag = T, upper = T))
      net <- dist
      net[dist == 0] <- 1
      net[dist != 0] <- 0
      net[is.na(dist)] <- 0
      diag(net) <- 0
      nets_partitions_sample[[s]][[o]] <- net
    }
  }
  
  numties_obs <- rep(0,n_obs)
  avgroupnumties_obs <- rep(0,n_obs)
  densities_obs <- rep(0,n_obs)
  indnumties_obs <- rep(0,n_obs)
  for(o in 1:n_obs){
    numties_obs[o] <- sum(net_attribute * nets_partitions_obs[[o]]) / 2
    avgroupnumties_obs[o] <- sum(net_attribute * nets_partitions_obs[[o]]) / (2*max(partitions_obs[,o], na.rm = T))
    densities_obs[o] <- computedensity(net_attribute, partitions_obs[,o])$average
    indnumties_obs[o] <- sum( rowSums(net_attribute * nets_partitions_obs[o]) > 0 ) 
  }
  
  allnumties <- matrix(0,length(nets_partitions_sample),n_obs)
  allavgroupnumties <- matrix(0,length(nets_partitions_sample),n_obs)
  alldensities  <- matrix(0,length(nets_partitions_sample),n_obs)
  allindnumties <- matrix(0,length(nets_partitions_sample),n_obs)
  for(s in 1:length(nets_partitions_sample)){
    for(o in 1:n_obs){
      allnumties[s,o] <- sum(net_attribute * nets_partitions_sample[[s]][[o]]) / 2
      allavgroupnumties[s,o] <- sum(net_attribute * nets_partitions_sample[[s]][[o]]) / (2*max(partitions_sample[[s]][,o], na.rm = T))
      alldensities[s,o] <- computedensity(net_attribute, partitions_sample[[s]][,o])$average
      allindnumties[s,o] <- sum( rowSums(net_attribute * nets_partitions_sample[[s]][[o]]) > 0 )
    }
  }
  
  minnum <- apply(allnumties,2, FUN = function(x){quantile(x, probs = 0.05)})
  maxnum <- apply(allnumties,2, FUN = function(x){quantile(x, probs = 0.95)})
  minavgroupnum <- apply(allavgroupnumties,2, FUN = function(x){quantile(x, probs = 0.05)})
  maxavgroupnum <- apply(allavgroupnumties,2, FUN = function(x){quantile(x, probs = 0.95)})
  mindensities <- apply(alldensities,2, FUN = function(x){quantile(x, probs = 0.05)})
  maxdensities <- apply(alldensities,2, FUN = function(x){quantile(x, probs = 0.95)})
  minindnum <- apply(allindnumties,2, FUN = function(x){quantile(x, probs = 0.05)})
  maxindnum <- apply(allindnumties,2, FUN = function(x){quantile(x, probs = 0.95)})
  
  p_num <- rep(0,n_obs)
  p_avgroupnum <- rep(0,n_obs)
  p_densities <- rep(0,n_obs)
  p_indnum <- rep(0,n_obs)
  
  for(o in 1:n_obs){
    p1 <- sum(allnumties[,o] <= numties_obs[o]) / length(nets_partitions_sample)
    p2 <- sum(allnumties[,o] >= numties_obs[o]) / length(nets_partitions_sample)
    if(p1 == 0) {
      p_num[o] <- p2
    } else if(p2 == 0) {
      p_num[o] <- p1
    } else {
      p_num[o] <- 2*min(p1,p2)
    }
    
    p1 <- sum(allavgroupnumties[,o] <= avgroupnumties_obs[o]) / length(nets_partitions_sample)
    p2 <- sum(allavgroupnumties[,o] >= avgroupnumties_obs[o]) / length(nets_partitions_sample)
    if(p1 == 0) {
      p_avgroupnum[o] <- p2
    } else if(p2 == 0) {
      p_avgroupnum[o] <- p1
    } else {
      p_avgroupnum[o] <- 2*min(p1,p2)
    }
    
    p1 <- sum(alldensities[,o] <= densities_obs[o]) / length(nets_partitions_sample)
    p2 <- sum(alldensities[,o] >= densities_obs[o]) / length(nets_partitions_sample)
    if(p1 == 0) {
      p_densities[o] <- p2
    } else if(p2 == 0) {
      p_densities[o] <- p1
    } else {
      p_densities[o] <- 2*min(p1,p2)
    }
    
    p1 <- sum(allindnumties[,o] <= indnumties_obs[o]) / length(nets_partitions_sample)
    p2 <- sum(allindnumties[,o] >= indnumties_obs[o]) / length(nets_partitions_sample)
    if(p1 == 0) {
      p_indnum[o] <- p2
    } else if(p2 == 0) {
      p_indnum[o] <- p1
    } else {
      p_indnum[o] <- 2*min(p1,p2)
    }
  }
  
  
  res <- list(number_ties = list(observed = numties_obs,
                                 mean_sample = apply(allnumties,2,mean),
                                 sd_sample = apply(allnumties,2,sd),
                                 min_CI_95 = minnum,
                                 max_CI_95 = maxnum,
                                 p = p_num),
              averagepergroup_number_ties = list(observed = avgroupnumties_obs,
                                                 mean_sample = apply(allavgroupnumties,2,mean),
                                                 sd_sample = apply(allavgroupnumties,2,sd),
                                                 min_CI_95 = minavgroupnum,
                                                 max_CI_95 = maxavgroupnum,
                                                 p = p_avgroupnum),
              intratie_densities = list(observed = densities_obs,
                                        mean_sample = apply(alldensities,2,mean),
                                        sd_sample = apply(alldensities,2,sd),
                                        min_CI_95 = mindensities,
                                        max_CI_95 = maxdensities,
                                        p = p_densities),
              individual_number_ties = list(observed = indnumties_obs,
                                            mean_sample = apply(allindnumties,2,mean),
                                            sd_sample = apply(allindnumties,2,sd),
                                            min_CI_95 = minindnum,
                                            max_CI_95 = maxindnum,
                                            p = p_indnum))
  
  return(res)
}


# CUG test for multiple partitions: continuous attribute
CUG_partitions_continuous <- function(partitions_sample, partitions_obs, attribute) {
  
  n_obs <- ncol(partitions_obs)
  
  nets_partitions_obs <- list()
  for(o in 1:n_obs){
    dist<- as.matrix(dist(partitions_obs[,o], diag = T, upper = T))
    net <- dist
    net[dist == 0] <- 1
    net[dist != 0] <- 0
    net[is.na(dist)] <- 0
    diag(net) <- 0
    nets_partitions_obs[[o]] <- net
  }
  
  nets_partitions_sample <- list()
  for(s in 1:length(partitions_sample)) {
    nets_partitions_sample[[s]] <- list()
    for(o in 1:n_obs){
      dist<- as.matrix(dist(partitions_sample[[s]][,o], diag = T, upper = T))
      net <- dist
      net[dist == 0] <- 1
      net[dist != 0] <- 0
      net[is.na(dist)] <- 0
      diag(net) <- 0
      nets_partitions_sample[[s]][[o]] <- net
    }
  }
  
  distances <- as.matrix(dist(attribute, diag = T, upper = T))
  
  diffs_obs <- rep(0,n_obs)
  avdiffs_obs <- rep(0,n_obs)
  icc_obs <- rep(0,n_obs)
  inddiffs_obs <- rep(0,n_obs)
  ranges_obs <- rep(0,n_obs)
  avranges_obs <- rep(0,n_obs)
  for(o in 1:n_obs){
    diffs_obs[o] <- sum(distances * nets_partitions_obs[[o]])
    avdiffs_obs[o] <- sum(distances * nets_partitions_obs[[o]]) / max(partitions_obs[,o], na.rm = T)
    icc_obs[o] <- computeicc(attribute, partitions_obs[,o])
    inddiffs_obs[o] <- sum( apply(distances * nets_partitions_obs[[o]], 1, FUN=min2) )
    ranges_obs[o] <- computerange(attribute, partitions_obs[,o])$sum
    avranges_obs[o] <- computerange(attribute, partitions_obs[,o])$average
  }
  
  alldiffs <- matrix(0,length(net_partition_sample),n_obs)
  allavdiffs <- matrix(0,length(net_partition_sample),n_obs)
  allicc <- matrix(0,length(net_partition_sample),n_obs)
  allinddiffs <- matrix(0,length(net_partition_sample),n_obs)
  allranges <- matrix(0,length(net_partition_sample),n_obs)
  allavranges <- matrix(0,length(net_partition_sample),n_obs)
  for(s in 1:length(net_partition_sample)){
    for(o in 1:n_obs){
      alldiffs[s,o] <- sum(distances * nets_partitions_sample[[s]][[o]])
      allavdiffs[s,o] <- sum(distances * nets_partitions_sample[[s]][[o]]) / max(partitions_sample[[s]][,o], na.rm = T)
      allicc[s,o] <- computeicc(attribute, partitions_sample[[s]][,o])
      allinddiffs[s,o] <- sum( apply(distances * nets_partitions_sample[[s]][[o]], 1,FUN=min2) )
      allranges[s,o] <- computerange(attribute, partitions_sample[[s]][,o])$sum
      allavranges[s,o] <- computerange(attribute, partitions_sample[[s]][,o])$average
    }
  }
  
  mindiffs <- apply(alldiffs,2,FUN=function(x){quantile(x, probs = 0.05)})
  maxdiffs <- apply(alldiffs,2,FUN=function(x){quantile(x, probs = 0.95)})
  minavdiffs <- apply(allavdiffs,2,FUN=function(x){quantile(x, probs = 0.05)})
  maxavdiffs <- apply(allavdiffs,2,FUN=function(x){quantile(x, probs = 0.95)})
  minicc <- apply(allicc,2,FUN=function(x){quantile(x, probs = 0.05)})
  maxicc <- apply(allicc,2,FUN=function(x){quantile(x, probs = 0.95)})
  mininddiffs <- apply(allinddiffs,2,FUN=function(x){quantile(x, probs = 0.05)})
  maxinddiffs <- apply(allinddiffs,2,FUN=function(x){quantile(x, probs = 0.95)})
  minranges <- apply(allranges,2,FUN=function(x){quantile(x, probs = 0.05)})
  maxranges <- apply(allranges,2,FUN=function(x){quantile(x, probs = 0.95)})
  minavranges <- apply(allavranges,2,FUN=function(x){quantile(x, probs = 0.05)})
  maxavranges <- apply(allavranges,2,FUN=function(x){quantile(x, probs = 0.95)})
  
  p_diffs <- rep(0,n_obs)
  p_avdiffs <- rep(0,n_obs)
  p_icc <- rep(0,n_obs)
  p_inddiffs <- rep(0,n_obs)
  p_ranges <- rep(0,n_obs)
  p_avranges <- rep(0,n_obs)
  for(o in 1:n_obs){
    p1 <- sum(alldiffs[,o] <= diffs_obs[o]) / length(nets_partitions_sample)
    p2 <- sum(alldiffs[,o] >= diffs_obs[o]) / length(nets_partitions_sample)
    if(p1 == 0) {
      p_diffs[o] <- p2
    } else if(p2 == 0) {
      p_diffs[o] <- p1
    } else {
      p_diffs[o] <- 2*min(p1,p2)
    }
    
    p1 <- sum(allavdiffs[,o] <= avdiffs_obs[o]) / length(nets_partitions_sample)
    p2 <- sum(allavdiffs[,o] >= avdiffs_obs[o]) / length(nets_partitions_sample)
    if(p1 == 0) {
      p_avdiffs[o] <- p2
    } else if(p2 == 0) {
      p_avdiffs[o] <- p1
    } else {
      p_avdiffs[o] <- 2*min(p1,p2)
    }
    
    p1 <- sum(allicc[,o] <= icc_obs[o]) / length(partitions_sample)
    p2 <- sum(allicc[,o] >= icc_obs[o]) / length(partitions_sample)
    if(p1 == 0) {
      p_icc[o] <- p2
    } else if(p2 == 0) {
      p_icc[o] <- p1
    } else {
      p_icc[o] <- 2*min(p1,p2)
    }
    
    p1 <- sum(allinddiffs[,o] <= inddiffs_obs[o]) / length(nets_partitions_sample)
    p2 <- sum(allinddiffs[,o] >= inddiffs_obs[o]) / length(nets_partitions_sample)
    if(p1 == 0) {
      p_inddiffs[o] <- p2
    } else if(p2 == 0) {
      p_inddiffs[o] <- p1
    } else {
      p_inddiffs[o] <- 2*min(p1,p2)
    }
    
    p1 <- sum(allranges[,o] <= ranges_obs[o]) / length(partitions_sample)
    p2 <- sum(allranges[,o] >= ranges_obs[o]) / length(partitions_sample)
    if(p1 == 0) {
      p_ranges[o] <- p2
    } else if(p2 == 0) {
      p_ranges[o] <- p1
    } else {
      p_ranges[o] <- 2*min(p1,p2)
    }
    
    p1 <- sum(allavranges[,o] <= avranges_obs[o]) / length(partitions_sample)
    p2 <- sum(allavranges[,o] >= avranges_obs[o]) / length(partitions_sample)
    if(p1 == 0) {
      p_avranges[o] <- p2
    } else if(p2 == 0) {
      p_avranges[o] <- p1
    } else {
      p_avranges[o] <- 2*min(p1,p2)
    }
  }
  
  res <- list(sum_absolute_differences = list(observed = diffs_obs,
                                              mean_sample = apply(alldiffs,2,mean),
                                              sd_sample = apply(alldiffs,2,sd),
                                              min_CI_95 = mindiffs,
                                              max_CI_95 = maxdiffs,
                                              p = p_diffs),
              averagepergroup_absolute_differences = list(observed = avdiffs_obs,
                                                          mean_sample = apply(allavdiffs,2,mean),
                                                          sd_sample = apply(allavdiffs,2,sd),
                                                          min_CI_95 = minavdiffs,
                                                          max_CI_95 = maxavdiffs,
                                                          p = p_avdiffs),
              icc = list(observed = icc_obs,
                         mean_sample = apply(allicc,2,mean),
                         sd_sample = apply(allicc,2,sd),
                         min_CI_95 = minicc,
                         max_CI_95 = maxicc,
                         p = p_icc),
              individual_absolute_differences = list(observed = inddiffs_obs,
                                                     mean_sample = apply(allinddiffs,2,mean),
                                                     sd_sample = apply(allinddiffs,2,sd),
                                                     min_CI_95 = mininddiffs,
                                                     max_CI_95 = maxinddiffs,
                                                     p = p_inddiffs),
              sum_ranges = list(observed = ranges_obs,
                                mean_sample = apply(allranges,2,mean),
                                sd_sample = apply(allranges,2,sd),
                                min_CI_95 = minranges,
                                max_CI_95 = maxranges,
                                p = p_ranges),
              averagepergroup_ranges = list(observed = avranges_obs,
                                            mean_sample = apply(allavranges,2,mean),
                                            sd_sample = apply(allavranges,2,sd),
                                            min_CI_95 = minavranges,
                                            max_CI_95 = maxavranges,
                                            p = p_avranges))
  
  return(res)
}


# CUG test for multiple partitions: continuous attribute with threshold
CUG_partitions_continuous_threshold <- function(partitions_sample, partitions_obs, attribute, threshold) {
  
  n_obs <- ncol(partitions_obs)
  
  nets_partitions_obs <- list()
  for(o in 1:n_obs){
    dist<- as.matrix(dist(partitions_obs[,o], diag = T, upper = T))
    net <- dist
    net[dist == 0] <- 1
    net[dist != 0] <- 0
    net[is.na(dist)] <- 0
    diag(net) <- 0
    nets_partitions_obs[[o]] <- net
  }
  
  nets_partitions_sample <- list()
  for(s in 1:length(partitions_sample)) {
    nets_partitions_sample[[s]] <- list()
    for(o in 1:n_obs){
      dist<- as.matrix(dist(partitions_sample[[s]][,o], diag = T, upper = T))
      net <- dist
      net[dist == 0] <- 1
      net[dist != 0] <- 0
      net[is.na(dist)] <- 0
      diag(net) <- 0
      nets_partitions_sample[[s]][[o]] <- net
    }
  }
  
  dist <- as.matrix(dist(attribute, diag = T, upper = T))
  distances <- dist
  distances[dist <= threshold] <- 1
  distances[dist > threshold] <- 0
  diag(distances) <- 0
  
  sameties_obs <- rep(0,n_obs)
  avgroupsameties_obs <- rep(0,n_obs)
  densities_obs <- rep(0,n_obs)
  indsameties_obs <- rep(0,n_obs)
  for(o in 1:n_obs){
    sameties_obs[o] <- sum(distances * nets_partitions_obs[[o]]) / 2
    avgroupsameties_obs[o] <- sum(distances * nets_partitions_obs[[o]]) / (2*max(partitions_obs[,o], na.rm = T))
    densities_obs[o] <- computedensity(distances, partitions_obs[,o])$average
    indsameties_obs[o] <- sum( rowSums(distances * nets_partitions_obs[[o]]) > 0 ) 
  }
  
  allsameties <- matrix(0,length(nets_partitions_sample),n_obs)
  allavgroupsameties <- matrix(0,length(nets_partitions_sample),n_obs)
  alldensities  <- matrix(0,length(nets_partitions_sample),n_obs)
  allindsameties <- matrix(0,length(nets_partitions_sample),n_obs)
  for(s in 1:length(nets_partitions_sample)){
    for(o in 1:n_obs){
      allsameties[s,o] <- sum(distances * nets_partitions_sample[[s]][[o]]) / 2
      allavgroupsameties[s,o] <- sum(distances * nets_partitions_sample[[s]][[o]]) / (2*max(partitions_sample[[s]][,o], na.rm = T))
      alldensities[s,o] <- computedensity(distances, partitions_sample[[s]][,o])$average
      allindsameties[s,o] <- sum( rowSums(distances * nets_partitions_sample[[s]][[o]]) > 0 )
    }
  }
  
  minsame <- apply(allsameties,2, FUN = function(x){quantile(x, probs = 0.05)})
  maxsame <- apply(allsameties,2, FUN = function(x){quantile(x, probs = 0.95)})
  minavgroupsame <- apply(allavgroupsameties,2, FUN = function(x){quantile(x, probs = 0.05)})
  maxavgroupsame <- apply(allavgroupsameties,2, FUN = function(x){quantile(x, probs = 0.95)})
  mindensities <- apply(alldensities,2, FUN = function(x){quantile(x, probs = 0.05)})
  maxdensities <- apply(alldensities,2, FUN = function(x){quantile(x, probs = 0.95)})
  minindsame <- apply(allindsameties,2, FUN = function(x){quantile(x, probs = 0.05)})
  maxindsame <- apply(allindsameties,2, FUN = function(x){quantile(x, probs = 0.95)})
  
  p_same <- rep(0,n_obs)
  p_avgroupsame <- rep(0,n_obs)
  p_densities <- rep(0,n_obs)
  p_indsame <- rep(0,n_obs)
  
  for(o in 1:n_obs){
    p1 <- sum(allsameties[,o] <= sameties_obs[o]) / length(nets_partitions_sample)
    p2 <- sum(allsameties[,o] >= sameties_obs[o]) / length(nets_partitions_sample)
    if(p1 == 0) {
      p_same[o] <- p2
    } else if(p2 == 0) {
      p_same[o] <- p1
    } else {
      p_same[o] <- 2*min(p1,p2)
    }
    
    p1 <- sum(allavgroupsameties[,o] <= avgroupsameties_obs[o]) / length(nets_partitions_sample)
    p2 <- sum(allavgroupsameties[,o] >= avgroupsameties_obs[o]) / length(nets_partitions_sample)
    if(p1 == 0) {
      p_avgroupsame[o] <- p2
    } else if(p2 == 0) {
      p_avgroupsame[o] <- p1
    } else {
      p_avgroupsame[o] <- 2*min(p1,p2)
    }
    
    p1 <- sum(alldensities[,o] <= densities_obs[o]) / length(nets_partitions_sample)
    p2 <- sum(alldensities[,o] >= densities_obs[o]) / length(nets_partitions_sample)
    if(p1 == 0) {
      p_densities[o] <- p2
    } else if(p2 == 0) {
      p_densities[o] <- p1
    } else {
      p_densities[o] <- 2*min(p1,p2)
    }
    
    p1 <- sum(allindsameties[,o] <= indsameties_obs[o]) / length(nets_partitions_sample)
    p2 <- sum(allindsameties[,o] >= indsameties_obs[o]) / length(nets_partitions_sample)
    if(p1 == 0) {
      p_indsame[o] <- p2
    } else if(p2 == 0) {
      p_indsame[o] <- p1
    } else {
      p_indsame[o] <- 2*min(p1,p2)
    }
    
  }
  
  
  res <- list(number_same_ties = list(observed = sameties_obs,
                                      mean_sample = apply(allsameties,2,mean),
                                      sd_sample = apply(allsameties,2,sd),
                                      min_CI_95 = minsame,
                                      max_CI_95 = maxsame,
                                      p = p_same),
              averagepergroup_number_same_ties = list(observed = avgroupsameties_obs,
                                                      mean_sample = apply(allavgroupsameties,2,mean),
                                                      sd_sample = apply(allavgroupsameties,2,sd),
                                                      min_CI_95 = minavgroupsame,
                                                      max_CI_95 = maxavgroupsame,
                                                      p = p_avgroupsame),
              intratie_densities = list(observed = densities_obs,
                                        mean_sample = apply(alldensities,2,mean),
                                        sd_sample = apply(alldensities,2,sd),
                                        min_CI_95 = mindensities,
                                        max_CI_95 = maxdensities,
                                        p = p_densities),
              individual_number_same_ties = list(observed = indsameties_obs,
                                                 mean_sample = apply(allindsameties,2,mean),
                                                 sd_sample = apply(allindsameties,2,sd),
                                                 min_CI_95 = minindsame,
                                                 max_CI_95 = maxindsame,
                                                 p = p_indsame))
  
  return(res)
}
