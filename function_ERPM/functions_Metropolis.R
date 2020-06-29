######################################################################
## Simulation and estimation of Exponential Random Partition Models ## 
## Function implementing the Metropolis HAstings algorithm to       ##
## through partitions given a certain model specification           ##
## Author: Marion Hoffman                                           ##
######################################################################



draw_Metropolis <- function(theta, 
                            first.partition, 
                            nodes, 
                            effects, 
                            objects, 
                            burnin, 
                            thining, 
                            num.steps, 
                            mini.steps, 
                            neighborhood, 
                            sizes.allowed,
                            sizes.simulated,
                            return.all.partitions = F) {
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)
  
  # instantiate with the starting network
  current.partition <- first.partition
  current.z <- computeStatistics(current.partition, nodes, effects, objects)
  current.logit <- theta * current.z
  
  # store the statistics collected for all networks simulated after the burn in
  all.z <-c()
  
  # store the partitions if needed (for GOF)
  if(return.all.partitions){
    all.partitions <- c()
  }
    
  end.walk <- FALSE
  cpt <- 0
  cpt2 <- 0
  cpt_thining <- 0
  
  while(!end.walk){
      
    #print("partition")
    #print(current.partition)
    
    ## IF NEIGHBORHOOD IS: MERGES OR SPLITS OF 2 (Pi_1)
    if(neighborhood == 1 && is.null(sizes.allowed)) {
      new.partition <- sample_new_partition_p1(current.partition, mini.steps)
    } else if(neighborhood == 2 && is.null(sizes.allowed)){
      new.partition <- sample_new_partition_p2(current.partition, mini.steps)
    }
    if(neighborhood == 1 && !is.null(sizes.allowed)) {
      new.partition <- sample_new_partition_p1_restricted(current.partition, mini.steps, sizes.simulated)
    } else if(neighborhood == 2 && !is.null(sizes.allowed)){
      new.partition <- sample_new_partition_p2_restricted(current.partition, mini.steps, sizes.simulated)
    }
       
    # compute new statistics only if it changed
    if(!all(current.partition == new.partition)) {
      new.z <- computeStatistics(new.partition, nodes, effects, objects)
      #new.z <- computeChangeStatistics(current.z, current.partition, new.partition, old_g1, old_g2, new_g1, new_g2, nodes, effects, objects)
    } else {
      new.z <- current.z
    }
    new.logit <- theta * new.z
      
    #print("stats")
    #print(new.z)
    
    # chose whether to change or not
    if(mini.steps == "normalized") {
      
      if(neighborhood == 1 && is.null(sizes.allowed)) {
        current.size <- compute_size_neighborhood_p1(current.partition)
        new.size <- compute_size_neighborhood_p1(new.partition)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else if(neighborhood == 2 && is.null(sizes.allowed)) {
        current.size <- compute_size_neighborhood_p2(current.partition)
        new.size <- compute_size_neighborhood_p2(new.partition)
        neighborhoods.ratio <- (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
      }
      if(neighborhood == 1 && !is.null(sizes.allowed)) {
        current.size <- compute_size_neighborhood_p1_restricted(current.partition, sizes.simulated)
        new.size <- compute_size_neighborhood_p1_restricted(new.partition, sizes.simulated)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else if(neighborhood == 2 && !is.null(sizes.allowed)) {
        current.size <- compute_size_neighborhood_p2_restricted(current.partition, sizes.simulated)
        new.size <- compute_size_neighborhood_p2_restricted(new.partition, sizes.simulated)
        neighborhoods.ratio <- (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
      }
      hastings.ratio <- (exp(sum(new.logit) - sum(current.logit))) * neighborhoods.ratio
   
    } else if(mini.steps == "selfloops"){
      
      hastings.ratio <- exp(sum(new.logit) - sum(current.logit))
      
    }
    proba.change <- min(1,hastings.ratio)
    change.made <- (runif(1,0,1) <= proba.change)

    #print("hastings proba.change and change.made")
    #print(hastings.ratio)
    #print(proba.change)
    #print(change.made)
    
    old.partition <- current.partition
    if(change.made) {
      current.partition <- new.partition
      current.z <- new.z
      current.logit <- new.logit
    }
    
    cpt <- cpt + 1
    if(cpt > burnin) cpt_thining <- cpt_thining + 1
    
    
    # storing results if all sizes are allowed
    if(is.null(sizes.allowed)){
      # store the results if we are out of burnin
      if(cpt >= burnin && cpt_thining == thining) {
        all.z <- rbind(all.z,current.z)
        cpt_thining <- 0
        
        # store all partitions if needed
        if(return.all.partitions){
          all.partitions <- rbind(all.partitions,new.partition)
        }
        
      }
      
      # stop the walk if number of steps reached 
      end.walk <- (cpt >= (burnin+thining*num.steps))
    }

  # storing results if sizes are constrained
  if(!is.null(sizes.allowed)){

    #print(check_sizes(current.partition,sizes.allowed))
    #print(current.partition)
    
    # if the partition is sampled with the right sizes
    # we keep the swapping strategy
    #if(check_sizes(current.partition,sizes.allowed)){
      #neighborhood <- 2
      
      # if wrong sizes, we continue searching, and if we are out of tolerated simulated partitions, we go back
      # when out of allowed partitions we keep the divide/merge strategy
    #} else {
      #neighborhood <- 2
    #  if(!check_sizes(current.partition,sizes.simulated)){
    #    current.partition <- old.partition
    #    if(check_sizes(current.partition,sizes.allowed)) neighborhood <- 1
    #  } 
    #}
  
    # store the results if we are out of burnin
    if(cpt >= burnin && cpt_thining == thining) {
      cpt_thining <- thining - 1
      if(check_sizes(current.partition,sizes.allowed)){
        cpt_thining <- 0
        all.z <- rbind(all.z,current.z)
        #if(nrow(all.z)>1) print(all.z[nrow(all.z),] - all.z[nrow(all.z)-1,])
        cpt2 <- cpt2 + 1
        
        # store all partitions if needed
        if(return.all.partitions){
          all.partitions <- rbind(all.partitions,new.partition)
        }
      } 
    }
    
    # stop the walk if number of steps reached 
    end.walk <- (cpt2 >= num.steps)
  }

  }
  
  # TODO: change this hack
  row.names(all.z) <- rep("", dim(all.z)[1])
  
  # compute the average statistics and the final network generated
  if(!return.all.partitions) {
    return( list("draws" = all.z, 
                 "last.partition" = current.partition))
  } else {
    return( list("draws" = all.z, 
               "last.partition" = current.partition,
               "all.partitions" = all.partitions)) 
  }
  
}


## NEIGHBORHOOD PI 1: only merges and divisions of 2 groups
compute_size_neighborhood_p1 <- function(partition){
  
  # find isolates, pairs and groups>2
  isolates <- c()
  pairs <- c()
  others <- c()
  num.groups <- max(partition)
  num.nodes <- length(partition)
  
  for(g in 1:num.groups){
    if(length(which(partition == g)) == 1) {
      isolates <- c(isolates,g)
    } else if(length(which(partition == g)) == 2) {
      pairs <- c(pairs,g)
    } else {
      others <- c(others,g)
    }
  }
  
  nums.swaps <- rep(0,num.nodes)
  done.pairs <- rep(0,num.groups)
  for(i in 1:num.nodes){
    g <- partition[i]
    if(g %in% isolates){
      nums.swaps[i] <- num.groups - sum(isolates<=i) #can join any other group, except other isolates before (because already done)
    }
    if(g %in% pairs){
      nums.swaps[i] <- num.groups -1 # can join any other group, except if the possibilty of splitting the pair is already counted
      if(done.pairs[g] == 0){
        done.pairs[g] <- 1
        nums.swaps[i] <- nums.swaps[i] + 1
      }
    }
    if(g %in% others){
      nums.swaps[i] <- num.groups # can join any other group or create its own isolate
    }
  }

  num.swaps <- sum(nums.swaps)
  
  return(list(isolates = isolates,
              pairs = pairs,
              others = others,
              nums.swaps = nums.swaps,
              num.swaps = num.swaps))
  
}

sample_new_partition_p1 <- function(current.partition, mini.steps){
  
  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- compute_size_neighborhood_p1(current.partition)
  isolates <- size$isolates
  pairs <- size$pairs
  others <- size$others
  nums.swaps <- size$nums.swaps
  num.swaps <- size$num.swaps
  
  if(mini.steps == "normalized") {
    total <- num.swaps
  } 

  pick.1 <- sample(total,1)
  #print(pick.1)
    
  # decide which actor to swap
  all.groups <- 1:num.groups
  new.partition <- current.partition
  
  for(i in 1:num.nodes){
    
    if(i == 1) {start <- 0} else {start <- sum(nums.swaps[1:i-1])}
    if(i == num.nodes) {end <- num.swaps} else {end <- sum(nums.swaps[1:i])}
    
    if(pick.1 > start && pick.1 <= end) {
      g <- current.partition[i]
      
      # an isolate is randomly joining another group except isolates before i
      if(g %in% isolates){
        previousisolates <- isolates[isolates<=g]
        newg <- sample(all.groups[!all.groups %in% previousisolates],1)
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }
      
      # a member of a pair is randomly joining another group 
      # or creating its isolate if it's the first of the pair
      if(g %in% pairs){
        
        members <- which(current.partition == g)
        isfirst <- members[1] == i
        
        if(isfirst){
          newg <- sample(all.groups,1)
          if(newg == g){
            newg <- num.groups + 1
          }
        } else {
          newg <- sample(all.groups[all.groups != g],1)
        }
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }
      
      # a member of a bigger group is randomly joining another group
      # or creating its isolate
      if(g %in% others){
        newg <- sample(all.groups,1)
        if(newg == g){
          newg <- num.groups + 1
        }
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }
    }
  }
  
  return(new.partition)
  
}


## NEIGHBORHOOD PI 1: only merges and divisions of 2 groups
compute_size_neighborhood_p2 <- function(partition){
  
  # calculate the number of neighbor partitions
  num.nodes <- length(partition)
  num.groups <- max(partition)
  num.merges <- num.groups*(num.groups-1)/2
  nums.divisions <- rep(0,num.groups)
  num.divisions <- 0
  
  for(k in 1:num.groups){
    sk <- length(which(partition == k))
    nums.divisions[k] <- 2^(sk-1) - 1
    num.divisions <- num.divisions + nums.divisions[k]
  }
  
  return(list(num.merges = num.merges,
              num.divisions = num.divisions,
              nums.divisions = nums.divisions))
}


sample_new_partition_p2 <- function(current.partition, mini.steps){
  
  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- compute_size_neighborhood_p2(current.partition)
  num.merges <- size$num.merges
  num.divisions <- size$num.divisions
  nums.divisions <- size$nums.divisions

  if(mini.steps == "normalized") {
    total <- num.merges + num.divisions 
  } else if(mini.steps == "selfloops") {
    total <- 2^(num.nodes-1)-1
  }
  
  # decide between merge or division (or self loop)
  pick.1 <- sample(total,1)
  all.groups <- 1:num.groups
  new.partition <- current.partition
  
  # if merge
  if(pick.1 <= num.merges) {
    
    # pick 2 groups to merge
    old_g1 <- sample(num.groups,1)
    other.groups <- all.groups[all.groups != old_g1]
    old_g2 <- sample(other.groups,1)
    gmin <- min(old_g1,old_g2)
    gmax <- max(old_g1,old_g2)
    old_g1 <- gmin
    old_g2 <- gmax
    new_g1 <- old_g1
    new_g2 <- 0
    
    # reassign one of the groups and remove useless ids
    new.partition[which(new.partition == old_g2)] <- old_g1
    new.partition <- order_groupids(new.partition)
    
  }
  # if division
  else if(pick.1 <= (num.merges+num.divisions)){
    
    # pick group to divide
    pick.2 <- sample(num.divisions,1)
    old_g1 <- 1
    cpt2 <- 1
    g <- 1
    found <- FALSE
    while(g<=num.groups && !found){
      if(pick.2 >= cpt2 && pick.2 <= sum(nums.divisions[1:g])){
        old_g1 <- g
        found <- TRUE
      }
      cpt2 <- cpt2 + nums.divisions[g]
      g <- g+1
    }
    old_g2 <- 0
    new_g1 <- old_g1
    new_g2 <- num.groups + 1
    
    # reassign one of the groups and remove useless ids
    new.groups <- sample(c(old_g1,num.groups+1),length(which(new.partition==old_g1)),replace=T)
    found <- (length(unique(new.groups)) > 1)
    while(!found){
      new.groups <- sample(c(old_g1,num.groups+1),length(which(new.partition==old_g1)),replace=T)
      found <- (length(unique(new.groups)) > 1)
    }
    new.partition[which(new.partition == old_g1)] <- new.groups
  }
  
  return(new.partition)
  
}



## NEIGHBORHOOD PI 1: only merges and divisions of 2 groups
# BUT REMOVE OUT OF SAMPLE PARTITIONS
# WARNING!!!!: for now it only works if sizes allowed are an interval ([size_min,size_max])
compute_size_neighborhood_p1_restricted <- function(partition, sizes.simulated){
  
  # check if current partition is allowed
  if(!check_sizes(partition, sizes.simulated)) stop("The partition we are in is not allowed.")
  
  # find maximum size allowed
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)
  
  # find isolates, pairs and groups>2
  isolates <- c()
  pairs <- c()
  others_allowed <- c()
  others_notallowed <- c()
  num.groups <- max(partition)
  num.nodes <- length(partition)
  
  for(g in 1:num.groups){
    if(length(which(partition == g)) == 1) {
      isolates <- c(isolates,g)
    } else if(length(which(partition == g)) == 2) {
      pairs <- c(pairs,g)
    } else if(length(which(partition == g)) < smax){
      others_allowed <- c(others_allowed,g)
    } else {
      others_notallowed <- c(others_notallowed,g)
    }
  }
  
  nums.swaps <- rep(0,num.nodes)
  done.pairs <- rep(0,num.groups)
  for(i in 1:num.nodes){
    g <- partition[i]
    if(g %in% isolates){
      # can join any other group (<smax), except other isolates before (because already done)
      nums.swaps[i] <- length(others_allowed) + length(pairs) + sum(isolates>i) 
    }
    if(g %in% pairs && smin == 1){
      # if isolates are allowed, then it can join any other group (smax), and splitting if it's not already counted
      nums.swaps[i] <- length(others_allowed) + length(pairs) + length(isolates) - 1 
      if(done.pairs[g] == 0){
        done.pairs[g] <- 1
        nums.swaps[i] <- nums.swaps[i] + 1
      }
    }
    if(g %in% others_allowed){
      # can join any other group (<smax) or create its own isolate (if it's allowed)
      nums.swaps[i] <- length(others_allowed) + length(pairs) + length(isolates) - 1
      if(smin == 1) nums.swaps[i] <- nums.swaps[i] + 1
    }
    if(g %in% others_notallowed){
      # can join any other group (<smax) or create its own isolate (if it's allowed)
      nums.swaps[i] <- length(others_allowed) + length(pairs) + length(isolates)
      if(smin == 1) nums.swaps[i] <- nums.swaps[i] + 1
    }
  }
  
  num.swaps <- sum(nums.swaps)
  
  return(list(isolates = isolates,
              pairs = pairs,
              others_allowed = others_allowed,
              others_notallowed = others_notallowed,
              nums.swaps = nums.swaps,
              num.swaps = num.swaps))
  
}

sample_new_partition_p1_restricted <- function(current.partition, mini.steps, sizes.simulated){
  
  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- compute_size_neighborhood_p1_restricted(current.partition, sizes.simulated)
  isolates <- size$isolates
  pairs <- size$pairs
  others_allowed <- size$others_allowed
  others_notallowed <- size$others_notallowed
  nums.swaps <- size$nums.swaps
  num.swaps <- size$num.swaps
  
  # allowed sizes
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)
  
  if(mini.steps == "normalized") {
    total <- num.swaps
  } 
  
  pick.1 <- sample(total,1)
  #print(pick.1)
  
  # decide which actor to swap
  all.groups <- 1:num.groups
  all.groups_allowed <- all.groups[!all.groups %in% others_notallowed]
  new.partition <- current.partition
  
  for(i in 1:num.nodes){
    
    if(i == 1) {start <- 0} else {start <- sum(nums.swaps[1:i-1])}
    if(i == num.nodes) {end <- num.swaps} else {end <- sum(nums.swaps[1:i])}
    
    if(pick.1 > start && pick.1 <= end) {
      g <- current.partition[i]
      
      # an isolate is randomly joining another group except isolates before i
      if(g %in% isolates){
        previousisolates <- isolates[isolates<=g]
        newg <- sample(all.groups_allowed[!all.groups_allowed %in% previousisolates],1)
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }
      
      # a member of a pair is randomly joining another group 
      # or creating its isolate if it's the first of the pair
      if(g %in% pairs){
        
        members <- which(current.partition == g)
        isfirst <- members[1] == i
        
        if(isfirst){
          newg <- sample(all.groups_allowed,1)
          if(newg == g){
            newg <- num.groups + 1
          }
        } else {
          newg <- sample(all.groups_allowed[all.groups_allowed != g],1)
        }
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }
      
      # a member of a bigger group is randomly joining another group
      # or creating its isolate
      if(g %in% others_allowed){
        if(smin == 1){
          newg <- sample(all.groups_allowed,1)
          if(newg == g) newg <- num.groups + 1
        } else {
          newg <- sample(all.groups_allowed[all.groups_allowed != g],1)
        }
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }
    }
  }
  
  return(new.partition)
  
}


## NEIGHBORHOOD PI 1: only merges and divisions of 2 groups
# BUT REMOVE OUT OF SAMPLE PARTITIONS
# WARNING!!!!: for now it only works if sizes allowed are an interval ([size_min,size_max])
compute_size_neighborhood_p2_restricted <- function(partition, sizes.simulated){
  
  # check if current partition is allowed
  if(!check_sizes(partition, sizes.simulated)) stop("The partition we are in is not allowed.")
  
  # find maximum size allowed
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)
  
  # calculate the number of neighbor partitions
  num.nodes <- length(partition)
  num.groups <- max(partition)
  num.merges <- 0
  merges <- matrix(0,num.groups,num.groups)
  nums.divisions <- rep(0,num.groups)
  num.divisions <- 0
  
  # merges
  if(num.groups>1){
    for(g1 in 1:(num.groups-1)){
      for(g2 in (g1+1):num.groups){
        if(length(which(partition == g1)) + length(which(partition == g2)) <= smax) {
          num.merges <- num.merges + 1
          merges[g1,g2] <- 1
        }
      }
    }
  }
  
  
  # divisions
  for(k in 1:num.groups){
    sg <- length(which(partition == k))
    if(smin>1){
      extras <- sum(unlist(lapply(1:(smin-1),FUN=function(x){choose(sg,x)})))
      if(sg%%2 == 0 && sg/2 < smin) extras <- extras - choose(sg,sg/2)/2
    } else {
      extras <- 0
    }
    nums.divisions[k] <- 2^(sg-1) - 1 - extras
    num.divisions <- num.divisions + nums.divisions[k]
  }
  
  return(list(num.merges = num.merges,
              merges = merges,
              num.divisions = num.divisions,
              nums.divisions = nums.divisions))
}


sample_new_partition_p2_restricted <- function(current.partition, mini.steps, sizes.simulated){
  
  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- compute_size_neighborhood_p2_restricted(current.partition, sizes.simulated)
  
  num.merges <- size$num.merges
  merges <- size$merges
  
  num.divisions <- size$num.divisions
  nums.divisions <- size$nums.divisions
  
  # sizes allowed
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)
  
  if(mini.steps == "normalized") {
    total <- num.merges + num.divisions 
  } else if(mini.steps == "selfloops") {
    total <- 2^(num.nodes-1)-1
  }
  
  # decide between merge or division (or self loop)
  pick.1 <- sample(total,1)
  all.groups <- 1:num.groups
  new.partition <- current.partition
  
  # if merge
  if(pick.1 <= num.merges) {
    
    # pick 2 groups to merge
    allindexes <- which(merges == 1)
    p <- allindexes[pick.1]
    indexes <- arrayInd(p, dim(merges))
    old_g1 <- indexes[,1]
    old_g2 <- indexes[,2]
    
    
    if(merges[old_g1,old_g2] == 0) {
      print("oproblem")
    }
    
    # reassign one of the groups and remove useless ids
    new.partition[which(new.partition == old_g2)] <- old_g1
    new.partition <- order_groupids(new.partition)
    
  }
  # if division
  else if(pick.1 <= (num.merges+num.divisions)){
    
    # pick group to divide
    pick.2 <- sample(num.divisions,1)
    old_g1 <- 1
    cpt2 <- 1
    g <- 1
    found <- FALSE
    while(g<=num.groups && !found){
      if(pick.2 >= cpt2 && pick.2 <= sum(nums.divisions[1:g])){
        old_g1 <- g
        found <- TRUE
      }
      cpt2 <- cpt2 + nums.divisions[g]
      g <- g+1
    }
    old_g2 <- 0
    new_g1 <- old_g1
    new_g2 <- num.groups + 1
    
    # reassign one of the groups and remove useless ids
    new.groups <- sample(c(old_g1,num.groups+1),length(which(new.partition==old_g1)),replace=T)
    found <- (length(unique(new.groups)) > 1) && min(table(new.groups)) >= smin
    while(!found){
      new.groups <- sample(c(old_g1,num.groups+1),length(which(new.partition==old_g1)),replace=T)
      found <- (length(unique(new.groups)) > 1) && min(table(new.groups)) >= smin
    }
    new.partition[which(new.partition == old_g1)] <- new.groups
  }
  
  return(new.partition)
  
}

