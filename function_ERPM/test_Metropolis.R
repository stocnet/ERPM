######################################################################
## Simulation and estimation of Exponential Random Partition Models ## 
## Function implementing the Metropolis HAstings algorithm to       ##
## through partitions given a certain model specification           ##
## Author: Marion Hoffman                                           ##
######################################################################


# SINGLE PARTITION PROCEDURE
draw_Metropolis_single <- function(theta, # model parameters
                                   first.partition, # starting partition for the Markov chain
                                   nodes, # nodeset (data frame)
                                   effects, # effects/sufficient statistics (list with a vector "names", and a vector "objects")
                                   objects, # objects used for statistics calculation (list with a vector "name", and a vector "object")
                                   burnin, # integer for the number of burn-in steps before sampling
                                   thining, # integer for the number of thining steps between sampling
                                   num.steps, # number of samples
                                   mini.steps, # type of transition, either "normalized", either "self-loops" (take "normalized")
                                   neighborhood, # way of choosing partitions: probability vector (proba actors swap, proba merge/division, proba single actor move)
                                   sizes.allowed, # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                                   sizes.simulated,# vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
                                   return.all.partitions = F) # option to return the sampled partitions on top of their statistics (for GOF)
{
  
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
    
    new.step <- draw_step_single(theta,
                                 current.partition, 
                                 current.logit,
                                 current.z,
                                 nodes, 
                                 effects, 
                                 objects, 
                                 mini.steps, 
                                 neighborhood, 
                                 sizes.allowed, 
                                 sizes.simulated)
    
    proba.change <- min(1,new.step$hastings.ratio)
    change.made <- (runif(1,0,1) <= proba.change)
    
    # update partition if needed
    old.partition <- current.partition
    if(change.made) {
      current.partition <- new.step$new.partition
      current.z <- new.step$new.z
      current.logit <- new.step$new.logit
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
          all.partitions <- rbind(all.partitions,current.partition)
        }
        
      }
      
      # stop the walk if number of steps reached 
      end.walk <- (cpt >= (burnin+thining*num.steps))
    }
    
    # storing results if sizes are constrained
    if(!is.null(sizes.allowed)){
      
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
            all.partitions <- rbind(all.partitions,current.partition)
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
                 "last.partition" = current.partition,
                 "all.partitions" = NULL))
  } else {
    return( list("draws" = all.z, 
                 "last.partition" = current.partition,
                 "all.partitions" = all.partitions)) 
  }
  
}


# MULTIPLE PARTITIONS PROCEDURE
draw_Metropolis_multiple <- function(theta, # model parameters
                                     first.partitions, # starting partition for the Markov chain
                                     presence.tables,  # matrix indicating which actors were present for each observations (mandatory)
                                     nodes, # nodeset (data frame)
                                     effects, # effects/sufficient statistics (list with a vector "names", and a vector "objects")
                                     objects, # objects used for statistics calculation (list with a vector "name", and a vector "object")
                                     burnin, # integer for the number of burn-in steps before sampling
                                     thining, # integer for the number of thining steps between sampling
                                     num.steps, # number of samples
                                     mini.steps, # type of transition, either "normalized", either "self-loops" (take "normalized")
                                     neighborhood, # way of choosing partitions: probability vector (proba actors swap, proba merge/division, proba single actor move)
                                     sizes.allowed, # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                                     sizes.simulated,# vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
                                     return.all.partitions = F) # option to return the sampled partitions on top of their statistics (for GOF)
{
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)
  num.obs <- ncol(presence.tables)
  
  # instantiate with the starting network
  current.partitions <- first.partitions
  current.z.contributions <- computeStatistics_multiple(current.partitions, presence.tables, nodes, effects, objects)
  current.z <- rowSums(current.z.contributions)
  current.logit <- theta * current.z
  
  # store the statistics collected for all networks simulated after the burn in
  all.z <-c()
  
  # store the partitions if needed (for GOF)
  if(return.all.partitions){
    all.partitions <- list()
    cpt_all <- 0
  }
  
  end.walk <- FALSE
  cpt <- 0
  cpt2 <- 0
  cpt_thining <- 0
  
  accepted <- rep(0,10)
  tried <- rep(0,10)
  
  while(!end.walk){
    
    start <- Sys.time()
    new.step <- draw_step_multiple(theta,
                                   current.partitions, 
                                   current.logit,
                                   current.z.contributions,
                                   current.z,
                                   presence.tables,  
                                   nodes, 
                                   effects, 
                                   objects, 
                                   mini.steps, 
                                   neighborhood, 
                                   sizes.allowed, 
                                   sizes.simulated)
    end <- Sys.time()
    
    proba.change <- min(1,new.step$hastings.ratio)
    change.made <- (runif(1,0,1) <= proba.change)
    
    if(change.made) print("change made")
    
    move <- new.step$move
    if(move == "swap") {
      tried[1] <- tried[1] + 1
      if(change.made) accepted[1] <- accepted[1] + 1
    } else if(move == "mergediv"){
      tried[2] <- tried[2] + 1
      if(change.made) accepted[2] <- accepted[2] + 1
    } else if(move == "single"){
      tried[3] <- tried[3] + 1
      if(change.made) accepted[3] <- accepted[3] + 1
    } else if(move == "double"){
      tried[4] <- tried[4] + 1
      if(change.made) accepted[4] <- accepted[4] + 1
    } else if(move == "swap2"){
      tried[5] <- tried[5] + 1
      if(change.made) accepted[5] <- accepted[5] + 1
    } else if(move == "exchange"){
      tried[6] <- tried[6] + 1
      if(change.made) accepted[6] <- accepted[6] + 1
    } else if(move == "partition_exchange"){
      tried[7] <- tried[7] + 1
      if(change.made) accepted[7] <- accepted[7] + 1
    } else if(move == "partition_next"){
      tried[8] <- tried[8] + 1
      if(change.made) accepted[8] <- accepted[8] + 1
    } else if(move == "group_exchange"){
      tried[9] <- tried[9] + 1
      if(change.made) accepted[9] <- accepted[9] + 1
    } else if(move == "group_next"){
      tried[10] <- tried[10] + 1
      if(change.made) accepted[10] <- accepted[10] + 1
    }
    
    #print(move)
    #print(end-start)
    
    # update partition if needed
    old.partitions <- current.partitions
    if(change.made) {
      current.partitions <- new.step$new.partitions
      current.z.contributions <- new.step$new.z.contributions
      current.z <- new.step$new.z
      current.logit <- new.step$new.logit
    }
    
    cpt <- cpt + 1
    if(cpt > burnin) cpt_thining <- cpt_thining + 1
    
    print(cpt)
    
    # storing results if all sizes are allowed
    if(is.null(sizes.allowed)){
      
      # store the results if we are out of burnin
      if(cpt >= burnin && cpt_thining == thining) {
        all.z <- rbind(all.z,current.z)
        cpt_thining <- 0
        
        # store all partitions if needed
        if(return.all.partitions){
          cpt_all <- cpt_all + 1
          all.partitions[[cpt_all]] <- current.partitions
        }
        
      }
      
      # stop the walk if number of steps reached 
      end.walk <- (cpt >= (burnin+thining*num.steps))
    }
    
    # storing results if sizes are constrained
    if(!is.null(sizes.allowed)){
      
      # store the results if we are out of burnin
      if(cpt >= burnin && cpt_thining == thining) {
        
        cpt_thining <- thining - 1
        
        # check all partitions one by one
        check_all <- T
        for(o in 1:num.obs) check_all <- check_all && check_sizes(current.partitions[as.logical(presence.tables[,o]),o],sizes.allowed)
        
        if(check_all){
          
          cpt_thining <- 0
          all.z <- rbind(all.z,current.z)
          cpt2 <- cpt2 + 1
          
          print(nrow(all.z))
          
          # store all partitions if needed
          if(return.all.partitions){
            cpt_all <- cpt_all + 1
            all.partitions[[cpt_all]] <- current.partitions
          }
        } 
      }
      
      # stop the walk if number of steps reached 
      end.walk <- (cpt2 >= num.steps)
    }
    
  }
  
  # TODO: change this hack
  row.names(all.z) <- rep("", dim(all.z)[1])
  
  print("acceptance")
  print(accepted / tried)
  print("correlations")
  autocors <- rep(0,num.effects)
  for(e in 1:num.effects) autocors[e] <- cor(all.z[1:(num.steps-1),e],all.z[2:num.steps,e])
  print(autocors)
  
  # compute the average statistics and the final network generated
  if(!return.all.partitions) {
    return( list("draws" = all.z, 
                 "last.partitions" = current.partitions,
                 "all.partitions" = NULL))
  } else {
    return( list("draws" = all.z, 
                 "last.partitions" = current.partitions,
                 "all.partitions" = all.partitions)) 
  }
  
}


# function to draw next partition and calculate HAstings ratio (one step in the Metropolis algorithm)
draw_step_single <- function(theta,
                             current.partition, 
                             current.logit, 
                             current.z,
                             nodes, 
                             effects, 
                             objects, 
                             mini.steps, 
                             neighborhood, 
                             sizes.allowed, 
                             sizes.simulated){
  
  move <- sample(c("swap","mergediv","single","double","swap2","exchange"),1,prob=neighborhood)
  
  ## IF MODE IS "swap" = swap two actors, "mergediv" = merge 2 groups or split a group in two, "single" = move one actor
  if(move == "swap" && is.null(sizes.allowed)) {
    current.size <- compute_size_neighborhood_p1(current.partition)
    if(current.size$num.swaps > 0) {
      new.partition <- sample_new_partition_p1(current.partition, mini.steps, current.size)
    } else {
      new.partition <- current.partition
    }
  } else if(move == "mergediv" && is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p2(current.partition)
    new.partition <- sample_new_partition_p2(current.partition, mini.steps, current.size)
  } else if(move == "single" && is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p3(current.partition)
    new.partition <- sample_new_partition_p3(current.partition, mini.steps, current.size)
  } else if(move == "double" && is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p4(current.partition)
    new.partition <- sample_new_partition_p4(current.partition, mini.steps, current.size)
  } else if(move == "swap2" && is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p5(current.partition)
    new.partition <- sample_new_partition_p5(current.partition, mini.steps, current.size)
  } else if(move == "exchange" && is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p6(current.partition)
    new.partition <- sample_new_partition_p6(current.partition, mini.steps, current.size)
  }
  # Alternative if sizes are restricted
  if(move == "swap" && !is.null(sizes.allowed)) {
    current.size <- compute_size_neighborhood_p1_restricted(current.partition, sizes.simulated)
    if(current.size$num.swaps > 0) {
      new.partition <- sample_new_partition_p1_restricted(current.partition, mini.steps, sizes.simulated, current.size)
    } else {
      new.partition <- current.partition
    }
  } else if(move == "mergediv" && !is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p2_restricted(current.partition, sizes.simulated)
    new.partition <- sample_new_partition_p2_restricted(current.partition, mini.steps, sizes.simulated, current.size)
  } else if(move == "single" && !is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p3_restricted(current.partition, sizes.simulated)
    new.partition <- sample_new_partition_p3_restricted(current.partition, mini.steps, sizes.simulated, current.size)
  } else if(move == "double" && !is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p4_restricted(current.partition, sizes.simulated)
    new.partition <- sample_new_partition_p4_restricted(current.partition, mini.steps, sizes.simulated, current.size)
  } else if(move == "swap2" && !is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p5_restricted(current.partition, sizes.simulated)
    new.partition <- sample_new_partition_p5_restricted(current.partition, mini.steps, sizes.simulated, current.size)
  } else if(move == "exchange" && !is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p6_restricted(current.partition, sizes.simulated)
    new.partition <- sample_new_partition_p6_restricted(current.partition, mini.steps, sizes.simulated, current.size)
  }
  
  if(!check_sizes(new.partition, sizes.simulated)) {
    print("old partition")
    print(current.partition)
    print("new partition")
    print(new.partition)
    print("neighborhood")
    print(move)
    stop("The partition we are in is not allowed.")
  }
  
  # compute new statistics only if it changed
  if(!all(current.partition == new.partition, na.rm=T)) {
    new.z <- computeStatistics(new.partition, nodes, effects, objects)
  } else {
    new.z <- current.z
  }
  new.logit <- theta * new.z
  
  # chose whether to change or not
  if(mini.steps == "normalized") {
    
    if(move == "swap" && is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p1(new.partition)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratio <- 1
      }
    } else if(move == "mergediv" && is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p2(new.partition)
      neighborhoods.ratio <- (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
    } else if(move == "single" && is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p3(new.partition)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratio <- 1
      }
    } else if(move == "double" && is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p4(new.partition)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratio <- 1
      }
    } else if(move == "swap2" && is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p5(new.partition)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratio <- 1
      }
    } else if(move == "exchange" && is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p6(new.partition)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratio <- 1
      }
    }
    
    # alternative if sizes are restricted
    if(move == "swap" && !is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p1_restricted(new.partition, sizes.simulated)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratio <- 1
      }
    } else if(move == "mergediv" && !is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p2_restricted(new.partition, sizes.simulated)
      neighborhoods.ratio <- (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
    } else if(move == "single" && !is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p3_restricted(new.partition, sizes.simulated)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratio <- 1
      }
    } else if(move == "double" && !is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p4_restricted(new.partition, sizes.simulated)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratio <- 1
      }
    } else if(move == "swap2" && !is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p5_restricted(new.partition, sizes.simulated)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratio <- 1
      }
    } else if(move == "exchange" && !is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p6_restricted(new.partition, sizes.simulated)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratio <- 1
      }
    }
    
    hastings.ratio <- (exp(sum(new.logit) - sum(current.logit))) * neighborhoods.ratio
    
  } else if(mini.steps == "selfloops"){
    
    hastings.ratio <- exp(sum(new.logit) - sum(current.logit))
    
  }
  
  # check if current partition is allowed
  #if(!check_sizes(new.partition, sizes.simulated)) {
  #  print("old partition")
  #  print(current.partition)
  #  print("neighborhood size")
  #  print(current.size)
  #  print("new partition")
  #  print(new.partition)
  #  stop("The partition we are in is not allowed.")
  #} 
  
  return(list("hastings.ratio" = hastings.ratio,
              "new.partition" = new.partition,
              "new.z" = new.z,
              "new.logit" = new.logit))
}


# function to draw next partition and calculate HAstings ratio (one step in the Metropolis algorithm)
draw_step_multiple <- function(theta,
                               current.partitions, 
                               current.logit,
                               current.z.contributions,
                               current.z,
                               presence.tables,  
                               nodes, 
                               effects, 
                               objects, 
                               mini.steps, 
                               neighborhood, 
                               sizes.allowed, 
                               sizes.simulated){
  
  num.nodes <- nrow(nodes)
  num.obs <- ncol(presence.tables)
  new.partitions <- current.partitions
  
  rand.o <- sample(1:num.obs,1)
  nodes.rand.o <- as.logical(presence.tables[,rand.o])
  
  move <- sample(c("swap","mergediv","single","double","swap2","exchange","partition_exchange","partition_next","group_exchange","group_next","pair_exchange","pair_next"),1,prob=neighborhood)
  #move <- sample(c("swap","mergediv","single","double","swap2","exchange","partition_exchange","partition_next","group_exchange","group_next"),1,prob=neighborhood)
  
  #start <- Sys.time()
  
  ## IF MODE IS "swap" = swap two actors, "mergediv" = merge 2 groups or split a group in two, "single" = move one actor
  if(move == "swap" && is.null(sizes.allowed)) {
    current.size <- compute_size_neighborhood_p1(current.partitions[nodes.rand.o,rand.o])
    if(current.size$num.swaps > 0){
      sampled <- sample_new_partition_p1(current.partitions[nodes.rand.o,rand.o], mini.steps, current.size)
      new.p <- sampled$new.partition
      nodes.changed <- ((1:num.nodes)[nodes.rand.o])[sampled$nodes.changed]
    } else {
      new.p <- current.partitions[nodes.rand.o,rand.o]
      nodes.changed <- c()
    }
  } else if(move == "mergediv" && is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p2(current.partitions[nodes.rand.o,rand.o])
    sampled <- sample_new_partition_p2(current.partitions[nodes.rand.o,rand.o], mini.steps, current.size)
    new.p <- sampled$new.partition
    nodes.changed <- ((1:num.nodes)[nodes.rand.o])[sampled$nodes.changed]
  } else if(move == "single" && is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p3(current.partitions[nodes.rand.o,rand.o])
    sampled <- sample_new_partition_p3(current.partitions[nodes.rand.o,rand.o], mini.steps, current.size)
    new.p <- sampled$new.partition
    nodes.changed <- ((1:num.nodes)[nodes.rand.o])[sampled$nodes.changed]
  } else if(move == "double" && is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p4(current.partitions[nodes.rand.o,rand.o])
    new.p <- sample_new_partition_p4(current.partitions[nodes.rand.o,rand.o], mini.steps, current.size)
  } else if(move == "swap2" && is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p5(current.partitions[nodes.rand.o,rand.o])
    new.p <- sample_new_partition_p5(current.partitions[nodes.rand.o,rand.o], mini.steps, current.size)
  } else if(move == "exchange" && is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p6(current.partitions[nodes.rand.o,rand.o])
    sampled <- sample_new_partition_p6(current.partitions[nodes.rand.o,rand.o], mini.steps, current.size)
    new.p <- sampled$new.partition
    nodes.changed <- ((1:num.nodes)[nodes.rand.o])[sampled$nodes.changed]
  } else if(move == "partition_exchange" && is.null(sizes.allowed)){
    rand.o2 <- sample((1:num.obs)[-rand.o],1)
    nodes.rand.o2 <- as.logical(presence.tables[,rand.o2])
    current.size <- compute_size_neighborhood_p7(current.partitions[,rand.o],current.partitions[,rand.o2])
    new.ps <- sample_new_partition_p7(current.partitions[,rand.o], current.partitions[,rand.o2])
    new.p <- new.ps$new.partition1
    new.p2 <- new.ps$new.partition2
  } else if(move == "partition_next" && is.null(sizes.allowed)){
    if(rand.o == 1) rand.o2 <- 2
    if(rand.o == num.obs) rand.o2 <- num.obs-1
    if(rand.o > 1 && rand.o < num.obs) rand.o2 <- sample(c(rand.o-1,rand.o+1),1)
    nodes.rand.o2 <- as.logical(presence.tables[,rand.o2])
    current.size <- compute_size_neighborhood_p7(current.partitions[,rand.o],current.partitions[,rand.o2])
    new.ps <- sample_new_partition_p7(current.partitions[,rand.o], current.partitions[,rand.o2])
    new.p <- new.ps$new.partition1
    new.p2 <- new.ps$new.partition2
  } else if(move == "group_exchange" && is.null(sizes.allowed)){
    rand.o2 <- sample((1:num.obs)[-rand.o],1)
    nodes.rand.o2 <- as.logical(presence.tables[,rand.o2])
    current.size <- compute_size_neighborhood_p8(current.partitions[,rand.o],current.partitions[,rand.o2])
    new.p <- sample_new_partition_p8(current.partitions[,rand.o], current.partitions[,rand.o2], current.size)
  } else if(move == "group_next" && is.null(sizes.allowed)){
    if(rand.o == 1) rand.o2 <- 2
    if(rand.o == num.obs) rand.o2 <- num.obs-1
    if(rand.o > 1 && rand.o < num.obs) rand.o2 <- sample(c(rand.o-1,rand.o+1),1)
    nodes.rand.o2 <- as.logical(presence.tables[,rand.o2])
    current.size <- compute_size_neighborhood_p8(current.partitions[,rand.o],current.partitions[,rand.o2])
    new.p <- sample_new_partition_p8(current.partitions[,rand.o], current.partitions[,rand.o2], current.size)
  } else if(move == "pair_exchange" && is.null(sizes.allowed)){
    rand.o2 <- sample((1:num.obs)[-rand.o],1)
    nodes.rand.o2 <- as.logical(presence.tables[,rand.o2])
    current.size <- compute_size_neighborhood_p9(current.partitions[,rand.o],current.partitions[,rand.o2])
    sampled <- sample_new_partition_p9(current.partitions[,rand.o], current.partitions[,rand.o2], current.size)
    new.p <- sampled$new.partition
    nodes.changed <- sampled$nodes.changed
  } else if(move == "pair_next" && is.null(sizes.allowed)){
    if(rand.o == 1) rand.o2 <- 2
    if(rand.o == num.obs) rand.o2 <- num.obs-1
    if(rand.o > 1 && rand.o < num.obs) rand.o2 <- sample(c(rand.o-1,rand.o+1),1)
    nodes.rand.o2 <- as.logical(presence.tables[,rand.o2])
    current.size <- compute_size_neighborhood_p9(current.partitions[,rand.o],current.partitions[,rand.o2])
    sampled <- sample_new_partition_p9(current.partitions[,rand.o], current.partitions[,rand.o2], current.size)
    new.p <- sampled$new.partition
    nodes.changed <- sampled$nodes.changed
  }
  
  # Alternative if sizes are restricted
  if(move == "swap" && !is.null(sizes.allowed)) {
    current.size <- compute_size_neighborhood_p1_restricted(current.partitions[nodes.rand.o,rand.o], sizes.simulated)
    if(current.size$num.swaps > 0){
      sampled <- sample_new_partition_p1_restricted(current.partitions[nodes.rand.o,rand.o], mini.steps, sizes.simulated, current.size)
      new.p <- sampled$new.partition
      nodes.changed <- ((1:num.nodes)[nodes.rand.o])[sampled$nodes.changed]
    } else {
      new.p <- current.partitions[nodes.rand.o,rand.o]
      nodes.changed <- c()
    }
  } else if(move == "mergediv" && !is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p2_restricted(current.partitions[nodes.rand.o,rand.o], sizes.simulated)
    sampled <- sample_new_partition_p2_restricted(current.partitions[nodes.rand.o,rand.o], mini.steps, sizes.simulated, current.size)
    new.p <- sampled$new.partition
    nodes.changed <- ((1:num.nodes)[nodes.rand.o])[sampled$nodes.changed]
  } else if(move == "single" && !is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p3_restricted(current.partitions[nodes.rand.o,rand.o], sizes.simulated)
    sampled <- sample_new_partition_p3_restricted(current.partitions[nodes.rand.o,rand.o], mini.steps, sizes.simulated, current.size)
    new.p <- sampled$new.partition
    nodes.changed <- ((1:num.nodes)[nodes.rand.o])[sampled$nodes.changed]
  } else if(move == "double" && !is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p4_restricted(current.partitions[nodes.rand.o,rand.o], sizes.simulated)
    new.p <- sample_new_partition_p4_restricted(current.partitions[nodes.rand.o,rand.o], mini.steps, sizes.simulated, current.size)
  } else if(move == "swap2" && !is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p5_restricted(current.partitions[nodes.rand.o,rand.o], sizes.simulated)
    new.p <- sample_new_partition_p5_restricted(current.partitions[nodes.rand.o,rand.o], mini.steps, sizes.simulated, current.size)
  } else if(move == "exchange" && !is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p6_restricted(current.partitions[nodes.rand.o,rand.o], sizes.simulated)
    sampled <- sample_new_partition_p6_restricted(current.partitions[nodes.rand.o,rand.o], mini.steps, sizes.simulated, current.size)
    new.p <- sampled$new.partition
    nodes.changed <- ((1:num.nodes)[nodes.rand.o])[sampled$nodes.changed]
  } else if(move == "partition_exchange" && !is.null(sizes.allowed)){
    rand.o2 <- sample((1:num.obs)[-rand.o],1)
    nodes.rand.o2 <- as.logical(presence.tables[,rand.o2])
    current.size <- compute_size_neighborhood_p7(current.partitions[,rand.o],current.partitions[,rand.o2])
    new.ps <- sample_new_partition_p7(current.partitions[,rand.o], current.partitions[,rand.o2])
    new.p <- new.ps$new.partition1
    new.p2 <- new.ps$new.partition2
  } else if(move == "partition_next" && !is.null(sizes.allowed)){
    if(rand.o == 1) rand.o2 <- 2
    if(rand.o == num.obs) rand.o2 <- num.obs-1
    if(rand.o > 1 && rand.o < num.obs) rand.o2 <- sample(c(rand.o-1,rand.o+1),1)
    nodes.rand.o2 <- as.logical(presence.tables[,rand.o2])
    current.size <- compute_size_neighborhood_p7(current.partitions[,rand.o],current.partitions[,rand.o2])
    new.ps <- sample_new_partition_p7(current.partitions[,rand.o], current.partitions[,rand.o2])
    new.p <- new.ps$new.partition1
    new.p2 <- new.ps$new.partition2
  } else if(move == "group_exchange" && !is.null(sizes.allowed)){
    rand.o2 <- sample((1:num.obs)[-rand.o],1)
    nodes.rand.o2 <- as.logical(presence.tables[,rand.o2])
    current.size <- compute_size_neighborhood_p8(current.partitions[,rand.o],current.partitions[,rand.o2])
    new.p <- sample_new_partition_p8(current.partitions[,rand.o], current.partitions[,rand.o2], current.size)
  } else if(move == "group_next" && !is.null(sizes.allowed)){
    if(rand.o == 1) rand.o2 <- 2
    if(rand.o == num.obs) rand.o2 <- num.obs-1
    if(rand.o > 1 && rand.o < num.obs) rand.o2 <- sample(c(rand.o-1,rand.o+1),1)
    nodes.rand.o2 <- as.logical(presence.tables[,rand.o2])
    current.size <- compute_size_neighborhood_p8(current.partitions[,rand.o],current.partitions[,rand.o2])
    new.p <- sample_new_partition_p8(current.partitions[,rand.o], current.partitions[,rand.o2], current.size)
  } else if(move == "pair_exchange" && !is.null(sizes.allowed)){
    rand.o2 <- sample((1:num.obs)[-rand.o],1)
    nodes.rand.o2 <- as.logical(presence.tables[,rand.o2])
    current.size <- compute_size_neighborhood_p9_restricted(current.partitions[,rand.o],current.partitions[,rand.o2], sizes.simulated)
    sampled <- sample_new_partition_p9_restricted(current.partitions[,rand.o], current.partitions[,rand.o2], current.size, sizes.simulated)
    new.p <- sampled$new.partition
    nodes.changed <- sampled$nodes.changed
  } else if(move == "pair_next" && !is.null(sizes.allowed)){
    if(rand.o == 1) rand.o2 <- 2
    if(rand.o == num.obs) rand.o2 <- num.obs-1
    if(rand.o > 1 && rand.o < num.obs) rand.o2 <- sample(c(rand.o-1,rand.o+1),1)
    nodes.rand.o2 <- as.logical(presence.tables[,rand.o2])
    current.size <- compute_size_neighborhood_p9_restricted(current.partitions[,rand.o],current.partitions[,rand.o2], sizes.simulated)
    sampled <- sample_new_partition_p9_restricted(current.partitions[,rand.o], current.partitions[,rand.o2], current.size, sizes.simulated)
    new.p <- sampled$new.partition
    nodes.changed <- sampled$nodes.changed
  }
  
  #end <- Sys.time()
  #print(end-start)
  
  if(move == "partition_exchange" || move == "partition_next") {
    new.partitions[,rand.o] <- new.p
    new.partitions[,rand.o2] <- new.p2
  } else if(move == "group_exchange" || move == "group_next" || move == "pair_exchange" || move == "pair_next") {
    new.partitions[,rand.o] <- new.p
  } else {
    new.partitions[,rand.o] <- rep(NA,num.nodes)
    new.partitions[nodes.rand.o,rand.o] <- new.p
  }
  
  # IF NO CHANGE, same statistics
  if(all(current.partitions == new.partitions, na.rm=T)) {
    new.z.contributions <- current.z.contributions
    new.z <- current.z
  }
  
  #start <- Sys.time()
  
  # IF CHANGE: compute new statistics
  if(!all(current.partitions == new.partitions, na.rm=T)) {
    if(move == "partition_exchange" || move == "partition_next"){
      recalculated.stats <- recalculate_statistics(new.partitions, rand.o, nodes.rand.o, nodes, effects, objects, current.z.contributions) 
      new.z.contributions <- recalculated.stats$new.z.contributions
      new.z <- recalculated.stats$new.z
      recalculated.stats <- recalculate_statistics(new.partitions, rand.o2, nodes.rand.o2, nodes, effects, objects, new.z.contributions) 
      new.z.contributions <- recalculated.stats$new.z.contributions
      new.z <- recalculated.stats$new.z
    } else {
      recalculated.stats <- recalculate_statistics(current.partitions, new.partitions, rand.o, nodes.changed, nodes, effects, objects, current.z.contributions) 
      new.z.contributions <- recalculated.stats$new.z.contributions
      new.z <- recalculated.stats$new.z
    }
    
  } 
  
  #end <- Sys.time()
  #print(end-start)
  
  
  new.logit <- theta * new.z
  
  #start <- Sys.time()
  
  # chose whether to change or not
  if(mini.steps == "normalized") {
    
    neighborhoods.ratios <- 1
    
    if(move == "swap" && is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p1(new.partitions[nodes.rand.o,rand.o])
        neighborhoods.ratios <- neighborhoods.ratios * current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratios <- 1
      }
    } else if(move == "mergediv" && is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p2(new.partitions[nodes.rand.o,rand.o])
      neighborhoods.ratios <- neighborhoods.ratios * (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
    } else if(move == "single" && is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p3(new.partitions[nodes.rand.o,rand.o])
        neighborhoods.ratios <- neighborhoods.ratios * current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratios <- 1
      }
    } else if(move == "double" && is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p4(new.partitions[nodes.rand.o,rand.o])
        neighborhoods.ratios <- neighborhoods.ratios * current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratios <- 1
      }
    } else if(move == "swap2" && is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p5(new.partitions[nodes.rand.o,rand.o])
        neighborhoods.ratios <- neighborhoods.ratios * current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratios <- 1
      }
    } else if(move == "exchange" && is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p6(new.partitions[nodes.rand.o,rand.o])
        neighborhoods.ratios <- neighborhoods.ratios * current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratios <- 1
      }
    } else if(move == "partition_exchange" && is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p7(new.partitions[,rand.o],new.partitions[,rand.o2])
      neighborhoods.ratios <- neighborhoods.ratios * current.size / new.size
    } else if(move == "partition_next" && is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p7(new.partitions[,rand.o],new.partitions[,rand.o2])
      neighborhoods.ratios <- neighborhoods.ratios * current.size / new.size
    } else if(move == "group_exchange" && is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p8(new.partitions[,rand.o],new.partitions[,rand.o2])
      neighborhoods.ratios <- 1 * current.size$num.reorganizations / new.size$num.reorganizations
    } else if(move == "group_next" && is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p8(new.partitions[,rand.o],new.partitions[,rand.o2])
      neighborhoods.ratios <- 1 * current.size$num.reorganizations / new.size$num.reorganizations
    } else if(move == "pair_exchange" && is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p9(new.partitions[,rand.o],new.partitions[,rand.o2])
      neighborhoods.ratios <- 1 * current.size$num.swaps / new.size$num.swaps
    } else if(move == "pair_next" && is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p9(new.partitions[,rand.o],new.partitions[,rand.o2])
      neighborhoods.ratios <- 1 * current.size$num.swaps / new.size$num.swaps
    }
    
    # Alternative if sizes are restricted
    if(move == "swap" && !is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p1_restricted(new.partitions[nodes.rand.o,rand.o], sizes.simulated)
        neighborhoods.ratios <- neighborhoods.ratios * current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratios <- 1
      }
    } else if(move == "mergediv" && !is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p2_restricted(new.partitions[nodes.rand.o,rand.o], sizes.simulated)
      neighborhoods.ratios <- neighborhoods.ratios * (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
    } else if(move == "single" && !is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p3_restricted(new.partitions[nodes.rand.o,rand.o], sizes.simulated)
        neighborhoods.ratios <- neighborhoods.ratios * current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratios <- 1
      }
    } else if(move == "double" && !is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p4_restricted(new.partitions[nodes.rand.o,rand.o], sizes.simulated)
        neighborhoods.ratios <- neighborhoods.ratios * current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratios <- 1
      }
    } else if(move == "swap2" && !is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p5_restricted(new.partitions[nodes.rand.o,rand.o], sizes.simulated)
        neighborhoods.ratios <- neighborhoods.ratios * current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratios <- 1
      }
    } else if(move == "exchange" && !is.null(sizes.allowed)) {
      if(current.size$num.swaps > 0){
        new.size <- compute_size_neighborhood_p6_restricted(new.partitions[nodes.rand.o,rand.o], sizes.simulated)
        neighborhoods.ratios <- neighborhoods.ratios * current.size$num.swaps / new.size$num.swaps 
      } else {
        neighborhoods.ratios <- 1
      }
    } else if(move == "partition_exchange" && !is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p7(new.partitions[,rand.o],new.partitions[,rand.o2])
      neighborhoods.ratios <- neighborhoods.ratios * current.size / new.size
    } else if(move == "partition_next" && !is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p7(new.partitions[,rand.o],new.partitions[,rand.o2])
      neighborhoods.ratios <- neighborhoods.ratios * current.size / new.size
    } else if(move == "group_exchange" && !is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p8(new.partitions[,rand.o],new.partitions[,rand.o2])
      neighborhoods.ratios <- 1 * current.size$num.reorganizations / new.size$num.reorganizations
    } else if(move == "group_next" && !is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p8(new.partitions[,rand.o],new.partitions[,rand.o2])
      neighborhoods.ratios <- 1 * current.size$num.reorganizations / new.size$num.reorganizations
    } else if(move == "pair_exchange" && !is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p9_restricted(new.partitions[,rand.o],new.partitions[,rand.o2], sizes.simulated)
      neighborhoods.ratios <- 1 * current.size$num.swaps / new.size$num.swaps
    } else if(move == "group_next" && !is.null(sizes.allowed)) {
      new.size <- compute_size_neighborhood_p9_restricted(new.partitions[,rand.o],new.partitions[,rand.o2], sizes.simulated)
      neighborhoods.ratios <- 1 * current.size$num.swaps / new.size$num.swaps
    }
    
    hastings.ratio <- (exp(sum(new.logit) - sum(current.logit))) * neighborhoods.ratios
    
  } else if(mini.steps == "selfloops"){
    
    hastings.ratio <- exp(sum(new.logit) - sum(current.logit))
    
  }
  
  #end <- Sys.time()
  #print(end-start)
  
  return(list("hastings.ratio" = hastings.ratio,
              "new.partitions" = new.partitions,
              "new.z.contributions" = new.z.contributions,
              "new.z" = new.z,
              "new.logit" = new.logit,
              "move" = move))
}



# Recalculate statistics for a changed partition in multiple estimation
recalculate_statistics <- function(current.partitions, new.partitions, rand.o, nodes.rand.o, nodes, effects, objects, current.z.contributions) {
  
  # store new statistics
  new.z.contributions <- current.z.contributions
  
  # calculate separately each effect
  for(e in 1:length(effects$names)) {
    
    object.name <- effects$objects[e]
    
    # TIE EFFECT: keep only present nodes
    if(effects$names[e] == "tie"){
      
      effects.temp <- list(names="tie",objects="net.temp")
      for(ob in 1:length(objects)){
        if(objects[[ob]][[1]] == object.name){
          net <- objects[[ob]][[2]]
          objects.temp <- list(list(name="net.temp",object=net[nodes.rand.o,nodes.rand.o]))
        }
      }
      nodes.temp <- nodes[nodes.rand.o,]
      
      old.z.temp <- computeStatistics(order_groupids(current.partitions[nodes.rand.o,rand.o]), nodes.temp, effects.temp, objects.temp)
      new.z.temp <- computeStatistics(order_groupids(new.partitions[nodes.rand.o,rand.o]), nodes.temp, effects.temp, objects.temp)
      new.z.contributions[e,rand.o] <- new.z.contributions[e,rand.o] - old.z.temp + new.z.temp
      #new.z.contributions[e,rand.o] <- as.numeric(computeStatistics(new.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp))
      
      #  SAME_VAR EFFECT: adapt the varying covariate  
    } else if(effects$names[e] == "same_var"){
      
      effects.temp <- list(names="same",objects="var.temp")
      for(ob in 1:length(objects)){
        if(objects[[ob]][[1]] == object.name){
          atts <- objects[[ob]][[2]]
        }
      }
      objects.temp <- list()
      nodes.temp <- nodes[nodes.rand.o,]
      nodes.temp$var.temp <- atts[nodes.rand.o,rand.o]
      
      old.z.temp <- computeStatistics(order_groupids(current.partitions[nodes.rand.o,rand.o]), nodes.temp, effects.temp, objects.temp)
      new.z.temp <- computeStatistics(order_groupids(new.partitions[nodes.rand.o,rand.o]), nodes.temp, effects.temp, objects.temp)
      new.z.contributions[e,rand.o] <- new.z.contributions[e,rand.o] - old.z.temp + new.z.temp
      #new.z.contributions[e,rand.o] <- as.numeric(computeStatistics(new.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp))
      
      
      #  TIE_VAR EFFECT: adapt the varying network
    } else if(effects$names[e] == "tie_var"){
      
      effects.temp <- list(names="tie",objects="net.temp")
      for(ob in 1:length(objects)){
        if(objects[[ob]][[1]] == object.name){
          nets <- objects[[ob]][[2]]
          net <- nets[[rand.o]]
          objects.temp <- list(list(name="net.temp",object=net[nodes.rand.o,nodes.rand.o]))
        }
      }
      nodes.temp <- nodes[nodes.rand.o,]
      
      old.z.temp <- computeStatistics(order_groupids(current.partitions[nodes.rand.o,rand.o]), nodes.temp, effects.temp, objects.temp)
      new.z.temp <- computeStatistics(order_groupids(new.partitions[nodes.rand.o,rand.o]), nodes.temp, effects.temp, objects.temp)
      new.z.contributions[e,rand.o] <- new.z.contributions[e,rand.o] - old.z.temp + new.z.temp
      #new.z.contributions[e,rand.o] <- as.numeric(computeStatistics(new.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp))
      
      #  TIE_VAR_X_DIFF EFFECT: adapt the varying network
    } else if(effects$names[e] == "tie_var_X_diff"){
      
      effects.temp <- list(names="tie_X_diff",objects="net.temp",objects2=effects$objects2[e])
      for(ob in 1:length(objects)){
        if(objects[[ob]][[1]] == object.name){
          nets <- objects[[ob]][[2]]
          net <- nets[[rand.o]]
          objects.temp <- list(list(name="net.temp",object=net[nodes.rand.o,nodes.rand.o]))
        }
      }
      nodes.temp <- nodes[nodes.rand.o,]
      
      old.z.temp <- computeStatistics(order_groupids(current.partitions[nodes.rand.o,rand.o]), nodes.temp, effects.temp, objects.temp)
      new.z.temp <- computeStatistics(order_groupids(new.partitions[nodes.rand.o,rand.o]), nodes.temp, effects.temp, objects.temp)
      new.z.contributions[e,rand.o] <- new.z.contributions[e,rand.o] - old.z.temp + new.z.temp
      #new.z.contributions[e,rand.o] <- as.numeric(computeStatistics(new.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp))
      
      
      
    } else {
      
      effects.temp <- list(names=effects$names[e],objects=effects$objects[e])
      objects.temp <- objects
      nodes.temp <- nodes[nodes.rand.o,]
      
      old.z.temp <- computeStatistics(order_groupids(current.partitions[nodes.rand.o,rand.o]), nodes.temp, effects.temp, objects.temp)
      new.z.temp <- computeStatistics(order_groupids(new.partitions[nodes.rand.o,rand.o]), nodes.temp, effects.temp, objects.temp)
      new.z.contributions[e,rand.o] <- new.z.contributions[e,rand.o] - old.z.temp + new.z.temp
      #new.z.contributions[e,rand.o] <- as.numeric(computeStatistics(new.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp))
      
    }
  }
  
  new.z <- rowSums(new.z.contributions)
  
  return(list(new.z.contributions = new.z.contributions,
              new.z = new.z))
  
}



## NEIGHBORHOOD PI 1: only swaps of two nodes (careful, cannot be used alone)
# compute_size_neighborhood_p1 <- function(partition){
#   
#   # find isolates, pairs and groups>2
#   isolates <- c()
#   pairs <- c()
#   others <- c()
#   num.groups <- max(partition)
#   num.nodes <- length(partition)
#   
#   for(g in 1:num.groups){
#     if(length(which(partition == g)) == 1) {
#       isolates <- c(isolates,g)
#     } else if(length(which(partition == g)) == 2) {
#       pairs <- c(pairs,g)
#     } else {
#       others <- c(others,g)
#     }
#   }
#   
#   nums.swaps <- matrix(0,num.nodes,num.nodes)
#   for(i in 1:(num.nodes-1)){
#     for(j in (i+1):num.nodes){
#       g_i <- partition[i]
#       g_j <- partition[j]
#       
#       allowed <- T
#       if(g_i == g_j) allowed <- F # two members of the same group cannot swap, it wouldn't change anything
#       if(g_i %in% isolates && g_j %in% isolates){
#         allowed <- F # two isolates cannot swap, it wouldn't change anything
#       }
#       
#       if(allowed) nums.swaps[i,j] <- 1
#     }
#   }
#   
#   num.swaps <- sum(nums.swaps)
#   
#   return(list(isolates = isolates,
#               pairs = pairs,
#               others = others,
#               nums.swaps = nums.swaps,
#               num.swaps = num.swaps))
#   
# }
# 
# sample_new_partition_p1 <- function(current.partition, mini.steps, size_neighborhood){
#   
#   # calculate the number of neighbor partitions
#   num.nodes <- length(current.partition)
#   num.groups <- max(current.partition)
#   size <- size_neighborhood
#   isolates <- size$isolates
#   pairs <- size$pairs
#   others <- size$others
#   nums.swaps <- size$nums.swaps
#   num.swaps <- size$num.swaps
#   
#   if(mini.steps == "normalized") {
#     total <- num.swaps
#   } 
#   
#   pick.1 <- sample(total,1)
#   
#   # decide which actors to swap
#   new.partition <- current.partition
#   swap <- which(nums.swaps == 1)[pick.1]
#   j <- ((swap-1) %/% num.nodes) + 1
#   i <- (swap-1) %% num.nodes + 1
#   new.partition[i] <- current.partition[j]
#   new.partition[j] <- current.partition[i]
#   new.partition <- order_groupids(new.partition)
#   
#   if(all(new.partition == current.partition)){
#     print("test")
#   }
#   
#   return(new.partition)
#   
# }



compute_size_neighborhood_p1 <- function(partition){
  
  return(list(num.swaps =1))
  
}

sample_new_partition_p1 <- function(current.partition, mini.steps, size_neighborhood){
  
  n <- size_neighborhood
  num.nodes <- length(current.partition)
  picks <- sample(1:num.nodes,2)
  new.partition <- current.partition
  new.partition[picks[1]] <- current.partition[picks[2]]
  new.partition[picks[2]] <- current.partition[picks[1]]
  new.partition <- order_groupids(new.partition)
  
  nodes.changed <- which(current.partition == current.partition[picks[1]] | current.partition == current.partition[picks[2]])
  
  return(list(new.partition = new.partition,
              nodes.changed = nodes.changed))
  
}

## NEIGHBORHOOD PI 2: only merges and divisions of 2 groups
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


sample_new_partition_p2 <- function(current.partition, mini.steps, size_neighborhood){
  
  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- size_neighborhood
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
    new.partition <- order_groupids(new.partition)
  }
  
  return(new.partition)
  
}


## NEIGHBORHOOD PI 3: only swaps of one node 
compute_size_neighborhood_p3 <- function(partition){
  
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
  done.pairs <- rep(0,length(pairs))
  done.isolates <- rep(0,length(isolates))
  for(i in 1:num.nodes){
    g <- partition[i]
    if(g %in% isolates){
      done.isolates[isolates == g] <- 1
      nums.swaps[i] <- num.groups - sum(done.isolates) # can join any other group, except other isolates before (because already done) and itself
    }
    if(g %in% pairs){
      nums.swaps[i] <- num.groups - done.pairs[pairs == g] # can join any other group or create its isolate, except if the possibilty of splitting the pair in 2 isolates is already counted
      done.pairs[pairs == g] <- 1
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

sample_new_partition_p3 <- function(current.partition, mini.steps, size_neighborhood){
  
  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- size_neighborhood
  isolates <- size$isolates
  pairs <- size$pairs
  others <- size$others
  nums.swaps <- size$nums.swaps
  num.swaps <- size$num.swaps
  
  if(mini.steps == "normalized") {
    total <- num.swaps
  } 
  
  pick.1 <- sample(total,1)
  
  # decide which actor to swap
  all.groups <- 1:num.groups
  new.partition <- current.partition
  done.pairs <- rep(0,length(pairs))
  done.isolates <- rep(0,length(isolates))
  
  for(i in 1:num.nodes){
    
    if(i == 1) {start <- 0} else {start <- sum(nums.swaps[1:i-1])}
    if(i == num.nodes) {end <- num.swaps} else {end <- sum(nums.swaps[1:i])}
    
    g <- current.partition[i]
    if(g %in% isolates) done.isolates[isolates == g] <- 1
    
    if(pick.1 > start && pick.1 <= end) {
      g <- current.partition[i]
      
      # an isolate is randomly joining another group except isolates before i
      if(g %in% isolates){
        previousisolates <- isolates[done.isolates == 1]
        tosample <- all.groups[!all.groups %in% previousisolates]
        if(length(tosample) == 1) newg <- tosample
        if(length(tosample) >= 2) newg <- sample(tosample,1)
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }
      
      # a member of a pair is randomly joining another group 
      # or creating its isolate if it's the first of the pair to appear
      if(g %in% pairs){
        if(done.pairs[pairs == g] == 0){
          newg <- sample(all.groups,1)
          if(newg == g){
            newg <- num.groups + 1
          }
        } else {
          tosample <- all.groups[all.groups != g]
          if(length(tosample) == 1) newg <- tosample
          if(length(tosample) >= 2) newg <- sample(tosample,1)
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
    
    if(g %in% pairs) done.pairs[pairs == g] <- 1  
  }
  
  return(new.partition)
  
}


## NEIGHBORHOOD PI 4: swaps a pair of nodes
compute_size_neighborhood_p4 <- function(partition){
  
  # find isolates, pairs and groups>2
  isolates <- c()
  pairs <- c()
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
  
  nums.swaps <- rep(0,max(partition))
  done.pairs <- rep(0,length(pairs))
  for(g in 1:num.groups){
    
    if(g %in% isolates){
      # nothing happens, there is no pair to swap here
    }
    if(g %in% pairs){
      nums.swaps[g] <- length(pairs) + length(others) - 1 - sum(done.pairs) + done.pairs[pairs==g]
      done.pairs[pairs == g] <- 1
    } 
    if(g %in% others){
      sizeg <- sum(partition == g)
      nums.swaps[g] <- (length(pairs) + length(others)-1) * sizeg * (sizeg-1) / 2
    } 
  }
  
  num.swaps <- sum(nums.swaps)
  
  return(list(isolates = isolates,
              pairs = pairs,
              others = others,
              nums.swaps = nums.swaps,
              num.swaps = num.swaps))
  
}

sample_new_partition_p4 <- function(current.partition, mini.steps, size_neighborhood){
  
  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- size_neighborhood
  isolates <- size$isolates
  pairs <- size$pairs
  others <- size$others
  nums.swaps <- size$nums.swaps
  num.swaps <- size$num.swaps
  
  if(mini.steps == "normalized") {
    total <- num.swaps
  } 
  
  pick.1 <- sample(total,1)
  
  # decide which group to pick a pair from
  all.groups <- 1:num.groups
  new.partition <- current.partition
  done.pairs <- rep(0,length(pairs))
  
  for(g in 1:num.groups){
    
    if(g == 1) {start <- 0} else {start <- sum(nums.swaps[1:g-1])}
    if(g == num.groups) {end <- num.swaps} else {end <- sum(nums.swaps[1:g])}
    
    if(pick.1 > start && pick.1 <= end) {
      
      # a pair can join other available groups, unless it's a pair that has already been "joined"
      if(g %in% pairs){
        members <- which(current.partition == g)
        previous.pairs <- pairs[done.pairs == 1]
        tosample <- all.groups[(all.groups %in% c(pairs,others)) & (all.groups != g) & !(all.groups %in% previous.pairs)]
        if(length(tosample) == 1) newg <- tosample
        if(length(tosample) >= 2) newg <- sample(tosample,1)
        new.partition[members] <- newg
        new.partition <- order_groupids(new.partition)
      }
      
      # a pair in a bigger group can join any other group 
      if(g %in% others && g %in% abovemin_groups){
        members <- which(current.partition == g)
        sizeg <- length(members)
        pairg <- sample(members,2)
        tosample <- all.groups[(all.groups %in% c(pairs,others)) & (all.groups != g)]
        if(length(tosample) == 1) newg <- tosample
        if(length(tosample) >= 2) newg <- sample(tosample,1)
        new.partition[pairg] <- newg
        new.partition <- order_groupids(new.partition)
      }
      
    }
    
    if(g %in% pairs) done.pairs[pairs == g] <- 1
  }
  
  return(new.partition)
  
}


## NEIGHBORHOOD PI 5: only swaps of two pairs of nodes (careful, cannot be used alone)
compute_size_neighborhood_p5 <- function(partition){
  
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
  
  aff <- as.matrix(table(data.frame(actor = 1:num.nodes, group = partition)))
  net <- aff %*% t(aff)
  net[upper.tri(net)] <- 0
  diag(net) <- 0
  paired.actors <- which(net == 1)
  num.paired.actors <- sum(net)
  
  groups.paired.actors <- rep(0,num.paired.actors)
  pairs.paired.actors <- rep(0,num.paired.actors)
  
  impossible.swaps <- 0
  for(g in 1:num.groups){
    members <- which(partition == g)
    
    if(g %in% pairs){
      i <- members[1]
      j <- members[2]
      p <- (i-1)*num.nodes+j
      pairs.paired.actors[which(paired.actors == p)] <- 1
      groups.paired.actors[which(paired.actors == p)] <- g
    }
    
    if(g %in% others){
      for(i in 1:(length(members)-1)){
        for(j in (i+1):length(members)){
          p <- (members[i]-1)*num.nodes+members[j]
          groups.paired.actors[which(paired.actors == p)] <- g
        }
      }
      num.pairs.g <- length(members)*(length(members)-1)/2
      impossible.swaps <- impossible.swaps + num.pairs.g*(num.pairs.g-1)/2
    }
    
  }
  
  num.swaps <- num.paired.actors * (num.paired.actors-1) / 2 - impossible.swaps - sum(pairs.paired.actors) * (sum(pairs.paired.actors)-1) / 2
  
  return(list(isolates = isolates,
              pairs = pairs,
              others = others,
              paired.actors = paired.actors,
              groups.paired.actors = groups.paired.actors,
              pairs.paired.actors = pairs.paired.actors,
              num.swaps = num.swaps))
  
}

sample_new_partition_p5 <- function(current.partition, mini.steps, size_neighborhood){
  
  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- size_neighborhood
  isolates <- size$isolates
  pairs <- size$pairs
  others <- size$others
  paired.actors <- size$paired.actors
  groups.paired.actors <- size$groups.paired.actors
  pairs.paired.actors <- size$pairs.paired.actors
  num.swaps <- size$num.swaps
  
  if(mini.steps == "normalized") {
    total <- num.swaps
  } 
  
  # decide which actors to swap
  new.partition <- current.partition
  found <- F
  while(!found){
    pick.1 <- sample(paired.actors,1)
    pick.2 <- sample(paired.actors,1)
    if(groups.paired.actors[which(paired.actors == pick.1)] != groups.paired.actors[which(paired.actors == pick.2)] &&
       !(pairs.paired.actors[which(paired.actors == pick.1)] && pairs.paired.actors[which(paired.actors == pick.2)])) {
      found <- T
    }
  }
  
  j1 <- ((pick.1-1) %/% num.nodes) + 1
  i1 <- (pick.1-1) %% num.nodes + 1
  j2 <- ((pick.2-1) %/% num.nodes) + 1
  i2 <- (pick.2-1) %% num.nodes + 1
  new.partition[i1] <- current.partition[i2]
  new.partition[j1] <- current.partition[j2]
  new.partition[i2] <- current.partition[i1]
  new.partition[j2] <- current.partition[j1]
  new.partition <- order_groupids(new.partition)
  
  return(new.partition)
  
}

## NEIGHBORHOOD PI 6: pick two groups, re attribute nodes within them (careful, cannot be used alone)
compute_size_neighborhood_p6 <- function(partition){
  
  # find isolates, pairs and groups>2
  #isolates <- c()
  #pairs <- c()
  #others <- c()
  num.groups <- max(partition)
  num.nodes <- length(partition)
  sizes <- table(partition)
  
  # for(g in 1:num.groups){
  #   if(length(which(partition == g)) == 1) {
  #     isolates <- c(isolates,g)
  #   } else if(length(which(partition == g)) == 2) {
  #     pairs <- c(pairs,g)
  #   } else {
  #     others <- c(others,g)
  #   }
  # }
  isolates <- which(sizes == 1)
  pairs <- which(sizes == 2)
  others <- which(sizes > 2)
  
  group.pairs <- matrix(0,num.groups,num.groups)
  nums.swaps <- matrix(0,num.groups,num.groups)
  
  for(g1 in 1:(num.groups-1)){
    
    for(g2 in (g1+1):num.groups){
      size.1 <- sum(partition == g1)
      size.2 <- sum(partition == g2)
      
      k <- size.1
      n <- size.1 + size.2
      
      nums.swaps[g1,g2] <- factorial(n) / (factorial(k)*factorial(n-k))
    }
    
  }
  num.swaps <- sum(nums.swaps)
  
  return(list(isolates = isolates,
              pairs = pairs,
              others = others,
              nums.swaps = nums.swaps,
              num.swaps = num.swaps))
  
}

sample_new_partition_p6 <- function(current.partition, mini.steps, size_neighborhood){
  
  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- size_neighborhood
  isolates <- size$isolates
  pairs <- size$pairs
  others <- size$others
  nums.swaps <- size$nums.swaps
  num.swaps <- size$num.swaps
  
  if(mini.steps == "normalized") {
    total <- num.swaps
  } 
  
  # decide which actors to swap
  new.partition <- current.partition
  pick <- sample(1:num.swaps,1)
  start <- 0
  end <- 0
  for(g1 in 1:(num.groups-1)){
    for(g2 in (g1+1):num.groups){
      
      end <- end + nums.swaps[g1,g2]
      
      if(start < pick && pick <= end){
        tosample <- c(which(current.partition==g1),which(current.partition==g2))
        new.g1 <- sample(tosample,sum(current.partition==g1))
        new.g2 <- tosample[-which(tosample %in% new.g1)]
        
        new.partition[new.g1] <- g1
        new.partition[new.g2] <- g2
        new.partition <- order_groupids(new.partition)
        
        nodes.changed <- which(current.partition == g1 | current.partition == g2)
      }
      
      start <- start + nums.swaps[g1,g2]
    }
  }
  
  return(list(new.partition = new.partition,
              nodes.changed = nodes.changed))
  
}



## NEIGHBORHOOD PI 7: pick two partitions and exchange then (only for multiple partition estimation)
# since we force one change, and it's reversible (bijection), we don't have to care about the size of the neighborhood
# CAREFUL, only works when singletons are allowed
compute_size_neighborhood_p7 <- function(partition1, partition2){
  
  num.nodes <- length(partition1)
  present1 <- !is.na(partition1)
  present2 <- !is.na(partition2)
  
  # calculate the number of ways to reconstruct p1 after having done the exchange with p2
  p1withoutabsentp2 <- partition1[present1 & present2]
  max1 <- length(unique(p1withoutabsentp2))
  nodesunknown1 <- sum(present1 & !present2)
  if(nodesunknown1 == 0) {
    possibilities1 <- 1
  } else {
    possibilities1 <- 0
    for(k in 1:nodesunknown1) {
      possibilities1 <- possibilities1 + Stirling(nodesunknown1,k)*(max1+1)^k
    }
  }
  
  # same the other way around
  p2withoutabsentp1 <- partition2[present2 & present1]
  max2 <- length(unique(p2withoutabsentp1))
  nodesunknown2 <- sum(present2 & !present1)
  if(nodesunknown2 == 0) {
    possibilities2 <- 1
  } else {
    possibilities2 <- 0
    for(k in 1:nodesunknown2) {
      possibilities2 <- possibilities2 + Stirling(nodesunknown2,k)*(max2+1)^k
    }
  }
  
  return(possibilities1 * possibilities2)
}

sample_new_partition_p7 <- function(current.partition1, current.partition2){
  
  num.nodes <- length(current.partition1)
  present1 <- which(!is.na(current.partition1))
  present2 <- which(!is.na(current.partition2))
  
  new.partition1 <- rep(0,num.nodes)
  new.partition2 <- rep(0,num.nodes)
  
  # fill in the first with the second, and create isolates for people not present in the second
  new.partition1[present1] <- current.partition2[present1]
  max1 <- max(new.partition1, na.rm=T)
  singletons1 <- which(is.na(new.partition1))
  if(length(singletons1) > 0) { # assign randomly the ones who are not present in the other one
    for(s in singletons1){
      newgroup <- sample(1:(max1+1),1)
      new.partition1[s] <- newgroup
      if(newgroup == (max1+1)) max1 <- max1 + 1
    }
  }
  new.partition1[present1] <- order_groupids(new.partition1[present1])
  new.partition1[new.partition1 == 0] <- NA
  
  # same the other way around
  new.partition2[present2] <- current.partition1[present2]
  max2 <- max(new.partition2, na.rm=T)
  singletons2 <- which(is.na(new.partition2))
  if(length(singletons2) > 0) { # assign randomly the ones who are not present in the other one
    for(s in singletons2){
      newgroup <- sample(1:(max2+1),1)
      new.partition2[s] <- newgroup
      if(newgroup == (max2+1)) max2 <- max2 + 1
    }
  }
  new.partition2[present2] <- order_groupids(new.partition2[present2])
  new.partition2[new.partition2 == 0] <- NA
  
  return(list(new.partition1 = new.partition1,
              new.partition2 = new.partition2))
  
}


## NEIGHBORHOOD PI 8: pick another partition, a group in it, and reproduce it in the current partition, or reorganizes it
# CAREFUL, only works when singletons are allowed
compute_size_neighborhood_p8 <- function(partition1, partition2){
  
  num.nodes <- length(partition1)
  present1 <- !is.na(partition1)
  present2 <- !is.na(partition2)
  
  # count number of groups from the second partition
  groups.toreorganize <- unique(partition2[present1 & present2])
  nums.reorganizations <- rep(0,length(groups.toreorganize))
  
  # count for each of them number of moves
  for(g in 1:length(groups.toreorganize)){
    group <- groups.toreorganize[g]
    members <- which(partition2 == group)
    max1 <- length(unique(partition1[present1 & partition1 != group]))
    for(k in 1:length(members)) {
      nums.reorganizations[g] <- nums.reorganizations[g] + Stirling(length(members),k)*(max1+1)^k
    }
    nums.reorganizations[g] <- nums.reorganizations[g] + 1 # when we just copy the group
  }
  
  return(list(groups.toreorganize = groups.toreorganize,
              nums.reorganizations = nums.reorganizations,
              num.reorganizations = sum(nums.reorganizations)))
}

sample_new_partition_p8 <- function(current.partition1, current.partition2, size_neighborhood){
  
  groups.toreorganize <- size_neighborhood$groups.toreorganize
  nums.reorganizations <- size_neighborhood$nums.reorganizations
  num.reorganizations <- size_neighborhood$num.reorganizations
  
  num.nodes <- length(current.partition1)
  present1 <- !is.na(current.partition1)
  present2 <- !is.na(current.partition2)
  
  new.partition <- current.partition1
  
  # pick 
  probas_groups <- nums.reorganizations / num.reorganizations
  pick <- sample(c(0,groups.toreorganize),1,prob=c(1/ num.reorganizations, probas_groups))
  
  #randomize group
  if(pick > 0) {
    group <- pick
    members <- which(current.partition2 == group && present1)
    max1 <- length(unique(current.partition1[present1 & current.partition1 != group]))
    for(a in members){
      newgroup <- sample(1:(max1+1),1)
      new.partition[a] <- newgroup
      if(newgroup == (max1+1)) max1 <- max1 + 1
    }
    new.partition[present1] <- order_groupids(new.partition[present1])
  }
  
  # copy group
  if(pick ==0) {
    # reproduce the group
    new.partition[members] <- max(new.partition,na.rm=T) + 1
    #reorder
    new.partition[present1] <- order_groupids(new.partition[present1])
  }
  
  return(new.partition)
  
}

#sample_new_partition_p8 <- function(current.partition1, current.partition2){
#  
#  num.nodes <- length(current.partition1)
#  present1 <- !is.na(current.partition1)
#  present2 <- !is.na(current.partition2)
#  
#  new.partition <- current.partition1
#  
#  # pick a group to reproduce
#  topick <- unique(current.partition2[present1 & present2])
#  g <- sample(topick,1)
#  members <- which(current.partition2 == g && present1)
#  
#  # reproduce the group
#  new.partition[members] <- max(new.partition,na.rm=T) + 1
#  
#  #reorder
#  new.partition[present1] <- order_groupids(new.partition[present1])
#  
#  return(new.partition)
#  
#}


## NEIGHBORHOOD PI 9: take two nodes together in partition 2, bring them together (single move)
## otherwise separate then
## careful need isolates
compute_size_neighborhood_p9 <- function(partition1, partition2){
  
  # find maximum size allowed
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)
  
  num.nodes <- length(partition1)
  present1 <- !is.na(partition1)
  present2 <- !is.na(partition2)
  
  num.nodes <- length(partition1)
  num.groups <- max(partition1,na.rm=T)
  sizes <- table(partition1)
  
  affiliation1 <- as.matrix(table(data.frame(actor = 1:num.nodes, group= partition1)))
  affiliation2 <- as.matrix(table(data.frame(actor = 1:num.nodes, group= partition2)))
  adjacency1 <- affiliation1 %*% t(affiliation1)
  adjacency2 <- affiliation2 %*% t(affiliation2)
  
  nums.swaps <- matrix(0,num.nodes,num.nodes)
  for(i in 1:num.nodes) {
    for(j in 1:num.nodes){
      ad1 <- adjacency1[i,j]
      ad2 <- adjacency2[i,j]
      if(!is.na(ad1) && !is.na(ad2) && ad2 == 1) {
        if(ad1 == 0) {
          nums.swaps[i,j] <- 1
        }
        if(ad1 == 1) {
          nums.swaps[i,j] <- num.groups
        }
      }
    }
  } 
  
  return(list(adjacency1 = adjacency1,
              nums.swaps = nums.swaps,
              num.swaps = sum(nums.swaps)))
}

sample_new_partition_p9 <- function(current.partition1, current.partition2, size_neighborhood){
  
  # find maximum size allowed
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)
  
  num.nodes <- length(current.partition1)
  num.groups <- max(current.partition1,na.rm=T)
  adjacency1 <- size_neighborhood$adjacency1
  nums.swaps <- size_neighborhood$nums.swaps
  num.swaps <- size_neighborhood$num.swaps
  
  present1 <- !is.na(current.partition1)
  new.partition1 <- current.partition1
  
  probas_pairs <- as.vector(nums.swaps) / num.swaps
  pick <- sample(1:(num.nodes*num.nodes),1,prob=probas_pairs)
  
  i <- (pick-1) %% num.nodes + 1
  j <- ((pick-1) %/% num.nodes) + 1
  
  if(adjacency1[i,j] == 0) {
    new.partition1[i] <- current.partition1[j]
  }
  if(adjacency1[i,j] == 1) {
    pickgroup <- sample(1:num.groups,1)
    if(pickgroup == current.partition1[i]) new.partition1[i] <- num.groups + 1
    else new.partition1[i] <- pickgroup
  }
  new.partition1[present1] <- order_groupids(new.partition1[present1])
  
  return(new.partition1)
  
}



## NEIGHBORHOOD PI 1: only swaps of two nodes (careful, cannot be used alone)
# BUT REMOVE OUT OF SAMPLE PARTITIONS
# WARNING!!!!: for now it only works if sizes allowed are an interval ([size_min,size_max])
compute_size_neighborhood_p1_restricted <- function(partition, sizes.simulated){
  # check if current partition is allowed
  if(!check_sizes(partition, sizes.simulated)) stop("The partition we are in is not allowed.")
  return(compute_size_neighborhood_p1(partition))
}

sample_new_partition_p1_restricted <- function(current.partition, mini.steps, sizes.simulated, size_neighborhood){
  return(sample_new_partition_p1(current.partition, mini.steps, size_neighborhood))
}


## NEIGHBORHOOD PI 2: only merges and divisions of 2 groups
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
  sizes <- table(partition)
  #num.merges <- 0
  #merges <- matrix(0,num.groups,num.groups)
  #nums.divisions <- rep(0,num.groups)
  #num.divisions <- 0
  
  # merges
  if(num.groups>1){
    # for(g1 in 1:(num.groups-1)){
    #   for(g2 in (g1+1):num.groups){
    #     if(length(which(partition == g1)) + length(which(partition == g2)) <= smax) {
    #       num.merges <- num.merges + 1
    #       merges[g1,g2] <- 1
    #     }
    #   }
    # }
    sums <- outer(sizes,sizes,FUN="+")
    merges <- sums <= smax
    diag(merges) <- 0
    num.merges <- sum(merges)
  }
  
  
  # divisions
  # for(k in 1:num.groups){
  #   sg <- sizes[k]
  #   if(smin>1){
  #     extras <- sum(unlist(lapply(1:(smin-1),FUN=function(x){choose(sg,x)})))
  #     if(sg%%2 == 0 && sg/2 < smin) extras <- extras - choose(sg,sg/2)/2
  #   } else {
  #     extras <- 0
  #   }
  #   nums.divisions[k] <- 2^(sg-1) - 1 - extras
  #   num.divisions <- num.divisions + nums.divisions[k]
  # }
  
  nums.divisions <- 2^(sizes-1) - 1 
  if(smin > 1){
    for(k in 1:num.groups){
      sg <- sizes[k]
      extras <- sum(unlist(lapply(1:(smin-1),FUN=function(x){choose(sg,x)})))
      if(sg%%2 == 0 && sg/2 < smin) extras <- extras - choose(sg,sg/2)/2
      nums.divisions[k] <- nums.divisions[k] - extras
    }
  }
  num.divisions <- sum(nums.divisions)
  
  return(list(num.merges = num.merges,
              merges = merges,
              num.divisions = num.divisions,
              nums.divisions = nums.divisions))
}


sample_new_partition_p2_restricted <- function(current.partition, mini.steps, sizes.simulated, size_neighborhood){
  
  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- size_neighborhood
  
  num.merges <- size$num.merges
  merges <- size$merges
  
  num.divisions <- size$num.divisions
  nums.divisions <- size$nums.divisions
  
  # sizes allowed
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)
  
  total <- num.merges + num.divisions 
  
  
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
    
    # reassign one of the groups and remove useless ids
    new.partition[which(new.partition == old_g2)] <- old_g1
    new.partition <- order_groupids(new.partition)
    
    nodes.changed <- which(current.partition == old_g1 | current.partition == old_g2)
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
    new.partition <- order_groupids(new.partition)
    
    nodes.changed <- which(current.partition == old_g1)
  }
  
  return(list(new.partition = new.partition,
              nodes.changed = nodes.changed))
  
}




## NEIGHBORHOOD PI 3: only swaps one node
# BUT REMOVE OUT OF SAMPLE PARTITIONS
# WARNING!!!!: for now it only works if sizes allowed are an interval ([size_min,size_max])
compute_size_neighborhood_p3_restricted <- function(partition, sizes.simulated){
  
  # find maximum size allowed
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)
  sizes <- table(partition)
  
  # find isolates, pairs and groups>2
  #isolates <- c()
  #pairs <- c()
  #others <- c()
  #abovemin_groups <- c()
  #belowmax_groups <- c()
  num.groups <- max(partition)
  num.nodes <- length(partition)
  
  # for(g in 1:num.groups){
  #   if(length(which(partition == g)) == 1) {
  #     isolates <- c(isolates,g)
  #   } else if(length(which(partition == g)) == 2) {
  #     pairs <- c(pairs,g)
  #   } else {
  #     others <- c(others,g)
  #   }
  #   # find the ones > smin and < smax (the first cannot be left, the second cannot be joined)
  #   if(length(which(partition == g)) > smin) {
  #     abovemin_groups <- c(abovemin_groups,g)
  #   }
  #   if(length(which(partition == g)) < smax) {
  #     belowmax_groups <- c(belowmax_groups,g)
  #   }
  # }
  isolates <- which(sizes == 1)
  pairs <- which(sizes == 2)
  others <- which(sizes > 2)
  abovemin_groups <- which(sizes > smin)
  belowmax_groups <- which(sizes < smax)
  
  nums.swaps <- rep(0,num.nodes)
  # done.pairs <- rep(0,length(pairs))
  # done.isolates <- rep(0,length(isolates))
  # for(i in 1:num.nodes){
  #   g <- partition[i]
  #   
  #   if(g %in% isolates){
  #     # can join any other group (<smax), except other isolates done before (because already done) and this one (that would not change anything)
  #     # since it's an isolate, it means that isolates are allowed in that context so we don't worry about smin (=1)
  #     done.isolates[isolates == g] <- 1
  #     nums.swaps[i] <- length(belowmax_groups) - sum(done.isolates) 
  #   }
  #   if(g %in% pairs && g %in% abovemin_groups){
  #     # if isolates are allowed, then it can join any other group (smax), and splitting if it's not already counted
  #     # here we don't worry about smin (=1) because 2 > smin
  #     gnotbelowmax <- !(g %in% belowmax_groups)
  #     nums.swaps[i] <- length(belowmax_groups) + gnotbelowmax - done.pairs[pairs == g]
  #     done.pairs[pairs == g] <- 1
  #   } # if size 2 is the minimum allowed, then we can never break a pair
  #   if(g %in% others && g %in% abovemin_groups){
  #     # can join any other group (<smax) or create its own isolate (if it's allowed)
  #     # we don't care here about checking that other groups are above min size, if they exist it means they are already allowed
  #     gnotbelowmax <- !(g %in% belowmax_groups)
  #     if(smin == 1) nums.swaps[i] <- length(belowmax_groups) + gnotbelowmax
  #     if(smin > 1) nums.swaps[i] <- length(belowmax_groups) + gnotbelowmax - 1
  #   } # we can never break a group that is already of minimal size
  # }
  for(i in 1:num.nodes){
    nums.swaps[i] <- (partition[i] %in% abovemin_groups) * (num.groups - length(belowmax_groups) + (partition[i] %in% belowmax_groups))
  }
  
  num.swaps <- sum(nums.swaps)
  
  return(list(isolates = isolates,
              pairs = pairs,
              others = others,
              abovemin_groups = abovemin_groups,
              belowmax_groups = belowmax_groups,
              nums.swaps = nums.swaps,
              num.swaps = num.swaps))
  
}

sample_new_partition_p3_restricted <- function(current.partition, mini.steps, sizes.simulated, size_neighborhood){
  
  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- size_neighborhood
  isolates <- size$isolates
  pairs <- size$pairs
  others <- size$others
  abovemin_groups <- size$abovemin_groups
  belowmax_groups <- size$belowmax_groups
  nums.swaps <- size$nums.swaps
  num.swaps <- size$num.swaps
  
  # allowed sizes
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)
  
  if(mini.steps == "normalized") {
    total <- num.swaps
  } 
  
  pick.1 <- sample(total,1)
  
  # decide which actor to swap
  all.groups <- 1:num.groups
  new.partition <- current.partition
  # done.pairs <- rep(0,length(pairs))
  # done.isolates <- rep(0,length(isolates))
  # 
  # for(i in 1:num.nodes){
  #   
  #   if(i == 1) {start <- 0} else {start <- sum(nums.swaps[1:i-1])}
  #   if(i == num.nodes) {end <- num.swaps} else {end <- sum(nums.swaps[1:i])}
  #   
  #   g <- current.partition[i]
  #   if(g %in% isolates) done.isolates[isolates == g] <- 1
  #   
  #   if(pick.1 > start && pick.1 <= end) {
  #     
  #     # an isolate can join any other group (<smax), except other isolates done before (because already done) and this one (that would not change anything)
  #     if(g %in% isolates){
  #       previousisolates <- isolates[done.isolates == 1]
  #       tosample <- all.groups[(all.groups %in% belowmax_groups) & !(all.groups %in% previousisolates)]
  #       if(length(tosample) == 1) newg <- tosample
  #       if(length(tosample) >= 2) newg <- sample(tosample,1)
  #       new.partition[i] <- newg
  #       new.partition <- order_groupids(new.partition)
  #     }
  #     
  #     # a member of a pair can join any other group (<smax), and splitting if it's not already counted, if isolates are allowed
  #     if(g %in% pairs && g %in% abovemin_groups){
  #       if(done.pairs[pairs == g] == 0){
  #         tosample <- all.groups[(all.groups %in% belowmax_groups) | (all.groups == g)]
  #         if(length(tosample) == 1) newg <- tosample
  #         if(length(tosample) >= 2) newg <- sample(tosample,1)
  #         if(newg == g){
  #           newg <- num.groups + 1
  #         }
  #       } else {
  #         tosample <- all.groups[(all.groups %in% belowmax_groups) & (all.groups != g)]
  #         if(length(tosample) == 1) newg <- tosample
  #         if(length(tosample) >= 2) newg <- sample(tosample,1)
  #       }
  #       new.partition[i] <- newg
  #       new.partition <- order_groupids(new.partition)
  #     }
  #     
  #     # a member of a bigger group (>smin) can join any other group (<smax) or create its own isolate (if it's allowed)
  #     if(g %in% others && g %in% abovemin_groups){
  #       if(smin == 1) {
  #         tosample <- all.groups[(all.groups %in% belowmax_groups) | (all.groups == g)]
  #         if(length(tosample) == 1) newg <- tosample
  #         if(length(tosample) >= 2) newg <- sample(tosample,1)
  #         if(newg == g){
  #           newg <- num.groups + 1
  #         }
  #       } else if(smin > 1){
  #         tosample <- all.groups[(all.groups %in% belowmax_groups) & (all.groups != g)]
  #         if(length(tosample) == 1) newg <- tosample
  #         if(length(tosample) >= 2) newg <- sample(tosample,1)
  #       } 
  #       new.partition[i] <- newg
  #       new.partition <- order_groupids(new.partition)
  #     }
  #     
  #   }
  #   
  #   if(g %in% pairs) done.pairs[pairs == g] <- 1
  # }
  
  probas_nodes <- nums.swaps / num.swaps
  picknode <- sample(1:num.nodes,1,prob=probas_nodes)
  
  topick <- unique(belowmax_groups,current.partition[picknode])
  pickgroup <- sample(topick,1)
  if(pickgroup == current.partition[picknode]) {
    new.partition[picknode] <- num.groups + 1
    nodes.changed <- which(current.partition == current.partition[picknode])
  } else {
    new.partition[picknode] <- pickgroup
    nodes.changed <- which(current.partition == current.partition[picknode] | current.partition == pickgroup)
  }
  new.partition <- order_groupids(new.partition)
  
  return(list(new.partition = new.partition,
              nodes.changed = nodes.changed))
  
}



## NEIGHBORHOOD PI 4: swaps a pair of nodes
# WARNING!!!!: for now it only works if sizes allowed are an interval ([size_min,size_max])
compute_size_neighborhood_p4_restricted <- function(partition, sizes.simulated){
  
  # check if current partition is allowed
  if(!check_sizes(partition, sizes.simulated)) {
    print(partition)
    #stop("The partition we are in is not allowed.")
  } 
  
  # find maximum size allowed
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)
  
  # find isolates, pairs and groups>2
  isolates <- c()
  pairs <- c()
  others <- c()
  abovemin_groups <- c()
  belowmax_groups <- c()
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
    # find the ones > smin+1 and < smax-1 (the first cannot be left, the second cannot be joined)
    if(length(which(partition == g)) > (smin+1)) {
      abovemin_groups <- c(abovemin_groups,g)
    }
    if(length(which(partition == g)) < (smax-1)) {
      belowmax_groups <- c(belowmax_groups,g)
    }
  }
  
  nums.swaps <- rep(0,max(partition))
  done.pairs <- rep(0,length(pairs))
  for(g in 1:num.groups){
    
    if(g %in% isolates){
      # nothing happens, there is no pair to swap here
    }
    if(g %in% pairs){
      # then the pair can move to a group <smax-1, unless it's another pair and the combination of the two pairs is already counted
      # is 2 <smax-1 we have to remove the pair from the count of possible other groups
      if(2 < smax-1) nums.swaps[g] <- length(belowmax_groups) - 1 - sum(done.pairs) + done.pairs[pairs==g]
      if(2 >= smax-1) nums.swaps[g] <- length(belowmax_groups) # then it's only isolates
      done.pairs[pairs == g] <- 1
    } 
    if(g %in% others && g %in% abovemin_groups){
      # then any pair inside the group can be separated and added to an available group that is <smax-1, unless the group is not >smin+1
      sizeg <- sum(partition == g)
      gbelowmax <- (g %in% belowmax_groups)
      if(gbelowmax) nums.swaps[g] <- (length(belowmax_groups)-1) * sizeg * (sizeg-1) / 2
      if(!gbelowmax) nums.swaps[g] <- length(belowmax_groups) * sizeg * (sizeg-1) / 2
    } # we can never break a group that is already of minimal size
  }
  
  num.swaps <- sum(nums.swaps)
  
  return(list(isolates = isolates,
              pairs = pairs,
              others = others,
              abovemin_groups = abovemin_groups,
              belowmax_groups = belowmax_groups,
              nums.swaps = nums.swaps,
              num.swaps = num.swaps))
  
}

sample_new_partition_p4_restricted <- function(current.partition, mini.steps, sizes.simulated, size_neighborhood){
  
  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- size_neighborhood
  isolates <- size$isolates
  pairs <- size$pairs
  others <- size$others
  abovemin_groups <- size$abovemin_groups
  belowmax_groups <- size$belowmax_groups
  nums.swaps <- size$nums.swaps
  num.swaps <- size$num.swaps
  
  # allowed sizes
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)
  
  if(mini.steps == "normalized") {
    total <- num.swaps
  } 
  
  pick.1 <- sample(total,1)
  
  # decide which group to pick a pair from
  all.groups <- 1:num.groups
  new.partition <- current.partition
  done.pairs <- rep(0,length(pairs))
  
  for(g in 1:num.groups){
    
    if(g == 1) {start <- 0} else {start <- sum(nums.swaps[1:g-1])}
    if(g == num.groups) {end <- num.swaps} else {end <- sum(nums.swaps[1:g])}
    
    if(pick.1 > start && pick.1 <= end) {
      
      # a pair can join other available groups (<smax-1), unless it's a pair that has already been "joined"
      if(g %in% pairs){
        members <- which(current.partition == g)
        previous.pairs <- pairs[done.pairs == 1]
        if(2 < smax-1){
          tosample <- all.groups[(all.groups %in% belowmax_groups) & (all.groups != g) & !(all.groups %in% previous.pairs)]
          if(length(tosample) == 1) newg <- tosample
          if(length(tosample) >= 2) newg <- sample(tosample,1)
        }else{
          tosample <- all.groups[(all.groups %in% belowmax_groups)]
          if(length(tosample) == 1) newg <- tosample
          if(length(tosample) >= 2) newg <- sample(tosample,1)
        }
        new.partition[members] <- newg
        new.partition <- order_groupids(new.partition)
      }
      
      # a pair in a bigger group (>smin+1) can join any other group (<smax-1)
      if(g %in% others && g %in% abovemin_groups){
        members <- which(current.partition == g)
        sizeg <- length(members)
        pairg <- sample(members,2)
        tosample <- all.groups[(all.groups %in% belowmax_groups) & (all.groups != g)]
        if(length(tosample) == 1) newg <- tosample
        if(length(tosample) >= 2) newg <- sample(tosample,1)
        new.partition[pairg] <- newg
        new.partition <- order_groupids(new.partition)
      }
      
    }
    
    if(g %in% pairs) done.pairs[pairs == g] <- 1
  }
  
  return(new.partition)
  
}


## NEIGHBORHOOD PI 5: only swaps of two pairs of nodes (careful, cannot be used alone)
compute_size_neighborhood_p5_restricted <- function(partition, sizes.simulated){
  # check if current partition is allowed
  if(!check_sizes(partition, sizes.simulated)) stop("The partition we are in is not allowed.")
  return(compute_size_neighborhood_p5(partition))
}

sample_new_partition_p5_restricted <- function(current.partition, mini.steps, sizes.simulated, size_neighborhood){
  return(sample_new_partition_p5(current.partition, mini.steps, size_neighborhood))
}


## NEIGHBORHOOD PI 6: exchanges of two groups (careful, cannot be used alone)
compute_size_neighborhood_p6_restricted <- function(partition, sizes.simulated){
  # check if current partition is allowed
  if(!check_sizes(partition, sizes.simulated)) stop("The partition we are in is not allowed.")
  return(compute_size_neighborhood_p6(partition))
}

sample_new_partition_p6_restricted <- function(current.partition, mini.steps, sizes.simulated, size_neighborhood){
  return(sample_new_partition_p6(current.partition, mini.steps, size_neighborhood))
}


## NEIGHBORHOOD PI 9: take two nodes together in partition 2, bring them together (single move)
## otherwise separate then
## careful need isolates
compute_size_neighborhood_p9_restricted <- function(partition1, partition2, sizes.simulated){
  
  # find maximum size allowed
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)
  
  num.nodes <- length(partition1)
  present1 <- !is.na(partition1)
  present2 <- !is.na(partition2)
  
  num.nodes <- length(partition1)
  num.groups <- max(partition1,na.rm=T)
  sizes <- table(partition1)
  abovemin_groups <- which(sizes > smin)
  belowmax_groups <- which(sizes < smax)
  
  affiliation1 <- as.matrix(table(data.frame(actor = 1:num.nodes, group= partition1)))
  affiliation2 <- as.matrix(table(data.frame(actor = 1:num.nodes, group= partition2)))
  adjacency1 <- affiliation1 %*% t(affiliation1)
  adjacency2 <- affiliation2 %*% t(affiliation2)
  
  nums.swaps <- matrix(0,num.nodes,num.nodes)
  for(i in 1:num.nodes) {
    for(j in 1:num.nodes){
      ad1 <- adjacency1[i,j]
      ad2 <- adjacency2[i,j]
      if(!is.na(ad1) && !is.na(ad2) && ad2 == 1) {
        if(ad1 == 0) {
          nums.swaps[i,j] <- (partition1[i] %in% abovemin_groups) & (partition1[j] %in% belowmax_groups)
        }
        if(ad1 == 1) {
          if(partition1[i] %in% abovemin_groups) nums.swaps[i,j] <- num.groups - length(belowmax_groups) + (partition1[i] %in% belowmax_groups)
        }
      }
    }
  } 
  
  return(list(abovemin_groups = abovemin_groups,
              belowmax_groups = belowmax_groups,
              adjacency1 = adjacency1,
              nums.swaps = nums.swaps,
              num.swaps = sum(nums.swaps)))
}

sample_new_partition_p9_restricted <- function(current.partition1, current.partition2, size_neighborhood, sizes.simulated){
  
  # find maximum size allowed
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)
  
  num.nodes <- length(current.partition1)
  num.groups <- max(current.partition1,na.rm=T)
  abovemin_groups <- size_neighborhood$abovemin_groups
  belowmax_groups <- size_neighborhood$belowmax_groups
  adjacency1 <- size_neighborhood$adjacency1
  nums.swaps <- size_neighborhood$nums.swaps
  num.swaps <- size_neighborhood$num.swaps
  
  present1 <- !(is.na(current.partition1))
  new.partition1 <- current.partition1
  
  probas_pairs <- as.vector(nums.swaps) / num.swaps
  pick <- sample(1:(num.nodes*num.nodes),1,prob=probas_pairs)
  
  i <- (pick-1) %% num.nodes + 1
  j <- ((pick-1) %/% num.nodes) + 1
  
  if(adjacency1[i,j] == 0) {
    new.partition1[i] <- current.partition1[j]
    nodes.changed <- which(current.partition1 == current.partition1[i] | current.partition1 == current.partition1[j])
  }
  if(adjacency1[i,j] == 1) {
    topick <- unique(belowmax_groups,current.partition1[i])
    pickgroup <- sample(topick,1)
    if(pickgroup == current.partition1[i]) {
      new.partition1[i] <- num.groups + 1
      nodes.changed <- which(current.partition1 == current.partition1[i])
    } else {
      new.partition1[i] <- pickgroup
      nodes.changed <- which(current.partition1 == current.partition1[i] | current.partition1 == pickgroup)
    }
  }
  new.partition1[present1] <- order_groupids(new.partition1[present1])
  
  return(list(new.partition = new.partition1,
              nodes.changed = nodes.changed))
  
}