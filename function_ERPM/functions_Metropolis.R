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
  current.z <- computeStatistics_multiple(current.partitions, presence.tables, nodes, effects, objects)
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
  
  while(!end.walk){
    
    new.step <- draw_step_multiple(theta,
                                 current.partitions, 
                                 current.logit,
                                 current.z,
                                 presence.tables,  
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
    old.partitions <- current.partitions
    if(change.made) {
      current.partitions <- new.step$new.partitions
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
  
  move <- sample(c("swap","mergediv","single"),1,prob=neighborhood)
  
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
  
  move <- sample(c("swap","mergediv","single"),1,prob=neighborhood)
  
  ## IF MODE IS "swap" = swap two actors, "mergediv" = merge 2 groups or split a group in two, "single" = move one actor
  if(move == "swap" && is.null(sizes.allowed)) {
    current.size <- compute_size_neighborhood_p1(current.partitions[nodes.rand.o,rand.o])
    if(current.size$num.swaps > 0){
      new.p <- sample_new_partition_p1(current.partitions[nodes.rand.o,rand.o], mini.steps, current.size)
    } else {
      new.p <- current.partitions[nodes.rand.o,rand.o]
    }
  } else if(move == "mergediv" && is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p2(current.partitions[nodes.rand.o,rand.o])
    new.p <- sample_new_partition_p2(current.partitions[nodes.rand.o,rand.o], mini.steps, current.size)
  } else if(move == "single" && is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p3(current.partitions[nodes.rand.o,rand.o])
    new.p <- sample_new_partition_p3(current.partitions[nodes.rand.o,rand.o], mini.steps, current.size)
  }
  # Alternative if sizes are restricted
  if(move == "swap" && !is.null(sizes.allowed)) {
    current.size <- compute_size_neighborhood_p1_restricted(current.partitions[nodes.rand.o,rand.o], sizes.simulated)
    if(current.size$num.swaps > 0){
      new.p <- sample_new_partition_p1_restricted(current.partitions[nodes.rand.o,rand.o], mini.steps, sizes.simulated, current.size)
    } else {
      new.p <- current.partitions[nodes.rand.o,rand.o]
    }
  } else if(move == "mergediv" && !is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p2_restricted(current.partitions[nodes.rand.o,rand.o], sizes.simulated)
    new.p <- sample_new_partition_p2_restricted(current.partitions[nodes.rand.o,rand.o], mini.steps, sizes.simulated, current.size)
  } else if(move == "single" && !is.null(sizes.allowed)){
    current.size <- compute_size_neighborhood_p3_restricted(current.partitions[nodes.rand.o,rand.o], sizes.simulated)
    new.p <- sample_new_partition_p3_restricted(current.partitions[nodes.rand.o,rand.o], mini.steps, sizes.simulated, current.size)
  } 
  
  new.partitions[,rand.o] <- rep(NA,num.nodes)
  new.partitions[nodes.rand.o,rand.o] <- new.p
  
  # IF NO CHANGE, same statistics
  if(all(current.partitions == new.partitions, na.rm=T)) {
    new.z <- current.z
  }
  
  # IF CHANGE: compute new statistics
  if(!all(current.partitions == new.partitions, na.rm=T)) {
    
    # store new statistics
    new.z <- rep(0,length(current.z))
    
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
        old.z.temp <- computeStatistics(current.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp)
        new.z.temp <- computeStatistics(new.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp)
        new.z[e] <- current.z[e] - old.z.temp + new.z.temp
      
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
        nodes.temp$var.temp <- atts[,rand.o]
       
        old.z.temp <- computeStatistics(current.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp)
        new.z.temp <- computeStatistics(new.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp)
        new.z[e] <- current.z[e] - old.z.temp + new.z.temp  
        
      # INERTIA_1 EFFECT: need to recalculate the partition before and after
      } else if(effects$names[e] == "inertia_1"){
        
        # first partition
        if(rand.o > 1) {
          effects.temp <- list(names="tie",objects="net.temp")
          aff.temp <- as.matrix(table(data.frame(actor = 1:num.nodes, group= new.partitions[,rand.o-1])))
          adj.temp <- aff.temp %*% t(aff.temp)
          diag(adj.temp) <- 0
          objects.temp <- list(list(name="net.temp", object=adj.temp[nodes.rand.o,nodes.rand.o]))
          nodes.temp <- nodes[nodes.rand.o,]
          old.z.temp <- computeStatistics(current.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp)
          new.z.temp <- computeStatistics(new.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp)
          new.z[e] <- current.z[e] - old.z.temp + new.z.temp
        }
        
        # second partition
        if(rand.o < num.obs){
          nodes.rand.o2 <- as.logical(presence.tables[,rand.o+1])
          effects.temp <- list(names="tie",objects="net.temp")
          aff.temp <- as.matrix(table(data.frame(actor = 1:num.nodes, group= new.partitions[,rand.o])))
          adj.temp <- aff.temp %*% t(aff.temp)
          diag(adj.temp) <- 0
          objects.temp <- list(list(name="net.temp", object=adj.temp[nodes.rand.o2,nodes.rand.o2]))
          nodes.temp <- nodes[nodes.rand.o2,]
          old.z.temp <- computeStatistics(current.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp)
          new.z.temp <- computeStatistics(new.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp)
          new.z[e] <- current.z[e] - old.z.temp + new.z.temp
        }
        
      # INERTIA_TOTAL: need to recalculate everything
      } else if(effects$names[e] == "inertia_total"){

        effects.temp <- list(names="inertia_total",objects="partitions",objects2="")
        objects.temp <- list()
        new.z.temp <- computeStatistics_multiple(new.partitions, presence.tables, nodes, effects.temp, objects.temp)
        new.z[e] <- new.z.temp

      # INERTIA_TOTAL_X_DIFF: need to recalculate everything
      } else if(effects$names[e] == "inertia_total_X_diff"){
        
        effects.temp <- list(names="inertia_total_X_diff",objects="partitions",objects2=effects$objects2[e])
        objects.temp <- list()
        new.z.temp <- computeStatistics_multiple(new.partitions, presence.tables, nodes, effects.temp, objects.temp)
        new.z[e] <- new.z.temp
        
      # ANY OTHER EFFECT: the effect also exists in the single partition version
      } else {
        
        effects.temp <- list(names=effects$names[e],objects=effects$objects[e])
        objects.temp <- objects
        nodes.temp <- nodes[nodes.rand.o,]
        
        old.z.temp <- computeStatistics(current.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp)
        new.z.temp <- computeStatistics(new.partitions[nodes.rand.o,rand.o], nodes.temp, effects.temp, objects.temp)
        new.z[e] <- current.z[e] - old.z.temp + new.z.temp  
        
      }
    }
  } 
  
  
  
  new.logit <- theta * new.z
  
  
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
    }

    hastings.ratio <- (exp(sum(new.logit) - sum(current.logit))) * neighborhoods.ratios
    
  } else if(mini.steps == "selfloops"){
    
    hastings.ratio <- exp(sum(new.logit) - sum(current.logit))
    
  }
  
  return(list("hastings.ratio" = hastings.ratio,
              "new.partitions" = new.partitions,
              "new.z" = new.z,
              "new.logit" = new.logit))
}



## NEIGHBORHOOD PI 1: only swaps of two nodes (careful, cannot be used alone)
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
  
  nums.swaps <- matrix(0,num.nodes,num.nodes)
  for(i in 1:(num.nodes-1)){
    for(j in (i+1):num.nodes){
      g_i <- partition[i]
      g_j <- partition[j]
      
      allowed <- T
      if(g_i == g_j) allowed <- F # two members of the same group cannot swap, it wouldn't change anything
      if(g_i %in% isolates && g_j %in% isolates){
        allowed <- F # two isolates cannot swap, it wouldn't change anything
      }
      
      if(allowed) nums.swaps[i,j] <- 1
    }
  }
  
  num.swaps <- sum(nums.swaps)
  
  return(list(isolates = isolates,
              pairs = pairs,
              others = others,
              nums.swaps = nums.swaps,
              num.swaps = num.swaps))
  
}

sample_new_partition_p1 <- function(current.partition, mini.steps, size_neighborhood){
  
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
  
  # decide which actors to swap
  new.partition <- current.partition
  swap <- which(nums.swaps == 1)[pick.1]
  if(swap %% num.nodes == 0){
    i <- num.nodes
    j <- swap %/% num.nodes
  } else {
    i <- swap %% num.nodes
    j <- swap %/% num.nodes + 1
  }
  new.partition[i] <- current.partition[j]
  new.partition[j] <- current.partition[i]
  new.partition <- order_groupids(new.partition)
  
  if(all(new.partition == current.partition)){
    print("test")
  }
  
  return(new.partition)
  
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
    new.partition <- order_groupids(new.partition)
  }
  
  return(new.partition)
  
}




## NEIGHBORHOOD PI 3: only swaps of one node, careful, buggy (sometimes returns not allwoed group?)
# BUT REMOVE OUT OF SAMPLE PARTITIONS
# WARNING!!!!: for now it only works if sizes allowed are an interval ([size_min,size_max])
compute_size_neighborhood_p3_restricted <- function(partition, sizes.simulated){
  
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
    # find the ones > smin and < smax (the first cannot be left, the second cannot be joined)
    if(length(which(partition == g)) > smin) {
      abovemin_groups <- c(abovemin_groups,g)
    }
    if(length(which(partition == g)) < smax) {
      belowmax_groups <- c(belowmax_groups,g)
    }
  }
  
  nums.swaps <- rep(0,num.nodes)
  done.pairs <- rep(0,length(pairs))
  done.isolates <- rep(0,length(isolates))
  for(i in 1:num.nodes){
    g <- partition[i]
    
    if(g %in% isolates){
      # can join any other group (<smax), except other isolates done before (because already done) and this one (that would not change anything)
      # since it's an isolate, it means that isolates are allowed in that context so we don't worry about smin (=1)
      done.isolates[isolates == g] <- 1
      nums.swaps[i] <- length(belowmax_groups) - sum(done.isolates) 
    }
    if(g %in% pairs && g %in% abovemin_groups){
      # if isolates are allowed, then it can join any other group (smax), and splitting if it's not already counted
      # here we don't worry about smin (=1) because 2 > smin
      gnotbelowmax <- !(g %in% belowmax_groups)
      nums.swaps[i] <- length(belowmax_groups) + gnotbelowmax - done.pairs[pairs == g]
      done.pairs[pairs == g] <- 1
    } # if size 2 is the minimum allowed, then we can never break a pair
    if(g %in% others && g %in% abovemin_groups){
      # can join any other group (<smax) or create its own isolate (if it's allowed)
      # we don't care here about checking that other groups are above min size, if they exist it means they are already allowed
      gnotbelowmax <- !(g %in% belowmax_groups)
      if(smin == 1) nums.swaps[i] <- length(belowmax_groups) + gnotbelowmax
      if(smin > 1) nums.swaps[i] <- length(belowmax_groups) + gnotbelowmax - 1
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
  done.pairs <- rep(0,length(pairs))
  done.isolates <- rep(0,length(isolates))
  
  for(i in 1:num.nodes){
    
    if(i == 1) {start <- 0} else {start <- sum(nums.swaps[1:i-1])}
    if(i == num.nodes) {end <- num.swaps} else {end <- sum(nums.swaps[1:i])}
    
    g <- current.partition[i]
    if(g %in% isolates) done.isolates[isolates == g] <- 1
    
    if(pick.1 > start && pick.1 <= end) {
      
      # an isolate can join any other group (<smax), except other isolates done before (because already done) and this one (that would not change anything)
      if(g %in% isolates){
        previousisolates <- isolates[done.isolates == 1]
        tosample <- all.groups[(all.groups %in% belowmax_groups) & !(all.groups %in% previousisolates)]
        if(length(tosample) == 1) newg <- tosample
        if(length(tosample) >= 2) newg <- sample(tosample,1)
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }
      
      # a member of a pair can join any other group (<smax), and splitting if it's not already counted, if isolates are allowed
      if(g %in% pairs && g %in% abovemin_groups){
        if(done.pairs[pairs == g] == 0){
          tosample <- all.groups[(all.groups %in% belowmax_groups) | (all.groups == g)]
          if(length(tosample) == 1) newg <- tosample
          if(length(tosample) >= 2) newg <- sample(tosample,1)
          if(newg == g){
            newg <- num.groups + 1
          }
        } else {
          tosample <- all.groups[(all.groups %in% belowmax_groups) & (all.groups != g)]
          if(length(tosample) == 1) newg <- tosample
          if(length(tosample) >= 2) newg <- sample(tosample,1)
        }
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }
      
      # a member of a bigger group (>smin) can join any other group (<smax) or create its own isolate (if it's allowed)
      if(g %in% others && g %in% abovemin_groups){
        if(smin == 1) {
          tosample <- all.groups[(all.groups %in% belowmax_groups) | (all.groups == g)]
          if(length(tosample) == 1) newg <- tosample
          if(length(tosample) >= 2) newg <- sample(tosample,1)
          if(newg == g){
            newg <- num.groups + 1
          }
        } else if(smin > 1){
          tosample <- all.groups[(all.groups %in% belowmax_groups) & (all.groups != g)]
          if(length(tosample) == 1) newg <- tosample
          if(length(tosample) >= 2) newg <- sample(tosample,1)
        } 
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }
    
    }
    
    if(g %in% pairs) done.pairs[pairs == g] <- 1
  }
  
  return(new.partition)
  
}

