######################################################################
## Simulation and estimation of Exponential Random Partition Models ## 
## Functions used to find appropriate burnin and thining in phase 1 ##
## Author: Marion Hoffman                                           ##
######################################################################

# SINGLE PARTITION PROCEDURE
simulate_burninthining_single <- function(partition, # observed partition
                               theta, # initial model parameters
                               nodes, # nodeset (data frame)
                               effects, # effects/sufficient statistics (list with a vector "names", and a vector "objects")
                               objects, # objects used for statistics calculation (list with a vector "name", and a vector "object")
                               num.steps, # number of samples wanted in phase 1
                               neighborhood, # way of choosing partitions: probability vector (proba actors swap, proba merge/division, proba single actor move)
                               sizes.allowed, # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                               sizes.simulated, # vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max) 
                               max.thining)
{
  num.effects <- length(effects$names)
  
  chain <- draw_Metropolis_single(theta, partition, nodes, effects, objects, 1, 1, num.steps, neighborhood, sizes.allowed, sizes.simulated, return.all.partitions = F)                              
  
  # first look at autocorrelations find ok thining (either max thining reached or autocor < 0.4)
  ok <- F
  current.draws <- chain$draws
  current.last <- chain$last.partition
  current.thining <- 1
  
  allautocors <- c()
  
  while(!ok){
    
    # we check the auto for the current burnin
    draws <- current.draws[seq(1, nrow(current.draws), current.thining),]
    autocors <- rep(0,num.effects)
    for(e in 1:num.effects){
      autocors[e] <- cor(draws[1:(num.steps-1),e],draws[2:num.steps,e])
    }
    if(current.thining %% 50 == 0){
      print("thining")
      print(current.thining)
      print("autocorrelations")
      print(autocors)
    }
    
    allautocors <- rbind(allautocors,autocors)
    
    #if(max(autocors,na.rm=T) < 0.4) {
    if(current.thining >= max.thining ) {
      
      # we can stop
      ok <- T
      
    } else{
      
      # we continue the chain
      new.chain <- draw_Metropolis_single(theta, current.last, nodes, effects, objects, 1, 1, num.steps, neighborhood, sizes.allowed, sizes.simulated, return.all.partitions = F)
      current.draws <- rbind(current.draws,new.chain$draws)
      current.last <- new.chain$last.partition
      
      current.thining <- current.thining + 1
    }
    
    
  }
  
  # now we check the evolution of the mean for choosing the burnin
  allmeans <- matrix(0,nrow(draws),num.effects)
  allmeans[1,] <- draws[1,] 
  for(d in 2:nrow(draws)) {
    allmeans[d,] <- colMeans(draws[1:d,]) 
  }
  
  # we smoothen out the curves (for plots), can be useful when the sample is small
  smoothedautocors <- allautocors
  smoothedmeans <- allmeans
  for(eff in 1:num.effects) {
    lo <- loess(y ~ x, data.frame(x=1:max.thining, y=allautocors[,eff]))
    smoothedautocors[!is.na(smoothedautocors[,eff]),eff] <- lo$fitted
    lo <- loess(y ~ x, data.frame(x=1:nrow(draws), y=allmeans[,eff]))
    smoothedmeans[!is.na(smoothedmeans[,eff]),eff] <- lo$fitted
  }
  
  
  return(list("draws" = draws,
              "autocorrelations" = allautocors,
              "moving.means" = allmeans,
              "autocorrelations.smoothed" = smoothedautocors,
              "moving.means.smoothed" = smoothedmeans))
  
}


# MULTIPLE PARTITION PROCEDURE
simulate_burninthining_multiple <- function(partitions, # observed partitions
                                        presence.tables, # to indicate which nodes were present when
                                        theta, # initial model parameters
                                        nodes, # nodeset (data frame)
                                        effects, # effects/sufficient statistics (list with a vector "names", and a vector "objects")
                                        objects, # objects used for statistics calculation (list with a vector "name", and a vector "object")
                                        num.steps, # number of samples wanted in phase 1
                                        neighborhood, # way of choosing partitions: probability vector (proba actors swap, proba merge/division, proba single actor move)
                                        sizes.allowed, # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                                        sizes.simulated, # vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max) 
                                        max.thining)
{
  num.effects <- length(effects$names)
  
  chain <- draw_Metropolis_multiple(theta, partitions, presence.tables, nodes, effects, objects, 1, 1, num.steps, neighborhood, sizes.allowed, sizes.simulated, return.all.partitions = F)                              
  
  # first look at autocorrelations find ok thining (either max thining reached or autocor < 0.4)
  ok <- F
  current.draws <- chain$draws
  current.last <- chain$last.partitions
  current.thining <- 1
  
  allautocors <- c()
  
  while(!ok){
    
    # we check the auto for the current burnin
    draws <- current.draws[seq(1, nrow(current.draws), current.thining),]
    autocors <- rep(0,num.effects)
    for(e in 1:num.effects){
      autocors[e] <- cor(draws[1:(num.steps-1),e],draws[2:num.steps,e])
    }
    if(current.thining %% 50 == 0){
      print("thining")
      print(current.thining)
      print("autocorrelations")
      print(autocors)
    }
    
    allautocors <- rbind(allautocors,autocors)
    
    #if(max(autocors,na.rm=T) < 0.4) {
    if(current.thining >= max.thining ) {
      
      # we can stop
      ok <- T
      
    } else{
      
      # we continue the chain
      new.chain <- draw_Metropolis_multiple(theta, current.last, presence.tables, nodes, effects, objects, 1, 1, num.steps, neighborhood, sizes.allowed, sizes.simulated, return.all.partitions = F)
      current.draws <- rbind(current.draws,new.chain$draws)
      current.last <- new.chain$last.partitions
      
      current.thining <- current.thining + 1
    }
    
    
  }
  
  # now we check the evolution of the mean for choosing the burnin
  allmeans <- matrix(0,nrow(draws),num.effects)
  allmeans[1,] <- draws[1,] 
  for(d in 2:nrow(draws)) {
    allmeans[d,] <- colMeans(draws[1:d,]) 
  }
  
  # we smoothen out the curves (for plots), can be useful when the sample is small
  smoothedautocors <- allautocors
  smoothedmeans <- allmeans
  for(eff in 1:num.effects) {
    lo <- loess(y ~ x, data.frame(x=1:max.thining, y=allautocors[,eff]))
    smoothedautocors[!is.na(smoothedautocors[,eff]),eff] <- lo$fitted
    lo <- loess(y ~ x, data.frame(x=1:nrow(draws), y=allmeans[,eff]))
    smoothedmeans[!is.na(smoothedmeans[,eff]),eff] <- lo$fitted
  }
  
  
  return(list("draws"= draws,
              "autocorrelations" = allautocors,
              "moving.means" = allmeans,
              "autocorrelations.smoothed" = smoothedautocors,
              "moving.means.smoothed" = smoothedmeans))
  
}

# GRID SEARCH FOR SINGLE PARTITION
gridsearch_burninthining_single <- function(partition, # observed partition
                                      theta, # initial model parameters
                                      nodes, # nodeset (data frame)
                                      effects, # effects/sufficient statistics (list with a vector "names", and a vector "objects")
                                      objects, # objects used for statistics calculation (list with a vector "name", and a vector "object")
                                      num.steps, # number of samples wanted in phase 1
                                      neighborhoods, # list of probability vectors (proba actors swap, proba merge/division, proba single actor move)
                                      sizes.allowed, # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                                      sizes.simulated, # vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max) 
                                      max.thining, # where to stop adding thining
                                      parallel = F, # to run different neighborhoods in parallel
                                      cpus = 1) {
  
  # parallel procedure
  if(parallel){
    
    #distribute neighborhoods in each cpu
    n <- ceiling(length(neighborhoods) / cpus)
    subindexes <- list()
    for(c in 1:cpus){
      start <- (c-1)*n + 1
      end <- c*n
      if(c == cpus) end <- length(neighborhoods)
      subindexes[[c]] <- start:end
    }
    
    sfExport("partition", "theta", "nodes", "effects", "objects", "num.steps", "neighborhoods", "sizes.allowed", "sizes.simulated", "max.thining", "subindexes")
    res <- sfLapply(1:cpus, fun = function(k) {
      subres <- list()
      for(i in 1:length(subindexes[[k]])){
        index <- subindexes[[k]][i]
        subneighborhood <- neighborhoods[[index]]
        subres[[i]] <- simulate_burninthining_single(partition, theta, nodes, effects, objects, num.steps, subneighborhood, sizes.allowed, sizes.simulated, max.thining)
      }
      return(subres)
    }
    )
    
    # append all results
    allsimulations <- list()
    for(c in 1:cpus) allsimulations <- append(allsimulations, res[[c]])
  
  } else {
    
    # just go through all neighborhoods one by one
    allsimulations <- list()
    for(i in 1:length(neighborhoods)){
      allsimulations[[i]] <- simulate_burninthining_single(partition, theta, nodes, effects, objects, num.steps, neighborhoods[[i]], sizes.allowed, sizes.simulated, max.thining)
    }
    
  }
  
  # collect properties for all neighborhoods
  p1 <- c()
  p2 <- c()
  p3 <- c()
  reach0.4 <- c()
  firstthining0.4 <- c()
  reach0.5 <- c()
  firstthining0.5 <- c()
  finalmaxautocor <- c()
  for(i in 1:length(neighborhoods)){
    p1 <- c(p1, neighborhoods[[i]][1])
    p2 <- c(p2, neighborhoods[[i]][2])
    p3 <- c(p3, neighborhoods[[i]][3])
    
    allmaxautocors <- apply(allsimulations[[i]]$autocorrelations.smoothed, 1, FUN=max)
    finalmaxautocor <- c(finalmaxautocor, allmaxautocors[max.thining])
    
    reach0.4 <- c(reach0.4, sum(allmaxautocors <= 0.4, na.rm=T))
    reach0.5 <- c(reach0.5, sum(allmaxautocors <= 0.5, na.rm=T))
    
    if(length(which(allmaxautocors <= 0.4)) > 0) {
      firstthining0.4 <- c(firstthining0.4, min(which(allmaxautocors <= 0.4)))
    } else {
      firstthining0.4 <- c(firstthining0.4, -1)
    }
    if(length(which(allmaxautocors <= 0.5)) > 0) {
      firstthining0.5 <- c(firstthining0.5, min(which(allmaxautocors <= 0.5)))
    } else {
      firstthining0.5 <- c(firstthining0.5, -1)
    }
  }
  results.search <- data.frame(p_swap = p1,
                               p_mergediv = p2,
                               p_single = p3,
                               final_max_autocorr = finalmaxautocor,
                               reach_autocorr0.4 = reach0.4,
                               first_thining_autocorr0.4 = firstthining0.4,
                               reach_autocorr0.5 = reach0.5,
                               first_thining_autocorr0.5 = firstthining0.5)
  
  return(list("results.search" = results.search,
              "all.simulations" = allsimulations))
  
}


# GRID SEARCH FOR MULTIPLE PARTITIONS
gridsearch_burninthining_multiple <- function(partitions, # observed partitions
                                            presence.tables, # presence of nodes
                                            theta, # initial model parameters
                                            nodes, # nodeset (data frame)
                                            effects, # effects/sufficient statistics (list with a vector "names", and a vector "objects")
                                            objects, # objects used for statistics calculation (list with a vector "name", and a vector "object")
                                            num.steps, # number of samples wanted in phase 1
                                            neighborhoods, # list of probability vectors (proba actors swap, proba merge/division, proba single actor move)
                                            sizes.allowed, # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                                            sizes.simulated, # vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max) 
                                            max.thining, # where to stop adding thining
                                            parallel = F, # to run different neighborhoods in parallel
                                            cpus = 1) {
  
  # parallel procedure
  if(parallel){
    
    #distribute neighborhoods in each cpu
    n <- ceiling(length(neighborhoods) / cpus)
    subindexes <- list()
    for(c in 1:cpus){
      start <- (c-1)*n + 1
      end <- c*n
      if(c == cpus) end <- length(neighborhoods)
      subindexes[[c]] <- start:end
    }
    
    sfExport("partitions", "presence.tables", "theta", "nodes", "effects", "objects", "num.steps", "neighborhoods", "sizes.allowed", "sizes.simulated", "max.thining", "subindexes")
    res <- sfLapply(1:cpus, fun = function(k) {
      subres <- list()
      for(i in 1:length(subindexes[[k]])){
        index <- subindexes[[k]][i]
        subneighborhood <- neighborhoods[[index]]
        subres[[i]] <- simulate_burninthining_multiple(partitions, presence.tables, theta, nodes, effects, objects, num.steps, subneighborhood, sizes.allowed, sizes.simulated, max.thining)
      }
      return(subres)
    }
    )
    
    # append all results
    allsimulations <- list()
    for(c in 1:cpus) allsimulations <- append(allsimulations, res[[c]])
    
  } else {
    
    # just go through all neighborhoods one by one
    allsimulations <- list()
    for(i in 1:length(neighborhoods)){
      allsimulations[[i]] <- simulate_burninthining_multiple(partitions, presence.tables, theta, nodes, effects, objects, num.steps, neighborhoods[[i]], sizes.allowed, sizes.simulated, max.thining)
    }
    
  }
  
  # collect properties for all neighborhoods
  p1 <- c()
  p2 <- c()
  p3 <- c()
  reach0.4 <- c()
  firstthining0.4 <- c()
  reach0.5 <- c()
  firstthining0.5 <- c()
  finalmaxautocor <- c()
  for(i in 1:length(neighborhoods)){
    p1 <- c(p1, neighborhoods[[i]][1])
    p2 <- c(p2, neighborhoods[[i]][2])
    p3 <- c(p3, neighborhoods[[i]][3])
    
    allmaxautocors <- apply(allsimulations[[i]]$autocorrelations.smoothed, 1, FUN=max)
    finalmaxautocor <- c(finalmaxautocor, allmaxautocors[max.thining])
    
    reach0.4 <- c(reach0.4, sum(allmaxautocors <= 0.4, na.rm=T))
    reach0.5 <- c(reach0.5, sum(allmaxautocors <= 0.5, na.rm=T))
    
    if(length(which(allmaxautocors <= 0.4)) > 0) {
      firstthining0.4 <- c(firstthining0.4, min(which(allmaxautocors <= 0.4)))
    } else {
      firstthining0.4 <- c(firstthining0.4, -1)
    }
    if(length(which(allmaxautocors <= 0.5)) > 0) {
      firstthining0.5 <- c(firstthining0.5, min(which(allmaxautocors <= 0.5)))
    } else {
      firstthining0.5 <- c(firstthining0.5, -1)
    }
  }
  results.search <- data.frame(p_swap = p1,
                               p_mergediv = p2,
                               p_single = p3,
                               final_max_autocorr = finalmaxautocor,
                               reach_autocorr0.4 = reach0.4,
                               first_thining_autocorr0.4 = firstthining0.4,
                               reach_autocorr0.5 = reach0.5,
                               first_thining_autocorr0.5 = firstthining0.5)
  
  return(list("results.search" = results.search,
              "all.simulations" = allsimulations))
  
}