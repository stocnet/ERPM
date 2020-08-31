######################################################################
## Simulation and estimation of Exponential Random Partition Models ## 
## Functions used to run the phase 3 of the estimation algorithm    ##
## Author: Marion Hoffman                                           ##
######################################################################


# Phase 3 for a single partition
run_phase3_single <- function(partition,
                       estimates.phase2, 
                       z.obs, 
                       nodes, 
                       effects, 
                       objects, 
                       burnin, 
                       thining,
                       a.scaling,
                       mini.steps, 
                       length.p3, 
                       neighborhood,
                       sizes.allowed,
                       sizes.simulated,
                       fixed.estimates,
                       parallel = F,
                       cpus = 1) {
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)
  
  # find a good starting point
  #first.partition <- find_startingpoint_single(nodes,sizes.allowed)
  first.partition <- partition
  
  # simulate a large sample with the estimates found in phase 2 
  if(parallel){
    
    sfExport("startingestimates", "first.partition", "nodes", "effects", "objects", "burnin", "thining", "length.p3", "cpus", "mini.steps", "neighborhood", "sizes.allowed", "sizes.simulated")
    res <- sfLapply(1:cpus, fun = function(k) {
      set.seed(k)
      subres <- draw_Metropolis_single(estimates.phase2, first.partition, nodes, effects, objects, burnin, thining, ceiling(length.p3/cpus), mini.steps, neighborhood, sizes.allowed, sizes.simulated)
      return(subres)
    }
    )
    all.z <- c()
    for(k in 1:cpus) all.z <- rbind(all.z,res[[k]]$draws)
    length.p3 <- cpus * ceiling(length.p3/cpus)
    results.phase3 <- list("draws" = all.z, "last.partition" = res[[cpus]]$last.partition, "all.partitions" = NULL) 
    
  }else{
    
    results.phase3 <- draw_Metropolis_single(estimates.phase2, first.partition, nodes, effects, objects, burnin, thining, length.p3, mini.steps, neighborhood, sizes.allowed, sizes.simulated)
  }
  z.phase3 <- results.phase3$draws
  
  # calculate covariance and scaling
  inverted_matrices <- calculate_inverted_covariance_and_scaling(estimates.phase2, 
                                                        z.obs, 
                                                        nodes, 
                                                        effects, 
                                                        objects,
                                                        a.scaling,
                                                        length.phase = length.p3, 
                                                        z.phase = z.phase3,
                                                        fixed.estimates)
  
  # hack for size constraints
  if(!is.null(sizes.allowed)){
    length.p3 <- nrow(z.phase3)
    print("new length of phase 3")
    print(length.p3)
  }
  
  # calculations of phase 3: mean, sd, se, conv ratios
  res.phase3 <- phase3(estimates.phase2, z.phase3, z.obs, effects, length.p3)
  
  return(list("means" = res.phase3$finalmean, 
              "standard.deviations" = res.phase3$finalsd, 
              "standard.errors" = res.phase3$finalse, 
              "convergence.ratios" = res.phase3$finalconvratios,
              "inv.zcov" = inverted_matrices$inv.zcov,
              "inv.scaling" = inverted_matrices$inv.scaling))
  
}


# Phase 3 for multiple partitions
run_phase3_multiple <- function(partitions,
                              estimates.phase2, 
                              z.obs, 
                              presence.tables, 
                              nodes, 
                              effects, 
                              objects, 
                              burnin, 
                              thining,
                              a.scaling,
                              mini.steps, 
                              length.p3, 
                              neighborhood,
                              sizes.allowed,
                              sizes.simulated,
                              fixed.estimates,
                              parallel = F,
                              cpus = 1) {
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)
  
  # find a good starting point
  #first.partitions <- find_startingpoint_multiple(presence.tables,nodes,sizes.allowed)
  first.partitions <- partitions
  
  # simulate a large sample with the estimates found in phase 2 
  if(parallel){
    
    sfExport("startingestimates", "first.partitions", "presence.tables", "nodes", "effects", "objects", "burnin", "thining", "length.p3", "cpus", "mini.steps", "neighborhood", "sizes.allowed", "sizes.simulated")
    res <- sfLapply(1:cpus, fun = function(k) {
      set.seed(k)
      subres <- draw_Metropolis_multiple(estimates.phase2, first.partitions, presence.tables, nodes, effects, objects, burnin, thining, ceiling(length.p3/cpus), mini.steps, neighborhood, sizes.allowed, sizes.simulated)
      return(subres)
    }
    )
    all.z <- c()
    for(k in 1:cpus) all.z <- rbind(all.z,res[[k]]$draws)
    length.p3 <- cpus * ceiling(length.p3/cpus)
    results.phase3 <- list("draws" = all.z, "last.partitions" = res[[cpus]]$last.partitions, "all.partitions" = NULL) 
  
  }else{
    
    results.phase3 <- draw_Metropolis_multiple(estimates.phase2, first.partitions, presence.tables, nodes, effects, objects, burnin, thining, length.p3, mini.steps, neighborhood, sizes.allowed, sizes.simulated)
  
  }
  z.phase3 <- results.phase3$draws
  
  # calculate covariance and scaling
  inverted_matrices <- calculate_inverted_covariance_and_scaling(estimates.phase2, 
                                                                 z.obs, 
                                                                 nodes, 
                                                                 effects, 
                                                                 objects,
                                                                 a.scaling,
                                                                 length.phase = length.p3, 
                                                                 z.phase = z.phase3,
                                                                 fixed.estimates)
  
  # hack for size constraints
  if(!is.null(sizes.allowed)){
    length.p3 <- nrow(z.phase3)
    print("new length of phase 3")
    print(length.p3)
  }
  
  # calculations of phase 3: mean, sd, se, conv ratios
  res.phase3 <- phase3(estimates.phase2, z.phase3, z.obs, effects, length.p3)
  
  return(list("means" = res.phase3$finalmean, 
              "standard.deviations" = res.phase3$finalsd, 
              "standard.errors" = res.phase3$finalse, 
              "convergence.ratios" = res.phase3$finalconvratios,
              "inv.zcov" = inverted_matrices$inv.zcov,
              "inv.scaling" = inverted_matrices$inv.scaling))
  
}


# function for the core calculations of phase 3: estimate means, s.d., s.e., and convergence ratios
phase3 <- function(estimates.phase2,
                   z.phase3,
                   z.obs,
                   effects,
                   length.p3){
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)
  
  # mean statistics
  finalmean <- rep(0,num.effects)
  for(e in 1:num.effects) {
    finalmean[e] <- Reduce('+',z.phase3[,e]) / length.p3
  }
  
  # standard deviations
  finalsd <- rep(0, num.effects)
  for(e in 1:num.effects) {
    statse <- rep(0,length.p3)
    for(i in 1:length.p3) {
      statse[i] <- z.phase3[i,e]
    }
    finalsd[e] <- sd(statse)
  }
  finalse <- finalsd / sqrt(length.p3)
  
  # convergence ratios
  finalconvratios <- (finalmean - z.obs) / finalsd
  
  print("Estimated statistics after phase 3")
  print(finalmean)
  print("Estimates after phase 3")
  print(estimates.phase2)
  
  return(list("finalmean" = finalmean,
              "finalsd" = finalsd,
              "finalse" = finalse,
              "finalconvratios" = finalconvratios))
}
