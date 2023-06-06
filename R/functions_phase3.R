######################################################################
## Simulation and estimation of Exponential Random Partition Models ## 
## Functions used to run the phase 3 of the estimation algorithm    ##
## Author: Marion Hoffman                                           ##
######################################################################


#' Phase 3 wrapper for single observation
#'
#' @param partition XXX
#' @param estimates.phase2 XXX
#' @param z.obs XXX
#' @param nodes XXX
#' @param effects XXX
#' @param objects XXX
#' @param burnin XXX
#' @param thining XXX
#' @param a.scaling XXX
#' @param length.p3 XXX
#' @param neighborhood XXX
#' @param numgroups.allowed XXX
#' @param numgroups.simulated XXX
#' @param sizes.allowed XXX
#' @param sizes.simulated XXX
#' @param fixed.estimates XXX
#' @param parallel XXX
#' @param cpus XXX
#' @return XXX
#' @export
run_phase3_single <- function(partition,
                       estimates.phase2, 
                       z.obs, 
                       nodes, 
                       effects, 
                       objects, 
                       burnin, 
                       thining,
                       a.scaling,
                       length.p3, 
                       neighborhood,
                       numgroups.allowed,
                       numgroups.simulated,
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
    
    sfExport("startingestimates", "first.partition", "nodes", "effects", "objects", "burnin", "thining", "length.p3", "cpus", "neighborhood", "numgroups.allowed", "numgroups.simulated", "sizes.allowed", "sizes.simulated")
    res <- sfLapply(1:cpus, fun = function(k) {
      set.seed(k)
      subres <- draw_Metropolis_single(estimates.phase2, first.partition, nodes, effects, objects, burnin, thining, ceiling(length.p3/cpus), neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
      return(subres)
    }
    )
    all.z <- c()
    for(k in 1:cpus) all.z <- rbind(all.z,res[[k]]$draws)
    length.p3 <- cpus * ceiling(length.p3/cpus)
    results.phase3 <- list("draws" = all.z, "last.partition" = res[[cpus]]$last.partition, "all.partitions" = NULL) 
    
  }else{
    
    results.phase3 <- draw_Metropolis_single(estimates.phase2, first.partition, nodes, effects, objects, burnin, thining, length.p3, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
  }
  z.phase3 <- results.phase3$draws
  
  # calculate autocorrelation to check afterhand
  autocors <- rep(0,num.effects)
  for(e in 1:num.effects){
    autocors[e] <- cor(results.phase3$draws[1:(length.p3-1),e],results.phase3$draws[2:length.p3,e])
  }
  print("Autocorrelations in phase 3:")
  print(autocors)

  
  # hack for size constraints
  if(!is.null(sizes.allowed)){
    length.p3 <- nrow(z.phase3)
    print("new length of phase 3")
    print(length.p3)
  }
  
  # calculations of phase 3: mean, sd, se, conv ratios
  res.phase3 <- phase3(estimates.phase2, z.phase3, z.obs, nodes, effects, length.p3, fixed.estimates)
  
  return(list("means" = res.phase3$finalmean, 
              "standard.deviations" = res.phase3$finalsd, 
              "standard.errors" = res.phase3$finalse, 
              "convergence.ratios" = res.phase3$finalconvratios,
              "inv.zcov" = res.phase3$inv.zcov,
              "inv.scaling" = res.phase3$inv.scaling,
              "autocorrelations" = autocors))
  
}


#' Phase 3 wrapper for multiple observation
#'
#' @param partitions XXX
#' @param estimates.phase2 XXX
#' @param z.obs XXX
#' @param presence.tables XXX
#' @param nodes XXX
#' @param effects XXX
#' @param objects XXX
#' @param burnin XXX
#' @param thining XXX
#' @param a.scaling XXX
#' @param length.p3 XXX
#' @param neighborhood XXX
#' @param numgroups.allowed XXX
#' @param numgroups.simulated XXX
#' @param sizes.allowed XXX
#' @param sizes.simulated XXX
#' @param fixed.estimates XXX
#' @param parallel XXX
#' @param cpus XXX
#' @return XXX
#' @export
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
                              length.p3, 
                              neighborhood,
                              numgroups.allowed,
                              numgroups.simulated,
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
    
    sfExport("startingestimates", "first.partitions", "presence.tables", "nodes", "effects", "objects", "burnin", "thining", "length.p3", "cpus", "neighborhood", "numgroups.allowed", "numgroups.simulated", "sizes.allowed", "sizes.simulated")
    res <- sfLapply(1:cpus, fun = function(k) {
      set.seed(k)
      subres <- draw_Metropolis_multiple(estimates.phase2, first.partitions, presence.tables, nodes, effects, objects, burnin, thining, ceiling(length.p3/cpus), neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, return.all.partitions = T)
      return(subres)
    }
    )
    all.z <- c()
    all.partitions <- list()
    for(k in 1:cpus) {
      all.z <- rbind(all.z,res[[k]]$draws)
      all.partitions[[k]] <- res[[k]]$all.partitions
    }
    length.p3 <- cpus * ceiling(length.p3/cpus)
    results.phase3 <- list("draws" = all.z, "last.partitions" = res[[cpus]]$last.partitions, "all.partitions" = all.partitions) 
  
  }else{
    
    results.phase3 <- draw_Metropolis_multiple(estimates.phase2, first.partitions, presence.tables, nodes, effects, objects, burnin, thining, length.p3, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, return.all.partitions = T)
  
  }
  z.phase3 <- results.phase3$draws
  
  # calculate autocorrelation to check afterhand
  autocors <- rep(0,num.effects)
  for(e in 1:num.effects){
    autocors[e] <- cor(results.phase3$draws[1:(length.p3-1),e],results.phase3$draws[2:length.p3,e])
  }
  print("Autocorrelations in phase 3:")
  print(autocors)
  
  # hack for size constraints
  if(!is.null(sizes.allowed)){
    length.p3 <- nrow(z.phase3)
    print("new length of phase 3")
    print(length.p3)
  }
  
  # calculations of phase 3: mean, sd, se, conv ratios
  res.phase3 <- phase3(estimates.phase2, z.phase3, z.obs, nodes, effects, length.p3, fixed.estimates)
  
  return(list("draws" = z.phase3, 
              "means" = res.phase3$finalmean, 
              "standard.deviations" = res.phase3$finalsd, 
              "standard.errors" = res.phase3$finalse, 
              "convergence.ratios" = res.phase3$finalconvratios,
              "inv.zcov" = res.phase3$inv.zcov,
              "inv.scaling" = res.phase3$inv.scaling,
              "autocorrelations" = autocors,
              "all.partitions" = results.phase3$all.partitions))
  
}



#' Core function for Phase 3
#'
#' @param estimates.phase2 XXX
#' @param z.phase3 XXX
#' @param z.obs XXX
#' @param nodes XXX
#' @param effects XXX
#' @param length.p3 XXX
#' @param fixed.estimates XXX
#' @return XXX
#' @export
phase3 <- function(estimates.phase2,
                   z.phase3,
                   z.obs,
                   nodes,
                   effects,
                   length.p3,
                   fixed.estimates){
  
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
  
  #covariance matrix
  inverted_matrices <- calculate_inverted_covariance_and_scaling(estimates.phase2, 
                                                                 z.obs, 
                                                                 nodes, 
                                                                 effects, 
                                                                 objects, 
                                                                 a.scaling = 0,
                                                                 length.phase = length.p3, 
                                                                 z.phase = z.phase3,
                                                                 fixed.estimates)
  inv.zcov <- inverted_matrices$inv.zcov
  inv.scaling <- inverted_matrices$inv.scaling
  
  #finalse <- finalsd / sqrt(length.p3)
  finalse <- sqrt(diag(inv.zcov))
  
  # convergence ratios
  finalconvratios <- (finalmean - z.obs) / finalsd
  
  print("Estimated statistics after phase 3")
  print(finalmean)
  print("Estimates after phase 3")
  print(estimates.phase2)
  
  return(list("finalmean" = finalmean,
              "finalsd" = finalsd,
              "finalse" = finalse,
              "finalconvratios" = finalconvratios,
              "inv.zcov" = inverted_matrices$inv.zcov,
              "inv.scaling" = inverted_matrices$inv.scaling))
}
