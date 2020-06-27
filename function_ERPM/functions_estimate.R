## Models for partitions
## Events data: Partitioning Event Model
## Cross sectional data: Exponential Random Partition Model

## Estimation ERPM functions:

estimate_ERPM <- function(partition, 
                          nodes, 
                          objects, 
                          effects, 
                          startingestimates, 
                          multiplicationfactor = 30, 
                          gainfactor = 0.1, 
                          mini.steps = "normalized", 
                          burnin = 30, 
                          thining = 10,
                          length.p1 = 100, 
                          min.iter.p2 = 10, 
                          max.iter.p2 = 200, 
                          num.steps.p2 = 6, 
                          length.p3 = 1000,
                          neighborhood = 2,
                          fixed.estimates = NULL,
                          sizes.allowed = NULL,
                          sizes.simulated = NULL,
                          double.averaging = F,
                          inv.zcov = NULL,
                          inv.scaling = NULL) {
  
  z.obs <- computeStatistics(partition, nodes, effects, objects)
  #density.obs <- sum(adjacency)/(num.nodes*(num.nodes-1))

  # TODO: what is the burn-in here???
  #burnin <- multiplicationfactor * density.obs * (1-density.obs) * num.nodes^2
  
  # TODO: check also this number of steps
  #num.steps <- 6
  #gainfactors <- c(gainfactor, gainfactor/2, gainfactor/4, gainfactor/8, gainfactor/16, gainfactor/32)
  gainfactors <- rep(0,num.steps.p2)
  for(i in 1:num.steps.p2){
    gainfactors[i] <- gainfactor/(2^(i-1))
  }
  
  # replace the starting estimates with a fixed value
  num.effects <- length(effects$names)
  if(!is.null(fixed.estimates)) {
    for(e in 1:num.effects){
      if(!is.null(fixed.estimates[[e]])){
        startingestimates[e] <- fixed.estimates[[e]]
      }
    }
  }

  print("Observed statistics")
  print(z.obs)
  
  print("Burn-in")
  print(burnin)
  
  print("Thining")
  print(thining)
  
  # --------- PHASE 1 ---------
  if(!is.null(inv.zcov)) {
    estimates.phase1 <- startingestimates
  } else {
    results.phase1 <- run_phase1(startingestimates, z.obs, nodes, effects, objects, burnin, thining, gainfactor, mini.steps, length.p1, neighborhood, fixed.estimates, sizes.allowed, sizes.simulated)
    estimates.phase1 <- results.phase1$estimates
    inv.zcov <- results.phase1$inv.zcov
    inv.scaling <- results.phase1$inv.scaling
  }
  
  # --------- PHASE 2 ---------
  results.phase2 <- run_phase2(estimates.phase1, inv.zcov,inv.scaling, z.obs, nodes, effects, objects, burnin, num.steps.p2, gainfactors, mini.steps, min.iter.p2, max.iter.p2, neighborhood, fixed.estimates, sizes.allowed, sizes.simulated, double.averaging)
  estimates.phase2 <- results.phase2$final.estimates
  
  # --------- PHASE 3 ---------
  results.phase3 <- run_phase3(estimates.phase2, z.obs, nodes, effects, objects, burnin, thining, mini.steps, length.p3, neighborhood, sizes.allowed, sizes.simulated, fixed.estimates)
  means <- results.phase3$means
  standard.deviations <- results.phase3$standard.deviations
  standard.errors <- results.phase3$standard.errors
  convergence.ratios <- results.phase3$convergence.ratios
  
  
  # ------ PRINT RESULTS ------
  results <- data.frame(effect = effects$names, 
                        object = effects$objects,
                        est = as.vector(estimates.phase2), 
                        std.err = standard.errors, 
                        conv = convergence.ratios)
  print_results(results)
  
  # ------ KEEP IMPORTANT OBJECTS ------
  objects.phase2 <- results.phase2$all.estimates
  objects.phase3 <- list(inv.zcov = inv.zcov,
                         inv.scaling = inv.scaling)
  
  return(list(results = results,
              objects.phase2 = objects.phase2,
              objects.phase3 = objects.phase3))
}





# JUST PHASE 3
estimate_ERPM_p3 <- function(partition, 
                          nodes, 
                          objects, 
                          effects, 
                          startingestimates, 
                          mini.steps = "normalized", 
                          burnin = 30, 
                          thining = 10,
                          length.p3 = 1000,
                          neighborhood = 1,
                          fixed.estimates = NULL,
                          sizes.allowed = NULL,
                          sizes.simulated = NULL) {
  
  z.obs <- computeStatistics(partition, nodes, effects, objects)
  
  # replace the starting estimates with a fixed value
  num.effects <- length(effects$names)
  if(!is.null(fixed.estimates)) {
    for(e in 1:num.effects){
      if(!is.null(fixed.estimates[[e]])){
        startingestimates[e] <- fixed.estimates[[e]]
      }
    }
  }
  
  
  # --------- PHASE 3 ---------
  results.phase3 <- run_phase3(startingestimates, z.obs, nodes, effects, objects, burnin, thining, mini.steps, length.p3, neighborhood, sizes.allowed, sizes.simulated)
  means <- results.phase3$means
  standard.deviations <- results.phase3$standard.deviations
  standard.errors <- results.phase3$standard.errors
  convergence.ratios <- results.phase3$convergence.ratios
  
  
  # ------ PRINT RESULTS ------
  results <- data.frame(effect = effects$names, 
                        object = effects$objects,
                        est = as.vector(startingestimates), 
                        std.err = standard.errors, 
                        conv = convergence.ratios)
  print_results(results)
  
  return(results)
}
