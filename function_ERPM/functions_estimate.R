######################################################################
## Simulation and estimation of Exponential Random Partition Models ## 
## Main function of the estimation algorithm                        ##
## Author: Marion Hoffman                                           ##
######################################################################


## Estimation ERPM functions:

estimate_ERPM <- function(partition, # observed partition
                          nodes, # nodeset (data frame)
                          objects, # objects used for statistics calculation (list with a vector "name", and a vector "object")
                          effects, # effects/sufficient statistics (list with a vector "names", and a vector "objects")
                          startingestimates, # first guess for the model parameters
                          multiplicationfactor = 30, # for now, useless
                          gainfactor = 0.1, # numeric used to decrease the size of steps made in the Newton optimization
                          a.scaling = 0.2, # numeric used to reduce the influence of non-diagonal elements in the scaling matrix (for stability)
                          r.truncation.p1 = 2, # numeric used to limit extreme values in the covariance matrix (for stability)
                          r.truncation.p2 = 5, # numeric used to limit extreme values in the covariance matrix (for stability)
                          mini.steps = "normalized", # type of transition in the Metropolis Hastings algorithm, either "normalized", either "self-loops" (take "normalized")
                          burnin = 30, # integer for the number of burn-in steps before sampling
                          thining = 10, # integer for the number of thining steps between sampling
                          length.p1 = 100, # number of samples in phase 1
                          min.iter.p2 = 10, # minimum number of sub-steps in phase 2
                          max.iter.p2 = 200, # maximum number of sub-steps in phase 2
                          num.steps.p2 = 6, # number of optimisation steps in phase 2
                          length.p3 = 1000, # number of samples in phase 3
                          neighborhood = 2, # way of choosing partitions, either 1 (actor swaps) or 2 (merges and divisions)
                          fixed.estimates = NULL, # if some parameters are fixed, list with as many elements as effects, these elements equal a fixed value if needed, or NULL if they should be estimated
                          sizes.allowed = NULL, # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                          sizes.simulated = NULL, # vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
                          double.averaging = F, # option to average the statistics sampled in each sub-step of phase 2
                          inv.zcov = NULL, # initial value of the inverted covariance matrix (if a phase 3 was run before) to bypass the phase 1
                          inv.scaling = NULL) { # initial value of the inverted scaling matrix (if a phase 3 was run before) to bypass the phase 1
  
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
    results.phase1 <- run_phase1_single(startingestimates, z.obs, nodes, effects, objects, burnin, thining, gainfactor, a.scaling, r.truncation.p1, mini.steps, length.p1, neighborhood, fixed.estimates, sizes.allowed, sizes.simulated)
    estimates.phase1 <- results.phase1$estimates
    inv.zcov <- results.phase1$inv.zcov
    inv.scaling <- results.phase1$inv.scaling
  }
  
  # --------- PHASE 2 ---------
  results.phase2 <- run_phase2_single(estimates.phase1, inv.zcov,inv.scaling, z.obs, nodes, effects, objects, burnin, num.steps.p2, gainfactors, r.truncation.p2, mini.steps, min.iter.p2, max.iter.p2, neighborhood, fixed.estimates, sizes.allowed, sizes.simulated, double.averaging)
  estimates.phase2 <- results.phase2$final.estimates
  
  # --------- PHASE 3 ---------
  results.phase3 <- run_phase3_single(estimates.phase2, z.obs, nodes, effects, objects, burnin, thining, a.scaling, mini.steps, length.p3, neighborhood, sizes.allowed, sizes.simulated, fixed.estimates)
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





## Estimation ERPM for multiple observations:

estimate_multipleERPM <- function(partitions, # observed partitions
                          presence.tables, # matrix indicating which actors were present for each observations (mandatory)
                          nodes, # nodeset (data frame)
                          objects, # objects used for statistics calculation (list with a vector "name", and a vector "object")
                          effects, # effects/sufficient statistics (list with a vector "names", and a vector "objects")
                          startingestimates, # first guess for the model parameters
                          multiplicationfactor = 30, # for now, useless
                          gainfactor = 0.1, # numeric used to decrease the size of steps made in the Newton optimization
                          a.scaling = 0.2, # numeric used to reduce the influence of non-diagonal elements in the scaling matrix (for stability)
                          r.truncation.p1 = 2, # numeric used to limit extreme values in the covariance matrix (for stability)
                          r.truncation.p2 = 5, # numeric used to limit extreme values in the covariance matrix (for stability)
                          mini.steps = "normalized", # type of transition in the Metropolis Hastings algorithm, either "normalized", either "self-loops" (take "normalized")
                          burnin = 30, # integer for the number of burn-in steps before sampling
                          thining = 10, # integer for the number of thining steps between sampling
                          length.p1 = 100, # number of samples in phase 1
                          min.iter.p2 = 10, # minimum number of sub-steps in phase 2
                          max.iter.p2 = 200, # maximum number of sub-steps in phase 2
                          num.steps.p2 = 6, # number of optimisation steps in phase 2
                          length.p3 = 1000, # number of samples in phase 3
                          neighborhood = 2, # way of choosing partitions, either 1 (actor swaps) or 2 (merges and divisions)
                          fixed.estimates = NULL, # if some parameters are fixed, list with as many elements as effects, these elements equal a fixed value if needed, or NULL if they should be estimated
                          sizes.allowed = NULL, # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                          sizes.simulated = NULL, # vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
                          double.averaging = F, # option to average the statistics sampled in each sub-step of phase 2
                          inv.zcov = NULL, # initial value of the inverted covariance matrix (if a phase 3 was run before) to bypass the phase 1
                          inv.scaling = NULL) { # initial value of the inverted scaling matrix (if a phase 3 was run before) to bypass the phase 1
  
  z.obs <- computeStatistics_multiple(partitions, presence.tables, nodes, effects, objects)
 
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
    results.phase1 <- run_phase1_multiple(startingestimates, z.obs, presence.tables, nodes, effects, objects, burnin, thining, gainfactor, a.scaling, r.truncation.p1, mini.steps, length.p1, neighborhood, fixed.estimates, sizes.allowed, sizes.simulated)
    estimates.phase1 <- results.phase1$estimates
    inv.zcov <- results.phase1$inv.zcov
    inv.scaling <- results.phase1$inv.scaling
  }
  
  # --------- PHASE 2 ---------
  results.phase2 <- run_phase2_multiple(estimates.phase1, inv.zcov,inv.scaling, z.obs, presence.tables, nodes, effects, objects, burnin, num.steps.p2, gainfactors, r.truncation.p2, mini.steps, min.iter.p2, max.iter.p2, neighborhood, fixed.estimates, sizes.allowed, sizes.simulated, double.averaging)
  estimates.phase2 <- results.phase2$final.estimates
  
  # --------- PHASE 3 ---------
  results.phase3 <- run_phase3_multiple(estimates.phase2, z.obs, presence.tables, nodes, effects, objects, burnin, thining, a.scaling, mini.steps, length.p3, neighborhood, sizes.allowed, sizes.simulated, fixed.estimates)
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



