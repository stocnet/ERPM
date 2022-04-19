######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Main function of the estimation algorithm                        ##
## Author: Marion Hoffman                                           ##
######################################################################


## Estimation ERPM functions:

#' Estimate ERPM
#'
#'
#' @param partition observed partition
#' @param nodes nodeset (data frame)
#' @param objects objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param effects effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param startingestimates first guess for the model parameters
#' @param multiplicationfactor = 30,  for now, useless
#' @param gainfactor = 0.1, numeric used to decrease the size of steps made in the Newton optimization
#' @param a.scaling = 0.2,  numeric used to reduce the influence of non-diagonal elements in the scaling matrix (for stability)
#' @param r.truncation.p1 = 2 numeric used to limit extreme values in the covariance matrix (for stability)
#' @param r.truncation.p2 = 5,  numeric used to limit extreme values in the covariance matrix (for stability)
#' @param burnin = 30  integer for the number of burn-in steps before sampling
#' @param thining = 10,  integer for the number of thining steps between sampling
#' @param length.p1 = 100,  number of samples in phase 1
#' @param min.iter.p2 = NULL,  minimum number of sub-steps in phase 2
#' @param max.iter.p2 = NULL,  maximum number of sub-steps in phase 2
#' @param multiplication.iter.p2 = 100, value for the lengths of sub-steps in phase 2 (multiplied by  2.52^k)
#' @param num.steps.p2 = 6,  number of optimisation steps in phase 2
#' @param length.p3 = 1000, number of samples in phase 3
#' @param neighborhood = c(0.7,0.3,0),  way of choosing partitions: probability vector (actors swap, merge/division, single actor move)
#' @param fixed.estimates = NULL, if some parameters are fixed, list with as many elements as effects, these elements equal a fixed value if needed, or NULL if they should be estimated
#' @param sizes.allowed = NULL,  vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated = NULL,  vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param double.averaging = F,  option to average the statistics sampled in each sub-step of phase 2
#' @param inv.zcov = NULL,  initial value of the inverted covariance matrix (if a phase 3 was run before) to bypass the phase 1
#' @param inv.scaling = NULL,  initial value of the inverted scaling matrix (if a phase 3 was run before) to bypass the phase 1
#' @param parallel = F,  whether the phase 1 and 3 should be parallelized
#' @param parallel2 = F,  whether there should be several phases 2 run in parallel
#' @param cpus = 1 ,how many cores can be used
#' @return A list estimates 3 phases
#' @export
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
                          burnin = 30, # integer for the number of burn-in steps before sampling
                          thining = 10, # integer for the number of thining steps between sampling
                          length.p1 = 100, # number of samples in phase 1
                          min.iter.p2 = NULL, # minimum number of sub-steps in phase 2
                          max.iter.p2 = NULL, # maximum number of sub-steps in phase 2
                          multiplication.iter.p2 = 100, # value for the lengths of sub-steps in phase 2 (multiplied by  2.52^k)
                          num.steps.p2 = 6, # number of optimisation steps in phase 2
                          length.p3 = 1000, # number of samples in phase 3
                          neighborhood = c(0.7,0.3,0), # way of choosing partitions: probability vector (actors swap, merge/division, single actor move)
                          fixed.estimates = NULL, # if some parameters are fixed, list with as many elements as effects, these elements equal a fixed value if needed, or NULL if they should be estimated
                          sizes.allowed = NULL, # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                          sizes.simulated = NULL, # vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
                          double.averaging = F, # option to average the statistics sampled in each sub-step of phase 2
                          inv.zcov = NULL, # initial value of the inverted covariance matrix (if a phase 3 was run before) to bypass the phase 1
                          inv.scaling = NULL, # initial value of the inverted scaling matrix (if a phase 3 was run before) to bypass the phase 1
                          parallel = F, # whether the phase 1 and 3 should be parallelized
                          parallel2 = F, # whether there should be several phases 2 run in parallel
                          cpus = 1) { # how many cores can be used

  z.obs <- computeStatistics(partition, nodes, effects, objects)

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
    autocorrelations.phase1 <- NULL
  } else {
    results.phase1 <- run_phase1_single(partition, startingestimates, z.obs, nodes, effects, objects, burnin, thining, gainfactor, a.scaling, r.truncation.p1, length.p1, neighborhood, fixed.estimates, sizes.allowed, sizes.simulated, parallel, cpus)
    estimates.phase1 <- results.phase1$estimates
    inv.zcov <- results.phase1$inv.zcov
    inv.scaling <- results.phase1$inv.scaling
    autocorrelations.phase1 <- results.phase1$autocorrelations
  }

  # --------- PHASE 2 ---------
  if(parallel2){

    sfExport("partition", "estimates.phase1", "inv.zcov", "inv.scaling", "z.obs", "nodes", "effects", "objects", "burnin", "thining", "num.steps.p2", "gainfactors", "r.truncation.p2", "min.iter.p2", "max.iter.p2", "neighborhood", "fixed.estimates", "sizes.allowed", "sizes.simulated", "double.averaging")
    res <- sfLapply(1:cpus, fun = function(k) {
      set.seed(k)
      subres <- run_phase2_single(partition, estimates.phase1, inv.zcov,inv.scaling, z.obs, nodes, effects, objects, burnin, thining, num.steps.p2, gainfactors, r.truncation.p2, min.iter.p2, max.iter.p2, multiplication.iter.p2, neighborhood, fixed.estimates, sizes.allowed, sizes.simulated, double.averaging)
      return(subres)
    }
    )
    final.estimates <- c()
    for(k in 1:cpus) final.estimates <- final.estimates + res[[k]]$final.estimates
    estimates.phase2 <- final.estimates / cpus

  }else{

    results.phase2 <- run_phase2_single(partition, estimates.phase1, inv.zcov,inv.scaling, z.obs, nodes, effects, objects, burnin, thining, num.steps.p2, gainfactors, r.truncation.p2, min.iter.p2, max.iter.p2, multiplication.iter.p2, neighborhood, fixed.estimates, sizes.allowed, sizes.simulated, double.averaging)
    estimates.phase2 <- results.phase2$final.estimates
  }

  # --------- PHASE 3 ---------
  results.phase3 <- run_phase3_single(partition, estimates.phase2, z.obs, nodes, effects, objects, burnin, thining, a.scaling, length.p3, neighborhood, sizes.allowed, sizes.simulated, fixed.estimates, parallel, cpus)
  means <- results.phase3$means
  standard.deviations <- results.phase3$standard.deviations
  standard.errors <- results.phase3$standard.errors
  convergence.ratios <- results.phase3$convergence.ratios
  autocorrelations.phase3 <- results.phase3$autocorrelations


  # ------ PRINT RESULTS ------
  results <- data.frame(effect = effects$names,
                        object = effects$objects,
                        est = as.vector(estimates.phase2),
                        std.err = standard.errors,
                        conv = convergence.ratios)
  print_results(results)

  # ------ KEEP IMPORTANT OBJECTS ------
  objects.phase1 <- list(autocorrelations = autocorrelations.phase1)
  objects.phase2 <- list(estimates = results.phase2$all.estimates,
                         lengths.subphases = results.phase2$lengths.subphases)
  objects.phase3 <- list(inv.zcov = inv.zcov,
                         inv.scaling = inv.scaling,
                         autocorrelations = autocorrelations.phase3)

  return(list(results = results,
              objects.phase1 = objects.phase1,
              objects.phase2 = objects.phase2,
              objects.phase3 = objects.phase3))
}





# JUST PHASE 3


#' Estimate ERPM phase 3
#'
#'
#' @param partition observed partition
#' @param nodes nodeset (data frame)
#' @param objects objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param effects effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param startingestimates first guess for the model parameters
#' @param multiplicationfactor = 30,  for now, useless
#' @param burnin = 30  integer for the number of burn-in steps before sampling
#' @param thining = 10,  integer for the number of thining steps between sampling
#' @param length.p3 = 1000, number of samples in phase 3
#' @param neighborhood = c(0.7,0.3,0),  way of choosing partitions: probability vector (actors swap, merge/division, single actor move)
#' @param fixed.estimates = NULL, if some parameters are fixed, list with as many elements as effects, these elements equal a fixed value if needed, or NULL if they should be estimated
#' @param sizes.allowed = NULL,  vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated = NULL,  vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @return A estimates phase 3
#' @export
estimate_ERPM_p3 <- function(partition,
                          nodes,
                          objects,
                          effects,
                          startingestimates,
                          burnin = 30,
                          thining = 10,
                          length.p3 = 1000,
                          neighborhood = c(0.7,0.3,0),
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
  results.phase3 <- run_phase3(partition, startingestimates, z.obs, nodes, effects, objects, burnin, thining, length.p3, neighborhood, sizes.allowed, sizes.simulated)
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

#' Estimate ERPM for multiple observations
#'
#'
#' @param partitions observed partitions
#' @param presence.tables XXX
#' @param nodes nodeset (data frame)
#' @param objects objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param effects effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param startingestimates first guess for the model parameters
#' @param multiplicationfactor = 30,  for now, useless
#' @param gainfactor = 0.1, numeric used to decrease the size of steps made in the Newton optimization
#' @param a.scaling = 0.2,  numeric used to reduce the influence of non-diagonal elements in the scaling matrix (for stability)
#' @param r.truncation.p1 = 2 numeric used to limit extreme values in the covariance matrix (for stability)
#' @param r.truncation.p2 = 5,  numeric used to limit extreme values in the covariance matrix (for stability)
#' @param burnin = 30  integer for the number of burn-in steps before sampling
#' @param thining = 10,  integer for the number of thining steps between sampling
#' @param length.p1 = 100,  number of samples in phase 1
#' @param min.iter.p2 = NULL,  minimum number of sub-steps in phase 2
#' @param max.iter.p2 = NULL,  maximum number of sub-steps in phase 2
#' @param multiplication.iter.p2 = 100, value for the lengths of sub-steps in phase 2 (multiplied by  2.52^k)
#' @param num.steps.p2 = 6,  number of optimisation steps in phase 2
#' @param length.p3 = 1000, number of samples in phase 3
#' @param neighborhood = c(0.7,0.3,0),  way of choosing partitions: probability vector (actors swap, merge/division, single actor move)
#' @param fixed.estimates = NULL, if some parameters are fixed, list with as many elements as effects, these elements equal a fixed value if needed, or NULL if they should be estimated
#' @param sizes.allowed = NULL,  vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated = NULL,  vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param double.averaging = F,  option to average the statistics sampled in each sub-step of phase 2
#' @param inv.zcov = NULL,  initial value of the inverted covariance matrix (if a phase 3 was run before) to bypass the phase 1
#' @param inv.scaling = NULL,  initial value of the inverted scaling matrix (if a phase 3 was run before) to bypass the phase 1
#' @param parallel = F,  whether the phase 1 and 3 should be parallelized
#' @param parallel2 = F,  whether there should be several phases 2 run in parallel
#' @param cpus = 1 ,how many cores can be used
#' @return A list estimates 3 phases
#' @export
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
                          burnin = 30, # integer for the number of burn-in steps before sampling
                          thining = 10, # integer for the number of thining steps between sampling
                          length.p1 = 100, # number of samples in phase 1
                          min.iter.p2 = NULL, # minimum number of sub-steps in phase 2
                          max.iter.p2 = NULL, # maximum number of sub-steps in phase 2
                          multiplication.iter.p2 = 200, # value for the lengths of sub-steps in phase 2 (multiplied by  2.52^k)
                          num.steps.p2 = 6, # number of optimisation steps in phase 2
                          length.p3 = 1000, # number of samples in phase 3
                          neighborhood = c(0.7,0.3,0), # way of choosing partitions: probability vector (actors swap, merge/division, single actor move)
                          fixed.estimates = NULL, # if some parameters are fixed, list with as many elements as effects, these elements equal a fixed value if needed, or NULL if they should be estimated
                          sizes.allowed = NULL, # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                          sizes.simulated = NULL, # vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
                          double.averaging = F, # option to average the statistics sampled in each sub-step of phase 2
                          inv.zcov = NULL, # initial value of the inverted covariance matrix (if a phase 3 was run before) to bypass the phase 1
                          inv.scaling = NULL, # initial value of the inverted scaling matrix (if a phase 3 was run before) to bypass the phase 1
                          parallel = F, # whether the phase 1 and 3 should be parallelized
                          parallel2 = F, # whether there should be several phases 2 run in parallel
                          cpus = 1) { # how many cores can be used

  z.obs <- rowSums( computeStatistics_multiple(partitions, presence.tables, nodes, effects, objects) )

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
    autocorrelations.phase1 <- NULL
  } else {
    results.phase1 <- run_phase1_multiple(partitions, startingestimates, z.obs, presence.tables, nodes, effects, objects, burnin, thining, gainfactor, a.scaling, r.truncation.p1, length.p1, neighborhood, fixed.estimates, sizes.allowed, sizes.simulated, parallel, cpus)
    estimates.phase1 <- results.phase1$estimates
    inv.zcov <- results.phase1$inv.zcov
    inv.scaling <- results.phase1$inv.scaling
    autocorrelations.phase1 <- results.phase1$autocorrelations
  }

  # --------- PHASE 2 ---------
  if(parallel2){

    sfExport("partitions", "estimates.phase1", "inv.zcov", "inv.scaling", "z.obs", "presence.tables", "nodes", "effects", "objects", "burnin", "thining", "num.steps.p2", "gainfactors", "r.truncation.p2", "min.iter.p2", "max.iter.p2", "neighborhood", "fixed.estimates", "sizes.allowed", "sizes.simulated", "double.averaging")
    res <- sfLapply(1:cpus, fun = function(k) {
      set.seed(k)
      subres <- run_phase2_multiple(partitions, estimates.phase1, inv.zcov,inv.scaling, z.obs, presence.tables, nodes, effects, objects, burnin, thining, num.steps.p2, gainfactors, r.truncation.p2, min.iter.p2, max.iter.p2, multiplication.iter.p2, neighborhood, fixed.estimates, sizes.allowed, sizes.simulated, double.averaging)
      return(subres)
    }
    )
    final.estimates <- c()
    for(k in 1:cpus) final.estimates <- final.estimates + res[[k]]$final.estimates
    estimates.phase2 <- final.estimates / cpus

  }else{

    results.phase2 <- run_phase2_multiple(partitions, estimates.phase1, inv.zcov,inv.scaling, z.obs, presence.tables, nodes, effects, objects, burnin, thining, num.steps.p2, gainfactors, r.truncation.p2, min.iter.p2, max.iter.p2, multiplication.iter.p2, neighborhood, fixed.estimates, sizes.allowed, sizes.simulated, double.averaging)
    estimates.phase2 <- results.phase2$final.estimates
  }


  # --------- PHASE 3 ---------
  results.phase3 <- run_phase3_multiple(partitions, estimates.phase2, z.obs, presence.tables, nodes, effects, objects, burnin, thining, a.scaling, length.p3, neighborhood, sizes.allowed, sizes.simulated, fixed.estimates, parallel, cpus)
  draws <- results.phase3$draws
  means <- results.phase3$means
  standard.deviations <- results.phase3$standard.deviations
  standard.errors <- results.phase3$standard.errors
  convergence.ratios <- results.phase3$convergence.ratios
  autocorrelations.phase3 <- results.phase3$autocorrelations

  # ------ PRINT RESULTS ------
  results <- data.frame(effect = effects$names,
                        object = effects$objects,
                        est = as.vector(estimates.phase2),
                        std.err = standard.errors,
                        conv = convergence.ratios)
  print_results(results)

  # ------ KEEP IMPORTANT OBJECTS ------
  objects.phase1 <- list(autocorrelations = autocorrelations.phase1)
  objects.phase2 <- list(estimates = results.phase2$all.estimates,
                         lengths.subphases = results.phase2$lengths.subphases)
  objects.phase3 <- list(draws = draws,
                         means = means,
                         standard.deviations = standard.deviations,
                         inv.zcov = inv.zcov,
                         inv.scaling = inv.scaling,
                         autocorrelations = autocorrelations.phase3)

  return(list(results = results,
              objects.phase1 = objects.phase1,
              objects.phase2 = objects.phase2,
              objects.phase3 = objects.phase3))
}



## Estimation ERPM for multiple observations, other parallelization (by effect), careful, we need at least as many cpus as effects:

estimate_multipleERPM_secondparallel <- function(partitions, # observed partitions
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
                                  burnin = 30, # integer for the number of burn-in steps before sampling
                                  thining = 10, # integer for the number of thining steps between sampling
                                  length.p1 = 100, # number of samples in phase 1
                                  min.iter.p2 = NULL, # minimum number of sub-steps in phase 2
                                  max.iter.p2 = NULL, # maximum number of sub-steps in phase 2
                                  multiplication.iter.p2 = 200, # value for the lengths of sub-steps in phase 2 (multiplied by  2.52^k)
                                  num.steps.p2 = 6, # number of optimisation steps in phase 2
                                  length.p3 = 1000, # number of samples in phase 3
                                  neighborhood = c(0.7,0.3,0,0,0,0,0,0,0,0,0,0), # way of choosing partitions: probability vector (actors swap, merge/division, single actor move, pair move)
                                  fixed.estimates = NULL, # if some parameters are fixed, list with as many elements as effects, these elements equal a fixed value if needed, or NULL if they should be estimated
                                  sizes.allowed = NULL, # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                                  sizes.simulated = NULL, # vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
                                  double.averaging = F, # option to average the statistics sampled in each sub-step of phase 2
                                  inv.zcov = NULL, # initial value of the inverted covariance matrix (if a phase 3 was run before) to bypass the phase 1
                                  inv.scaling = NULL, # initial value of the inverted scaling matrix (if a phase 3 was run before) to bypass the phase 1
                                  parallel = F, # whether the calculation of stats should be parallelized
                                  cpus = 1) { # how many cores can be used

  z.obs <- rowSums( computeStatistics_multiple(partitions, presence.tables, nodes, effects, objects) )

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
    autocorrelations.phase1 <- NULL
  } else {
    results.phase1 <- run_phase1_multiple_secondparallel(partitions, startingestimates, z.obs, presence.tables, nodes, effects, objects, burnin, thining, gainfactor, a.scaling, r.truncation.p1, length.p1, neighborhood, fixed.estimates, sizes.allowed, sizes.simulated, parallel, cpus)
    estimates.phase1 <- results.phase1$estimates
    inv.zcov <- results.phase1$inv.zcov
    inv.scaling <- results.phase1$inv.scaling
    autocorrelations.phase1 <- results.phase1$autocorrelations
  }

  # --------- PHASE 2 ---------
  results.phase2 <- run_phase2_multiple_secondparallel(partitions, estimates.phase1, inv.zcov,inv.scaling, z.obs, presence.tables, nodes, effects, objects, burnin, thining, num.steps.p2, gainfactors, r.truncation.p2, min.iter.p2, max.iter.p2, multiplication.iter.p2, neighborhood, fixed.estimates, sizes.allowed, sizes.simulated, double.averaging, parallel, cpus)
  estimates.phase2 <- results.phase2$final.estimates


  # --------- PHASE 3 ---------
  results.phase3 <- run_phase3_multiple_secondparallel(partitions, estimates.phase2, z.obs, presence.tables, nodes, effects, objects, burnin, thining, a.scaling, length.p3, neighborhood, sizes.allowed, sizes.simulated, fixed.estimates, parallel, cpus)
  means <- results.phase3$means
  standard.deviations <- results.phase3$standard.deviations
  standard.errors <- results.phase3$standard.errors
  convergence.ratios <- results.phase3$convergence.ratios
  autocorrelations.phase3 <- results.phase3$autocorrelations


  # ------ PRINT RESULTS ------
  results <- data.frame(effect = effects$names,
                        object = effects$objects,
                        est = as.vector(estimates.phase2),
                        std.err = standard.errors,
                        conv = convergence.ratios)
  print_results(results)

  # ------ KEEP IMPORTANT OBJECTS ------
  objects.phase1 <- list(autocorrelations = autocorrelations.phase1)
  objects.phase2 <- list(estimates = results.phase2$all.estimates,
                         lengths.subphases = results.phase2$lengths.subphases)
  objects.phase3 <- list(inv.zcov = inv.zcov,
                         inv.scaling = inv.scaling,
                         autocorrelations = autocorrelations.phase3)

  return(list(results = results,
              objects.phase1 = objects.phase1,
              objects.phase2 = objects.phase2,
              objects.phase3 = objects.phase3))
}



## Estimation ERPM for multiple observations with Bayesian estimation:

estimate_multipleBERPM <- function(partitions, # observed partitions
                                  presence.tables, # matrix indicating which actors were present for each observations (mandatory)
                                  nodes, # nodeset (data frame)
                                  objects, # objects used for statistics calculation (list with a vector "name", and a vector "object")
                                  effects, # effects/sufficient statistics (list with a vector "names", and a vector "objects")
                                  mean.priors, # means of the normal distributions of prior parameters
                                  sd.priors, # standard deviations of the normal distributions of prior parameters
                                  start.chains = NULL, # define a list of starting values for parameters

                                  burnin.1 = 30, # integer for the number of burn-in steps before sampling in the main MCMC chain
                                  thining.1 = 10, # integer for the number of thining steps between sampling in the main MCMC chain
                                  num.chains = 3, # number of MChains
                                  length.chains = 1000, # number of samples of each chain

                                  burnin.2 = 30, # integer for the number of burn-in steps before sampling int the MCMC to sample partitions

                                  neighborhood.partition = c(0.7,0.3,0,0,0,0), # way of choosing partitions: probability vector (actors swap, merge/division, single actor move, pair move)

                                  neighborhood.augmentation = NULL, # standard deviations auround the parameters to draw the augmented distrubtion

                                  sizes.allowed = NULL, # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                                  sizes.simulated = NULL, # vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)

                                  parallel = F, # whether the chains are parallelized (possibly within chains too)
                                  cpus = 1) { # how many cores can be used

  num.effects <- length(effects$names)
  z.obs <- rowSUms( computeStatistics_multiple(partitions, presence.tables, nodes, effects, objects) )

  print("Observed statistics")
  print(z.obs)

  print("Burn-in")
  print(burnin.1)

  print("Thining")
  print(thining.1)

  # if the starts of the chains are not given
  if(is.null(start.chains)){
    start.chains <- list()
    for(p in 1:num.effects){
      start.chains[[p]] <- rnorm(num.effects, mean=mean.priors, sd=sd.priors)
    }
  }

  # if proposal for auxiliary distribution is not given
  if(is.null(neighborhood.augmentation)){
    neighborhood.augmentation <- rep(0.1,num.effects)
  }


  # --------- MAIN MCMC: EXCHANGE ALGORITHM ---------
  results_exchange <- draw_exchangealgorithm_multiple(partitions,
                                                      z.obs,
                                                      presence.tables,
                                                      nodes,
                                                      objects,
                                                      effects,
                                                      mean.priors,
                                                      sd.priors,
                                                      start.chains,
                                                      burnin.1,
                                                      thining.1,
                                                      num.chains,
                                                      length.chains,
                                                      burnin.2,
                                                      neighborhood.partition,
                                                      neighborhood.augmentation,
                                                      sizes.allowed,
                                                      sizes.simulated,
                                                      parallel,
                                                      cpus)

  # ------ PRINT RESULTS ------
  results <- data.frame(effect = effects$names,
                        object = effects$objects,
                        post.mean = results_exchange$post.mean,
                        post.sd = results_exchange$post.sd)
  print_results_bayesian(results)

  return(list(results = results,
              all.chains = results_exchange$all.chains))
}

