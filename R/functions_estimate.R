######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Main function of the estimation algorithm                        ##
## Author: Marion Hoffman                                           ##
######################################################################


## Estimation ERPM functions:

#' Estimate ERPM
#'
#' Function to estimate a given model for a given observed partition.
#' All options of the algorithm can be specified here.
#'
#'
#' @param partition observed partition
#' @param nodes nodeset (data frame)
#' @param objects objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param effects effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param startingestimates first guess for the model parameters
#' @param gainfactor numeric used to decrease the size of steps made in the Newton optimization
#' @param a.scaling numeric used to reduce the influence of non-diagonal elements in the scaling matrix (for stability)
#' @param r.truncation.p1 numeric used to limit extreme values in the covariance matrix (for stability)
#' @param r.truncation.p2 numeric used to limit extreme values in the covariance matrix (for stability)
#' @param burnin integer for the number of burn-in steps before sampling
#' @param thining integer for the number of thining steps between sampling
#' @param length.p1 number of samples in phase 1
#' @param min.iter.p2 minimum number of sub-steps in phase 2
#' @param max.iter.p2 maximum number of sub-steps in phase 2
#' @param multiplication.iter.p2 value for the lengths of sub-steps in phase 2 (multiplied by  2.52^k)
#' @param num.steps.p2 number of optimisation steps in phase 2
#' @param length.p3 number of samples in phase 3
#' @param neighborhood way of choosing partitions: probability vector (actors swap, merge/division, single actor move)
#' @param fixed.estimates if some parameters are fixed, list with as many elements as effects, these elements equal a fixed value if needed, or NULL if they should be estimated
#' @param numgroups.allowed vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated vector containing the number of groups simulated
#' @param sizes.allowed vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param double.averaging option to average the statistics sampled in each sub-step of phase 2
#' @param inv.zcov initial value of the inverted covariance matrix (if a phase 3 was run before) to bypass the phase 1
#' @param inv.scaling initial value of the inverted scaling matrix (if a phase 3 was run before) to bypass the phase 1
#' @param parallel whether the phase 1 and 3 should be parallelized
#' @param parallel2 whether there should be several phases 2 run in parallel
#' @param cpus how many cores can be used
#' @param verbose logical: should intermediate results during the estimation be printed or not? Defaults to FALSE.
#' @return A list with the outputs of the three different phases of the algorithm
#' @export
#' @examples
#' # define an arbitrary set of n = 6 nodes with attributes, and an arbitrary covariate matrix
#' n <- 6 
#' nodes <- data.frame(label = c("A","B","C","D","E","F"),
#'                     gender = c(1,1,2,1,2,2),
#'                     age = c(20,22,25,30,30,31)) 
#' friendship <- matrix(c(0, 1, 1, 1, 0, 0,
#'                        1, 0, 0, 0, 1, 0,
#'                        1, 0, 0, 0, 1, 0,
#'                        1, 0, 0, 0, 0, 0,
#'                        0, 1, 1, 0, 0, 1,
#'                        0, 0, 0, 0, 1, 0), 6, 6, TRUE)
#'
#' # choose the effects to be included (see manual for all effect names)
#' effects <- list(names = c("num_groups","same","diff","tie"),
#' objects = c("partition","gender","age","friendship"))
#' objects <- list()
#' objects[[1]] <- list(name = "friendship", object = friendship)
#' 
#' # define observed partition
#' partition <- c(1,1,2,2,2,3)
#' 
#' \donttest{
#' # estimate
#' startingestimates <- c(-2,0,0,0)
#' estimation <- estimate_ERPM(partition, 
#'                             nodes, 
#'                             objects, 
#'                             effects, 
#'                             startingestimates = startingestimates, 
#'                             burnin = 100, 
#'                             thining = 20,
#'                             length.p1 = 500, # number of samples in phase 1
#'                             
#'                             multiplication.iter.p2 = 20,  
#'                             # factor for the number of iterations in phase 2 subphases
#'                             
#'                             num.steps.p2 = 4, # number of phase 2 subphases
#'                             length.p3 = 1000) # number of samples in phase 3
#' 
#' # get results table
#' estimation
#' }
#' 
estimate_ERPM <- function(partition, 
                          nodes, 
                          objects, 
                          effects, 
                          startingestimates,
                          gainfactor = 0.1, 
                          a.scaling = 0.8,
                          r.truncation.p1 = -1, 
                          r.truncation.p2 = -1, 
                          burnin = 30, 
                          thining = 10, 
                          length.p1 = 100, 
                          min.iter.p2 = NULL,
                          max.iter.p2 = NULL,
                          multiplication.iter.p2 = 100,
                          num.steps.p2 = 6,
                          length.p3 = 1000, 
                          neighborhood = c(0.7,0.3,0), 
                          fixed.estimates = NULL, 
                          numgroups.allowed = NULL,
                          numgroups.simulated = NULL,
                          sizes.allowed = NULL,
                          sizes.simulated = NULL,
                          double.averaging = FALSE,
                          inv.zcov = NULL,
                          inv.scaling = NULL,
                          parallel = FALSE, 
                          parallel2 = FALSE, 
                          cpus = 1,
                          verbose = FALSE) { 

  # calculate observed statistics
  z.obs <- computeStatistics(partition, nodes, effects, objects)

  # set gain factors for phase 2
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

  if (verbose) {
    cat("Observed statistics\n")
    cat(z.obs, "\n\n")
    
    cat("Burn-in\n")
    cat(burnin, "\n\n")
    
    cat("Thining\n")
    cat(thining, "\n\n")
  }
  
  # --------- PHASE 1 ---------
  if(!is.null(inv.zcov)) {
    estimates.phase1 <- startingestimates
    autocorrelations.phase1 <- NULL
  } else {
    results.phase1 <- run_phase1_single(partition, startingestimates, z.obs, nodes, effects, objects, burnin, thining, gainfactor, a.scaling, r.truncation.p1, length.p1, neighborhood, fixed.estimates, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, parallel, cpus, verbose)
    estimates.phase1 <- results.phase1$estimates
    inv.zcov <- results.phase1$inv.zcov
    inv.scaling <- results.phase1$inv.scaling
    autocorrelations.phase1 <- results.phase1$autocorrelations
  }

  # --------- PHASE 2 ---------
  results.phase2 <- run_phase2_single(partition, estimates.phase1, inv.zcov,inv.scaling, z.obs, nodes, effects, objects, burnin, thining, num.steps.p2, gainfactors, r.truncation.p2, min.iter.p2, max.iter.p2, multiplication.iter.p2, neighborhood, fixed.estimates, numgroups.allowed, numgroups.simulated,sizes.allowed, sizes.simulated, double.averaging, parallel2, cpus, verbose)
  estimates.phase2 <- results.phase2$final.estimates

  # --------- PHASE 3 ---------
  results.phase3 <- run_phase3_single(partition, estimates.phase2, z.obs, nodes, effects, objects, burnin, thining, a.scaling, length.p3, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, fixed.estimates, parallel, cpus, verbose)
  means <- results.phase3$means
  standard.deviations <- results.phase3$standard.deviations
  standard.errors <- results.phase3$standard.errors
  convergence.ratios <- results.phase3$convergence.ratios
  autocorrelations.phase3 <- results.phase3$autocorrelations


  # ------ EXTRACT RESULTS ------
  results <- data.frame(effect = effects$names,
                        object = effects$objects,
                        est = as.vector(estimates.phase2),
                        std.err = standard.errors,
                        conv = convergence.ratios)

  # ------ KEEP IMPORTANT OBJECTS ------
  objects.phase1 <- list(autocorrelations = autocorrelations.phase1)
  objects.phase2 <- list(estimates = results.phase2$all.estimates,
                         lengths.subphases = results.phase2$lengths.subphases)
  objects.phase3 <- list(inv.zcov = inv.zcov,
                         inv.scaling = inv.scaling,
                         autocorrelations = autocorrelations.phase3)

  results.list <- (list(results = results,
                        objects.phase1 = objects.phase1,
                        objects.phase2 = objects.phase2,
                        objects.phase3 = objects.phase3))
  
  class(results.list) <- "results.list.erpm"
  
  return(results.list)
}





# # JUST PHASE 3 - FOR NOW DEPRECATED
# 
# #' Estimate ERPM phase 3
# #'
# #' Function to run only the phase 3 of the estimation algorithm, for a given model, a given observed partition, and estimates (from phase 2).
# #' All options of the algorithm can be specified here.
# #'
# #' @param partition observed partition
# #' @param nodes nodeset (data frame)
# #' @param objects objects used for statistics calculation (list with a vector "name", and a vector "object")
# #' @param effects effects/sufficient statistics (list with a vector "names", and a vector "objects")
# #' @param startingestimates first guess for the model parameters
# #' @param burnin integer for the number of burn-in steps before sampling
# #' @param thining integer for the number of thining steps between sampling
# #' @param length.p3 number of samples in phase 3
# #' @param neighborhood way of choosing partitions: probability vector (actors swap, merge/division, single actor move)
# #' @param fixed.estimates if some parameters are fixed, list with as many elements as effects, these elements equal a fixed value if needed, or NULL if they should be estimated
# #' @param numgroups.allowed vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
# #' @param numgroups.simulated vector containing the number of groups simulated
# #' @param sizes.allowed vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
# #' @param sizes.simulated vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
# #' @return A list with the outputs of the phase 3 of the algorithm
# #' @export
# estimate_ERPM_p3 <- function(partition,
#                           nodes,
#                           objects,
#                           effects,
#                           startingestimates,
#                           burnin = 30,
#                           thining = 10,
#                           length.p3 = 1000,
#                           neighborhood = c(0.7,0.3,0),
#                           fixed.estimates = NULL,
#                           numgroups.allowed = NULL,
#                           numgroups.simulated = NULL,
#                           sizes.allowed = NULL,
#                           sizes.simulated = NULL) {
#                           
#   z.obs <- computeStatistics(partition, nodes, effects, objects)
#   
#   # replace the starting estimates with a fixed value
#   num.effects <- length(effects$names)
#   if(!is.null(fixed.estimates)) {
#     for(e in 1:num.effects){
#       if(!is.null(fixed.estimates[[e]])){
#         startingestimates[e] <- fixed.estimates[[e]]
#       }
#     }
#   }
#   
#   
#   # --------- PHASE 3 ---------
#   results.phase3 <- run_phase3_single(partition, startingestimates, z.obs, nodes, effects, objects, burnin, thining, length.p3, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
#   means <- results.phase3$means
#   standard.deviations <- results.phase3$standard.deviations
#   standard.errors <- results.phase3$standard.errors
#   convergence.ratios <- results.phase3$convergence.ratios
#   
#   
#   # ------ EXTRACT RESULTS ------
#   results <- data.frame(effect = effects$names,
#                         object = effects$objects,
#                         est = as.vector(startingestimates),
#                         std.err = standard.errors,
#                         conv = convergence.ratios)
#                         
#   class(results) <- "results.p3.erpm"
#   
#   return(results)
# }





## Estimation ERPM for multiple observations:

#' Estimate ERPM for multiple observations
#'
#' Function to estimate a given model for given observed (multiple) partitions.
#' All options of the algorithm can be specified here.
#'
#' @param partitions observed partitions
#' @param presence.tables XXX
#' @param nodes nodeset (data frame)
#' @param objects objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param effects effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param startingestimates first guess for the model parameters
#' @param gainfactor numeric used to decrease the size of steps made in the Newton optimization
#' @param a.scaling numeric used to reduce the influence of non-diagonal elements in the scaling matrix (for stability)
#' @param r.truncation.p1 numeric used to limit extreme values in the covariance matrix (for stability)
#' @param r.truncation.p2 numeric used to limit extreme values in the covariance matrix (for stability)
#' @param burnin integer for the number of burn-in steps before sampling
#' @param thining integer for the number of thining steps between sampling
#' @param length.p1 number of samples in phase 1
#' @param min.iter.p2 minimum number of sub-steps in phase 2
#' @param max.iter.p2 maximum number of sub-steps in phase 2
#' @param multiplication.iter.p2 value for the lengths of sub-steps in phase 2 (multiplied by  2.52^k)
#' @param num.steps.p2 number of optimisation steps in phase 2
#' @param length.p3 number of samples in phase 3
#' @param neighborhood way of choosing partitions: probability vector (actors swap, merge/division, single actor move)
#' @param fixed.estimates if some parameters are fixed, list with as many elements as effects, these elements equal a fixed value if needed, or NULL if they should be estimated
#' @param numgroups.allowed vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated vector containing the number of groups simulated
#' @param sizes.allowed vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param double.averaging option to average the statistics sampled in each sub-step of phase 2
#' @param inv.zcov initial value of the inverted covariance matrix (if a phase 3 was run before) to bypass the phase 1
#' @param inv.scaling initial value of the inverted scaling matrix (if a phase 3 was run before) to bypass the phase 1
#' @param parallel whether the phase 1 and 3 should be parallelized
#' @param parallel2 whether there should be several phases 2 run in parallel
#' @param cpus how many cores can be used
#' @param verbose logical: should intermediate results during the estimation be printed or not? Defaults to FALSE.
#' @return A list with the outputs of the three different phases of the algorithm
#' @export
#' @examples
#' # define an arbitrary set of n = 6 nodes with attributes, and an arbitrary covariate matrix
#' n <- 6 
#' nodes <- data.frame(label = c("A","B","C","D","E","F"),
#'                     gender = c(1,1,2,1,2,2),
#'                     age = c(20,22,25,30,30,31)) 
#' friendship <- matrix(c(0, 1, 1, 1, 0, 0,
#'                        1, 0, 0, 0, 1, 0,
#'                        1, 0, 0, 0, 1, 0,
#'                        1, 0, 0, 0, 0, 0,
#'                        0, 1, 1, 0, 0, 1,
#'                        0, 0, 0, 0, 1, 0), 6, 6, TRUE) 
#' 
#' # specify whether nodes are present at different points of time
#' presence.tables <- matrix(c(1, 1, 1, 1, 1, 1,
#'                             0, 1, 1, 1, 1, 1,
#'                             1, 0, 1, 1, 1, 1), 6, 3)
#' 
#' # choose effects to be included in the estimated model
#' effects <- list(names = c("num_groups","same","diff","tie","inertia_1"),
#'                 objects = c("partitions","gender","age","friendship","partitions"),
#'                 objects2 = c("","","","",""))
#' objects <- list()
#' objects[[1]] <- list(name = "friendship", object = friendship)
#' 
#' # define the observation
#' partitions <- matrix(c(1, 1, 2, 2, 2, 3,
#'                        NA, 1, 1, 2, 2, 2,
#'                        1, NA, 2, 3, 3, 1), 6, 3) 
#' 
#' \donttest{
#' # estimate
#' startingestimates <- c(-2,0,0,0,0)
#' estimation <- estimate_multipleERPM(partitions,
#'                                     presence.tables,          
#'                                     nodes, 
#'                                     objects, 
#'                                     effects, 
#'                                     startingestimates = startingestimates, 
#'                                     burnin = 100, 
#'                                     thining = 50,
#'                                     gainfactor = 0.6,
#'                                     length.p1 = 200, 
#'                                     multiplication.iter.p2 = 20, 
#'                                     num.steps.p2 = 4, 
#'                                     length.p3 = 1000) 
#' 
#' # get results table
#' estimation
#' }
#' 
estimate_multipleERPM <- function(partitions, 
                                  presence.tables, 
                                  nodes, 
                                  objects, 
                                  effects, 
                                  startingestimates, 
                                  gainfactor = 0.1, 
                                  a.scaling = 0.8, 
                                  r.truncation.p1 = -1, 
                                  r.truncation.p2 = -1, 
                                  burnin = 30, 
                                  thining = 10, 
                                  length.p1 = 100, 
                                  min.iter.p2 = NULL, 
                                  max.iter.p2 = NULL, 
                                  multiplication.iter.p2 = 200, 
                                  num.steps.p2 = 6, 
                                  length.p3 = 1000, 
                                  neighborhood = c(0.7,0.3,0), 
                                  fixed.estimates = NULL, 
                                  numgroups.allowed = NULL, 
                                  numgroups.simulated = NULL,
                                  sizes.allowed = NULL, 
                                  sizes.simulated = NULL, 
                                  double.averaging = FALSE, 
                                  inv.zcov = NULL, 
                                  inv.scaling = NULL, 
                                  parallel = FALSE, 
                                  parallel2 = FALSE, 
                                  cpus = 1,
                                  verbose = FALSE) { 

  # calculate observed statistics
  z.obs <- rowSums( computeStatistics_multiple(partitions, presence.tables, nodes, effects, objects) )
  
  # set gain factors for phase 2
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

  if (verbose) {
    cat("Observed statistics\n")
    cat(z.obs, "\n\n")
    
    cat("Burn-in\n")
    cat(burnin, "\n\n")
    
    cat("Thining\n")
    cat(thining, "\n\n")
  }
  
  # --------- PHASE 1 ---------
  if(!is.null(inv.zcov)) {
    estimates.phase1 <- startingestimates
    autocorrelations.phase1 <- NULL
  } else {
    results.phase1 <- run_phase1_multiple(partitions, startingestimates, z.obs, presence.tables, nodes, effects, objects, burnin, thining, gainfactor, a.scaling, r.truncation.p1, length.p1, neighborhood, fixed.estimates, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, parallel, cpus, verbose)
    estimates.phase1 <- results.phase1$estimates
    inv.zcov <- results.phase1$inv.zcov
    inv.scaling <- results.phase1$inv.scaling
    autocorrelations.phase1 <- results.phase1$autocorrelations
  }

  # --------- PHASE 2 ---------
  results.phase2 <- run_phase2_multiple(partitions, estimates.phase1, inv.zcov,inv.scaling, z.obs, presence.tables, nodes, effects, objects, burnin, thining, num.steps.p2, gainfactors, r.truncation.p2, min.iter.p2, max.iter.p2, multiplication.iter.p2, neighborhood, fixed.estimates, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, double.averaging, parallel2, cpus, verbose)
  estimates.phase2 <- results.phase2$final.estimates


  # --------- PHASE 3 ---------
  results.phase3 <- run_phase3_multiple(partitions, estimates.phase2, z.obs, presence.tables, nodes, effects, objects, burnin, thining, a.scaling, length.p3, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, fixed.estimates, parallel, cpus, verbose)
  draws <- results.phase3$draws
  means <- results.phase3$means
  standard.deviations <- results.phase3$standard.deviations
  standard.errors <- results.phase3$standard.errors
  convergence.ratios <- results.phase3$convergence.ratios
  autocorrelations.phase3 <- results.phase3$autocorrelations

  # ------ EXTRACT RESULTS ------
  results <- data.frame(effect = effects$names,
                        object = effects$objects,
                        est = as.vector(estimates.phase2),
                        std.err = standard.errors,
                        conv = convergence.ratios)

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

  results.list <- (list(results = results,
                        objects.phase1 = objects.phase1,
                        objects.phase2 = objects.phase2,
                        objects.phase3 = objects.phase3))
  
  class(results.list) <- "results.list.erpm"
  
  return(results.list)
}





## Estimation ERPM for multiple observations with Bayesian estimation:
## Beta version

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

                                  neighborhood.partition = c(0.7,0.3,0), # way of choosing partitions: probability vector (actors swap, merge/division, single actor move, pair move)

                                  neighborhood.augmentation = NULL, # standard deviations auround the parameters to draw the augmented distrubtion

                                  numgroups.allowed = NULL, # vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
                                  numgroups.simulated = NULL, # vector containing the number of groups simulated
                                  sizes.allowed = NULL, # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                                  sizes.simulated = NULL, # vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)

                                  parallel = FALSE, # whether the chains are parallelized (possibly within chains too)
                                  cpus = 1, # how many cores can be used
                                  verbose = FALSE) {

  num.effects <- length(effects$names)
  z.obs <- rowSums( computeStatistics_multiple(partitions, presence.tables, nodes, effects, objects) )

  if (verbose) {
    cat("Observed statistics\n")
    cat(z.obs, "\n\n")
    
    cat("Burn-in\n")
    cat(burnin.1, "\n\n")
    
    cat("Thining\n")
    cat(thining.1, "\n\n")
  }
  
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
                                                      numgroups.allowed,
                                                      numgroups.simulated,
                                                      sizes.allowed,
                                                      sizes.simulated,
                                                      parallel,
                                                      cpus)

  # ------ EXTRACT RESULTS ------
  results <- data.frame(effect = effects$names,
                        object = effects$objects,
                        post.mean = results_exchange$post.mean,
                        post.sd = results_exchange$post.sd)

  results.bayesian <- (list(results = results,
                            all.chains = results_exchange$all.chains))
  
  class(results.bayesian) <- "results.bayesian.erpm"
  
  return(results.bayesian)
}

