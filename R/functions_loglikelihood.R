######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Functions used to estimate the log-likelihood and AIC of a       ##
## given model                                                      ##
## Author: Marion Hoffman                                           ##
######################################################################

#' Estimate log likelihood
#'
#' Function to estimate the log likelihood of a model for an observed partition
#'
#' @param partition observed partition
#' @param nodes node set (data frame)
#' @param effects effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param objects objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param theta estimated model parameters
#' @param theta_0 model parameters if all other effects than "num-groups" are fixed to 0 (basic Dirichlet partition model)
#' @param M number of steps in the path-sampling algorithm
#' @param num.steps number of samples in each step
#' @param burnin integer for the number of burn-in steps before sampling
#' @param thining integer for the number of thining steps between sampling
#' @param neighborhoods = c(0.7,0.3,0) way of choosing partitions
#' @param numgroups.allowed = NULL, # vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated = NULL, # vector containing the number of groups simulated
#' @param sizes.allowed = NULL,   vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated = NULL,  vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param logL_0 = NULL, if known, the value of the log likelihood of the basic dirichlet model
#' @param parallel = FALSE, indicating whether the code should be run in parallel
#' @param cpus = 1, number of cpus required for the parallelization
#' @param verbose = FALSE, to print the current step the algorithm is in
#' @return List with the log likelihood , AIC, lambda and the draws
#' @importFrom snowfall sfExport sfLapply
#' @export
estimate_logL <- function(partition, 
                          nodes, 
                          effects, 
                          objects, 
                          theta, 
                          theta_0, 
                          M, 
                          num.steps, 
                          burnin, 
                          thining, 
                          neighborhoods = c(0.7,0.3,0), 
                          numgroups.allowed = NULL, 
                          numgroups.simulated = NULL,
                          sizes.allowed = NULL, 
                          sizes.simulated = NULL,
                          logL_0 = NULL, 
                          parallel = FALSE, 
                          cpus = 1,
                          verbose = FALSE)
{

  if(parallel && cpus != M) warning("Please set the number of cpus equal to the number of steps in the path-sampling algorithm.")

  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)
  z.obs <- computeStatistics(partition, nodes, effects, objects)

  if(is.null(sizes.allowed)) {sizes.allowed <- 1:num.nodes}
  if(is.null(sizes.simulated)) {sizes.simulated <- 1:num.nodes}

  # find a good starting point
  first.partition <- partition

  # value of estimated denominators ratio
  lambda <- 0
  diff_vector <- theta - theta_0

  #store draws to check the "jumping frogs strategy"
  all_draws <- list()

  # path sampling
  if(parallel){

    sfExport("M", "theta", "theta_0", "diff_vector", "first.partition", "nodes", "effects", "objects", "burnin", "thining", "num.steps", "mini.steps", "neighborhoods", "numgroups.allowed", "numgroups.simulated", "sizes.allowed", "sizes.simulated")
    res <- sfLapply(1:M, fun = function(m) {
      theta_m <- m/M * theta + (1-m)/M * theta_0
      draws_m <- draw_Metropolis_single(theta_m, first.partition, nodes, effects, objects, burnin, thining, num.steps, neighborhoods, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
      z_m <- colMeans(draws_m$draws)
      subres <- list(contribution_lambda = diff_vector %*% z_m, draws = draws_m$draws)
      return(subres)
    }
    )
    for(m in 1:M) {
      lambda <- lambda + res[[m]]$contribution_lambda
      all_draws[[m]] <- res[[m]]$draws
    }
  }else{
    for(m in 1:M){
      if(verbose) cat(paste("step",m), "\n")
      theta_m <- m/M * theta + (1-m)/M * theta_0
      draws_m <- draw_Metropolis_single(theta_m, first.partition, nodes, effects, objects, burnin, thining, num.steps,neighborhoods, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
      all_draws[[m]] <- draws_m
      z_m <- colMeans(draws_m$draws)
      lambda <- lambda + diff_vector %*% z_m
    }
  }
  lambda <- 1/M * lambda

  # value of log likelihood for basic Dirichlet model (theta_0)
  if(is.null(logL_0)){
    index_0 <- which(effects$names == "num_groups")
    logL_0 <- log(calculate_proba_Dirichlet_restricted(theta_0[index_0],max(partition),length(partition),min(sizes.allowed),max(sizes.allowed)))
  }

  # estimated value of log likelihood for full model (theta)
  logL <- diff_vector %*% z.obs - lambda + logL_0

  # calculate resulting AIC
  AIC <- -2*logL -2*num.effects

  return(list("logL" = logL,
              "AIC" = AIC,
              "lambda" = lambda,
              "draws" = all_draws))
}


#calculate_logL_Dirichlet <- function(nodes, theta_0,sizes.allowed) {
#
#  num.nodes <- nrow(nodes)
#  denom <- 0
#  for(k in 1:num.nodes){
#    denom <- denom + Stirling2_constraints(num.nodes,k,min(sizes.allowed),max(sizes.allowed)) * exp(theta_0*k)
#  }
#  return(denom)
#}
