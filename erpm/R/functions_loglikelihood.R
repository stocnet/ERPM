######################################################################
## Simulation and estimation of Exponential Random Partition Models ## 
## Functions used to estimate the log-likelihood and AIC of a       ##
## given model                                                      ##
## Author: Marion Hoffman                                           ##
######################################################################



estimate_logL <- function(partition, # observed partition
                          nodes, # nodeset (data frame)
                          effects, # effects/sufficient statistics (list with a vector "names", and a vector "objects")
                          objects, # objects used for statistics calculation (list with a vector "name", and a vector "object")
                          theta, # estimated model parameters
                          theta_0, # model parameters if all other effects than "num-groups" are fixed to 0 (basic Dirichlet partition model)
                          M, # number of steps in the path-sampling algorithm
                          num.steps, # number of samples in each step
                          burnin, # integer for the number of burn-in steps before sampling
                          thining, # integer for the number of thining steps between sampling
                          neighborhoods = c(0.7,0.3,0), # way of choosing partitions
                          sizes.allowed = NULL,  # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                          sizes.simulated = NULL, # vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
                          logL_0 = NULL,  # if known, the value of the log likelihood of the basic dirichlet model
                          parallel = F, # whether each step is run in parallel
                          cpus = 1) # number of cpus (should be equal to M)
{
  
  if(parallel && cpus != M) print("Please set the number of cpus equal to the number of steps in the path-sampling algorithm.")
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)
  z.obs <- computeStatistics(partition, nodes, effects, objects)
  
  if(is.null(sizes.allowed)) {sizes.allowed <- 1:num.nodes}
  if(is.null(sizes.simulated)) {sizes.simulated <- 1:num.nodes}
  
  # find a good starting point
  first.partition <- find_startingpoint_single(nodes,sizes.allowed)
  
  # value of estimated denominators ratio
  lambda <- 0
  diff_vector <- theta - theta_0
  
  #store draws to check the "jumping frogs strategy"
  all_draws <- list()
  
  # path sampling
  if(parallel){
    
    sfExport("M", "theta", "theta_0", "diff_vector", "first.partition", "nodes", "effects", "objects", "burnin", "thining", "num.steps", "mini.steps", "neighborhoods", "sizes.allowed", "sizes.simulated")
    res <- sfLapply(1:M, fun = function(m) {
      theta_m <- m/M * theta + (1-m)/M * theta_0
      draws_m <- draw_Metropolis_single(theta_m, first.partition, nodes, effects, objects, burnin, thining, num.steps, neighborhoods, sizes.allowed, sizes.simulated)
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
      print(paste("step",m))
      theta_m <- m/M * theta + (1-m)/M * theta_0
      draws_m <- draw_Metropolis_single(theta_m, first.partition, nodes, effects, objects, burnin, thining, num.steps,neighborhoods, sizes.allowed, sizes.simulated)
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
