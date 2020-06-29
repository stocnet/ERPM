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
                          mini.steps = "normalized", # type of transition in the Metropolis Hastings algorithm, either "normalized", either "self-loops" (take "normalized")
                          neighborhood = 2, # way of choosing partitions, either 1 (actor swaps) or 2 (merges and divisions)
                          sizes.allowed = NULL,  # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                          sizes.simulated = NULL ) # vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
{
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)
  z.obs <- computeStatistics(partition, nodes, effects, objects)
  
  if(is.null(sizes.allowed)) {sizes.allowed <- 1:num.nodes}
  if(is.null(sizes.simulated)) {sizes.simulated <- 1:num.nodes}
  
  # find a good starting point
  if(is.null(sizes.allowed)){
    first.partition <- 1 + rbinom(num.nodes, as.integer(num.nodes/2), 0.5)
  } else {
    smin <- min(sizes.allowed)
    cpt <- 0
    g <- 1
    first.partition <- rep(0,num.nodes)
    for(i in 1:num.nodes){
      if(cpt == smin) {
        g <- g + 1
        first.partition[i] <- g
        cpt <- 1
      } else {
        first.partition[i] <- g
        cpt <- cpt + 1
      }
    }
  }
  first.partition <- order_groupids(first.partition)
  
  # value of estimated denominators ratio
  lambda <- 0
  diff_vector <- theta - theta_0
  
  #store draws to check the "jumping frogs strategy"
  all_draws <- list()
  
  # path sampling
  for(m in 1:M){
    print(paste("step",m))
    theta_m <- m/M * theta + (1-m)/M * theta_0
    draws_m <- draw_Metropolis(theta_m, first.partition, nodes, effects, objects, burnin, thining, num.steps, mini.steps, neighborhood, sizes.allowed, sizes.simulated)
    all_draws[[m]] <- draws_m
    z_m <- colMeans(draws_m$draws)
    lambda <- lambda + diff_vector %*% z_m 
  }
  lambda <- 1/M * lambda
  
  # value of log likelihood for basic Dirichlet model (theta_0)
  index_0 <- which(effects$names == "num_groups")
  logL_0 <- log(exp(theta_0[index_0]*max(partition))/calculate_logL_Dirichlet(nodes, theta_0[index_0],sizes.allowed))
  
  # estimated value of log likelihood for full model (theta)
  logL <- diff_vector %*% z.obs - lambda + logL_0
  
  # calculate resulting AIC 
  AIC <- -2*logL -2*num.effects
  
  return(list("logL" = logL,
              "AIC" = AIC,
              "draws" = all_draws))
}

calculate_logL_Dirichlet <- function(nodes, theta_0,sizes.allowed) {
 
  num.nodes <- nrow(nodes)
  denom <- 0
  for(k in 1:num.nodes){
    denom <- denom + Stirling2_constraints(num.nodes,k,min(sizes.allowed),max(sizes.allowed)) * exp(theta_0*k)
  }
  return(denom)
   
}
