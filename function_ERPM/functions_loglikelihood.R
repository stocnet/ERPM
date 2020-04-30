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
                          mini.steps,
                          neighborhood,
                          sizes.allowed,
                          sizes.simulated) {
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)
  z.obs <- computeStatistics(partition, nodes, effects, objects)
  
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
  
  # path sampling
  for(m in 1:M){
    print(paste("step",m))
    theta_m <- m/M * theta + (1-m)/M * theta_0
    draws_m <- draw_Metropolis(theta_m, first.partition, nodes, effects, objects, burnin, thining, num.steps, mini.steps, neighborhood, sizes.allowed, sizes.simulated)
    z_m <- colMeans(draws_m$draws)
    lambda <- lambda + diff_vector %*% t(z_m) 
  }
  lambda <- 1/M * lambda
  
  # value of log likelihood for basic Dirichlet model (theta_0)
  logL_0 <- calculate_logL_Dirichlet(nodes, theta_0,sizes.allowed)
  
  # estimated value of log likelihood for full model (theta)
  logL <- diff_vector %*% t(z.obs) - lambda + logL_0
  
  # calculate resulting AIC 
  AIC <- -2*logL -2*num.effects
  
  return(list("logL" = logL,
              "AIC" = AIC))
}

calculate_logL_Dirichlet <- function(nodes, theta_0,sizes.allowed) {
 
  num.nodes <- nrow(nodes)
  denom <- 0
  for(k in 1:num.nodes){
    denom <- denom + Stirling2_constraints(num.nodes,k,min(sizes.allowed),max(sizes.allowed)) * exp(theta_0*k)
  }
  return(denom)
   
}