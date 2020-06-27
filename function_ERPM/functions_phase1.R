run_phase1 <- function(startingestimates, 
                       z.obs, 
                       nodes, 
                       effects, 
                       objects, 
                       burnin, 
                       thining,
                       gainfactor, 
                       mini.steps, 
                       length.p1, 
                       neighborhood,
                       fixed.estimates,
                       sizes.allowed,
                       sizes.simulated) {
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)
  # TODO: figure out what is the length of phase 1 now?
  #length.p1 <- 7 + 3*num.effects
  #length.p1 <- 100
  
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
  
  # simulate the statisticis distribution
  results.phase1 <- draw_Metropolis(startingestimates, first.partition, nodes, effects, objects, burnin, thining, length.p1, mini.steps, neighborhood, sizes.allowed, sizes.simulated)
  z.phase1 <- results.phase1$draws
  
  # hack for size constraints
  if(!is.null(sizes.allowed)){
    length.p1 <- nrow(z.phase1)
    print("new length of phase 1")
    print(length.p1)
  }
  
  # calculate the covariance and scaling matrices
  inverted_matrices <- calculate_inverted_covariance_and_scaling(startingestimates, 
                                                        z.obs, 
                                                        nodes, 
                                                        effects, 
                                                        objects, 
                                                        length.phase = length.p1, 
                                                        z.phase = z.phase1,
                                                        fixed.estimates)
  inv.zcov <- inverted_matrices$inv.zcov
  inv.scaling <- inverted_matrices$inv.scaling
  
  # normal procedure ( no fixed estimates)
  if(is.null(fixed.estimates)){
    
    
    # compute a first rough estimate of statistics by averaging all results of phase
    z.mean <- rep(0, num.effects)
    for(i in 1:length.p1) {
      z.mean <- z.mean + z.phase1[i,]
    }
    z.mean <- z.mean / length.p1
    
    
    # compute truncating factor
    truncst <- 2
    r <- 1
    diff <- (z.mean - z.obs)
    maxratio <-  max(sqrt((t(diff) %*% inv.zcov %*% diff / num.effects)))
    if(maxratio > truncst) {
      r <- truncst / maxratio
    }
    
    # compute new estimates
    estimates.phase1 <- startingestimates - r*inv.scaling%*%(z.mean - z.obs)
    
    
  # fixed parameters procedure
  } else {
    
    # find the indexes to remove from the estimation
    fixed.indexes <- c()
    unfixed.indexes <- c()
    for(e in 1:num.effects){
      if(!is.null(fixed.estimates[[e]])) 
        fixed.indexes <- c(fixed.indexes,e)
      else
        unfixed.indexes <- c(unfixed.indexes,e)
    }
    
    # compute a first rough estimate of statistics by averaging all results of phase 1
    z.mean <- rep(0, length(unfixed.indexes))
    for(i in 1:length.p1) {
      z.mean <- z.mean + z.phase1[i,unfixed.indexes]
    }
    z.mean <- z.mean / length.phase1
    
    # compute truncating factor
    truncst <- 2
    r <- 1
    diff <- (z.mean - z.obs[unfixed.indexes])
    maxratio <-  max(sqrt((t(diff) %*% inv.zcov %*% diff / length(unfixed.indexes))))
    if(maxratio > truncst) {
      r <- truncst / maxratio
    }
    
    # compute new estimates
    estimates.phase1 <- startingestimates
    estimates.phase1[unfixed.indexes] <- estimates.phase1[unfixed.indexes] - r*inv.scaling%*%(z.mean - z.obs[unfixed.indexes])
    
  }
  
  print("Estimated statistics after phase 1")
  print(z.mean)
  print("Estimates after phase 1")
  print(estimates.phase1)
  
  return( list("estimates" = estimates.phase1, 
               "inv.zcov" = inv.zcov,
               "inv.scaling" = inv.scaling) )
  
}


calculate_inverted_covariance_and_scaling <- function(startingestimates, 
                                                      z.obs, 
                                                      nodes, 
                                                      effects, 
                                                      objects, 
                                                      length.phase, 
                                                      z.phase,
                                                      fixed.estimates){
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)
  
  
  # normal procedure ( no fixed estimates)
  if(is.null(fixed.estimates)){
    
    # compute a first rough estimate of statistics by averaging all results of phase
    z.mean <- rep(0, num.effects)
    for(i in 1:length.phase) {
      z.mean <- z.mean + z.phase[i,]
    }
    z.mean <- z.mean / length.phase
    
    # compute the covariance matrix of all results
    z.cov <- matrix(0, nrow=num.effects, ncol=num.effects)
    for(i in 1:length.phase) {
      z.cov <- z.cov + z.phase[i,]%*%t(z.phase[i,])
    }
    z.cov <- z.cov / length.phase
    z.cov <- z.cov - z.mean %*% t(z.mean)
    inv.zcov <- solve(z.cov)
    
    # compute scaling matrix
    scaling <- z.cov
    
    scaling[ row(scaling) != col(scaling) ] <- 0.2 * scaling[ row(scaling) != col(scaling) ]
    inv.scaling <- solve(scaling)

    
    # fixed parameters procedure
  } else {
    
    # find the indexes to remove from the estimation
    fixed.indexes <- c()
    unfixed.indexes <- c()
    for(e in 1:num.effects){
      if(!is.null(fixed.estimates[[e]])) 
        fixed.indexes <- c(fixed.indexes,e)
      else
        unfixed.indexes <- c(unfixed.indexes,e)
    }
    
    # compute a first rough estimate of statistics by averaging all results of phase 1
    z.mean <- rep(0, length(unfixed.indexes))
    for(i in 1:length.phase) {
      z.mean <- z.mean + z.phase[i,unfixed.indexes]
    }
    z.mean <- z.mean / length.phase
    
    # compute the covariance matrix of all results
    z.cov <- matrix(0, nrow=length(unfixed.indexes), 
                    ncol=length(unfixed.indexes))
    for(i in 1:length.phase) {
      z.cov <- z.cov + z.phase[i,length(unfixed.indexes)]%*%t(z.phase[i,length(unfixed.indexes)])
    }
    z.cov <- z.cov / length.phase
    z.cov <- z.cov - z.mean %*% t(z.mean)
    inv.zcov <- solve(z.cov)
    
    # compute scaling matrix
    scaling <- z.cov
    
    # # handle extreme cases
    # if(det(as.matrix(z.cov)) == 0){
    #   z.cov[ row(z.cov) == col(z.cov) ] <- 0.1
    #   scaling[ row(scaling) == col(scaling) ] <- 0.1
    # }
    
    scaling[ row(scaling) != col(scaling) ] <- 0.7 * scaling[ row(scaling) != col(scaling) ]
    inv.scaling <- solve(scaling)

    
  }
  
  return(list(inv.zcov = inv.zcov,
              inv.scaling = inv.scaling))
}

