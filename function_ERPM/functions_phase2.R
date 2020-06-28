######################################################################
## Simulation and estimation of Exponential Random Partition Models ## 
## Functions used to run the phase 2 of the estimation algorithm    ##
## Author: Marion Hoffman                                           ##
######################################################################



run_phase2 <- function(estimates.phase1, 
                       inv.zcov,
                       inv.scaling, 
                       z.obs, 
                       nodes, 
                       effects, 
                       objects, 
                       burnin, 
                       num.steps, 
                       gainfactors,
                       r.truncation.p2,
                       mini.steps, 
                       min.iter, 
                       max.iter, 
                       neighborhood,
                       fixed.estimates,
                       sizes.allowed,
                       sizes.simulated,
                       double.averaging) {

  num.effects <- length(effects$names)
  num.nodes <- nrow(nodes)
  estimates <- estimates.phase1
  all.estimates <- c()
  
  # iterate over all substeps of phase 2
  for(step in 1:num.steps) {
    
    # normal procedure: no fixed estimates
    if(is.null(fixed.estimates)) {
      
      stop.iterations <- FALSE
      i <- 1
      theta.i <- estimates
      all.estimates <- rbind(all.estimates,matrix(theta.i,nrow=1)) 
      
      # find a good starting point
      if(is.null(sizes.allowed)){
        partition.i <- 1 + rbinom(num.nodes, as.integer(num.nodes/2), 0.5)
      } else {
        smin <- min(sizes.allowed)
        cpt <- 0
        g <- 1
        partition.i <- rep(0,num.nodes)
        for(j in 1:num.nodes){
          if(cpt == smin) {
            g <- g + 1
            partition.i[j] <- g
            cpt <- 1
          } else {
            partition.i[j] <- g
            cpt <- cpt + 1
          }
        }
      }
      partition.i <- order_groupids(partition.i)
      sign.i_1 <- rep(0, num.effects) # used to check cross statistics
      
      # TODO: check whether this is right?
      crossedstats <- rep(FALSE,num.effects)
      
      # store the estimates collected for all networks simulated after the burn in
      mean.theta <- rep(0,num.effects)
      if(double.averaging) {
        mean.mean.z <- rep(0,num.effects) 
        mean.mean.theta <- rep(0,num.effects)
        mean.cpt <- 0
      }
      
      # Newton-Raphson procedure until generated statistics cross the observed ones
      while(!stop.iterations) {
        
        # draw chain
        results.i <- draw_Metropolis(theta.i, partition.i, nodes, effects, objects, burnin, 1, 1, mini.steps, neighborhood, sizes.allowed, sizes.simulated)
        z.i <- results.i$draws
        partition.i <- results.i$last.partition
        
        # compute new parameters
        if(double.averaging) {
          mean.mean.z <- mean.cpt / (mean.cpt+1) *  mean.mean.z + z.i / (mean.cpt+1)
          
          # compute truncating factor
          r <- 1
          diff <- t(mean.mean.z - z.obs)
          maxratio <-  max(sqrt((t(diff) %*% inv.zcov %*% diff / num.effects)))
          if(maxratio > r.truncation.p2) {
            r <- r.truncation.p2 / maxratio
          }
          
          # new theta
          theta.i <- theta.i - gainfactors[step] * mean.cpt * r * inv.scaling %*% t(mean.mean.z - z.obs)
          
          if(mean.cpt == 0) {
            mean.mean.theta <- theta.i
          } else {
            mean.mean.theta <- mean.cpt / (mean.cpt+1) * mean.mean.theta + theta.i / (mean.cpt+1)
          }
          mean.cpt <- mean.cpt + 1
        } else {
          
          # compute truncating factor
          r <- 1
          diff <- t(z.i - z.obs)
          maxratio <-  max(sqrt((t(diff) %*% inv.zcov %*% diff / num.effects)))
          if(maxratio > r.truncation.p2) {
            r <- r.truncation.p2 / maxratio
          }
          
          # new theta
          theta.i <- theta.i - gainfactors[step] * r * inv.scaling %*% t(z.i - z.obs)
        }
        sign.i <- (z.i > z.obs) - (z.i < z.obs)
        
        # add to the mean
        mean.theta <- mean.theta + theta.i
        
        # stopping criteria
        if( i > max.iter ){
          stop.iterations <- TRUE
        } else if (i > min.iter) {
          crossedstats <- crossedstats + (sign.i_1 != sign.i)
          stop.iterations <- all(crossedstats)
        }
        
        sign.i_1 <- sign.i
        i <- i+1
        
      }
      
      mean.theta <- mean.theta / (i-1)
      estimates <- mean.theta
      
    # fixed estimates procedure  
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
      
      stop.iterations <- FALSE
      i <- 1
      theta.i <- estimates[unfixed.indexes]
      all.estimates <- rbind(all.estimates,matrix(theta.i,nrow=1))
      
      # find a good starting point
      if(is.null(sizes.allowed)){
        partition.i <- 1 + rbinom(num.nodes, as.integer(num.nodes/2), 0.5)
      } else {
        smin <- min(sizes.allowed)
        cpt <- 0
        g <- 1
        partition.i <- rep(0,num.nodes)
        for(j in 1:num.nodes){
          if(cpt == smin) {
            g <- g + 1
            partition.i[j] <- g
            cpt <- 1
          } else {
            partition.i[j] <- g
            cpt <- cpt + 1
          }
        }
      }
      partition.i <- order_groupids(partition.i)
      sign.i_1 <- rep(0, length(unfixed.indexes)) # used to check cross statistics
      
      # TODO: check whether this is right?
      crossedstats <- rep(FALSE,length(unfixed.indexes))
      
      # store the estimates collected for all networks simulated after the burn in
      mean.theta <- rep(0,length(unfixed.indexes))
      
      # Newton-Raphson procedure until generated statistics cross the observed ones
      while(!stop.iterations) {
        
        # compute new parameters
        fulltheta.i <- estimates.phase1
        fulltheta.i[unfixed.indexes] <- theta.i
        results.i <- draw_Metropolis(fulltheta.i, partition.i, nodes, effects, objects, burnin, 1, 1, mini.steps, neighborhood, sizes.allowed, sizes.simulated)
        z.i <- results.i$draws[unfixed.indexes]
        partition.i <- results.i$last.partition
        
        # compute truncating factor
        r <- 1
        diff <- t(z.i - z.obs[unfixed.indexes])
        maxratio <-  max(sqrt((t(diff) %*% inv.zcov %*% diff / length(unfixed.indexes))))
        if(maxratio > r.truncation.p2) {
          r <- r.truncation.p2 / maxratio
        }
        
        theta.i <- theta.i - gainfactors[step] * r * inv.scaling %*% t(z.i - z.obs[unfixed.indexes])
        sign.i <- (z.i > z.obs[unfixed.indexes]) - (z.i < z.obs[unfixed.indexes])
        
        # add to the mean
        mean.theta <- mean.theta + theta.i
        
        # stopping criteria
        if( i > max.iter ){
          stop.iterations <- TRUE
        } else if (i > min.iter) {
          crossedstats <- crossedstats + (sign.i_1 != sign.i)
          stop.iterations <- all(crossedstats)
        }
        
        sign.i_1 <- sign.i
        i <- i+1
        
      }
      
      mean.theta <- mean.theta / (i-1)
      estimates[unfixed.indexes] <- mean.theta
      
      print(cat("step",step))
      print(cat("current estimate",estimates))

    }

    print(cat("Estimated statistics after phase 2, step",step))
    print(z.i)
    print(cat("Estimates after phase 2, step",step))
    print(estimates)
  }
  
  return(list(all.estimates = all.estimates,
              final.estimates = estimates))
  
}
