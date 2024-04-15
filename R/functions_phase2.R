######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Functions used to run the phase 2 of the estimation algorithm    ##
## Author: Marion Hoffman                                           ##
######################################################################



#' Phase 2 wrapper for single observation
#'
#' @param partition observed partition
#' @param estimates.phase1 vector containing parameter values after phase 1
#' @param inv.zcov inverted covariance matrix
#' @param inv.scaling scaling matrix
#' @param z.obs observed statistics
#' @param nodes node set (data frame)
#' @param effects effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param objects objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param burnin integer for the number of burn-in steps before sampling
#' @param thining integer for the number of thining steps between sampling
#' @param num.steps number of sub-phases in phase 2
#' @param gainfactors vector of gain factors
#' @param r.truncation.p2 truncation factor 
#' @param min.iter minimum numbers of steps in each subphase
#' @param max.iter maximum numbers of steps in each subphase
#' @param multiplication.iter used to calculate min.iter and max.iter if not specified
#' @param neighborhood vector for the probability of choosing a particular transition in the chain
#' @param fixed.estimates if some parameters are fixed, list with as many elements as effects, these elements equal a fixed value if needed, or NULL if they should be estimated
#' @param numgroups.allowed vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated vector containing the number of groups simulated
#' @param sizes.allowed vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param double.averaging boolean to indicate whether we follow the double-averaging procedure (often leads to better convergence)
#' @param parallel boolean to indicate whether the code should be run in parallel
#' @param cpus number of cpus if parallel = TRUE
#' @return a list
#' @importFrom snowfall sfExport sfLapply
#' @export
run_phase2_single <- function(partition,
                       estimates.phase1,
                       inv.zcov,
                       inv.scaling,
                       z.obs,
                       nodes,
                       effects,
                       objects,
                       burnin,
                       thining,
                       num.steps,
                       gainfactors,
                       r.truncation.p2,
                       min.iter,
                       max.iter,
                       multiplication.iter,
                       neighborhood,
                       fixed.estimates,
                       numgroups.allowed,
                       numgroups.simulated,
                       sizes.allowed,
                       sizes.simulated,
                       double.averaging,
                       parallel = FALSE,
                       cpus = 1) {

  num.effects <- length(effects$names)
  num.nodes <- nrow(nodes)
  estimates <- estimates.phase1
  all.estimates <- c()

  # length of subphases
  if(is.null(min.iter) && is.null(max.iter)){
    min.iter <- rep(0,num.steps)
    max.iter <- rep(0,num.steps)
    for(i in 1:num.steps){
      min.iter[i] <- multiplication.iter * (2.52)^(i-1)
      max.iter[i] <- multiplication.iter * (2.52)^(i-1) + 200
    }
  } else {
    min.iter <- rep(min.iter,num.steps)
    max.iter <- rep(max.iter,num.steps)
  }
  lengths.subphases <- rep(0,num.steps)

  # MAIN STEP: iterate over all substeps of phase 2
  for(step in 1:num.steps) {

    # normal procedure: no fixed estimates
    if(is.null(fixed.estimates)) {

      stop.iterations <- FALSE
      i <- 1
      theta.i <- estimates
      all.estimates <- rbind(all.estimates,matrix(theta.i,nrow=1))

      # find a good starting point
      #partition.i <- find_startingpoint_single(nodes, numgroups.allowed, sizes.allowed)
      partition.i <- partition
      sign.i_1 <- rep(0, num.effects) # used to check cross statistics

      # store all stats
      allz <- c()

      # TODO: check whether this is right?
      crossedstats <- rep(FALSE,num.effects)

      # store the estimates collected for all networks simulated after the burn in
      mean.theta <- rep(0,num.effects)
      if(double.averaging) {
        mean.mean.z <- rep(0,num.effects)
        mean.mean.theta <- rep(0,num.effects)
        mean.cpt <- 1
      }

      # SUB STEP: until generated statistics cross the observed ones
      while(!stop.iterations) {

        # draw one element from the chain
        if(parallel){
          
          sfExport("cpus", "theta.i", "partition.i", "nodes", "effects", "objects", "burnin", "neighborhood", "numgroups.allowed", "numgroups.simulated", "sizes.allowed", "sizes.simulated")
          res <- sfLapply(1:cpus, fun = function(k) {
            set.seed(k)
            subres <- draw_Metropolis_single(theta.i, partition.i, nodes, effects, objects, burnin, 1, 1, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
            return(subres)
          }
          )
          all.z <- c()
          for(k in 1:cpus) all.z <- rbind(all.z,res[[k]]$draws)
          z.i <- colMeans(all.z)
          partition.i <- res[[1]]$last.partition
          
        } else {
          
          results.i <- draw_Metropolis_single(theta.i, partition.i, nodes, effects, objects, burnin, 1, 1, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
          z.i <- results.i$draws
          partition.i <- results.i$last.partition
          
        }


        # store z
        allz <- rbind(allz, z.i)

        # compute new parameters
        if(double.averaging) {
          doubleaveraging.step <- compute_parameters_doubleaveraging(z.i, z.obs, theta.i, mean.mean.z, mean.mean.theta, mean.cpt, inv.zcov, inv.scaling, gainfactors[step], r.truncation.p2)
          theta.i <- doubleaveraging.step$theta.i
          mean.mean.z <- doubleaveraging.step$mean.mean.z
          mean.mean.theta <- doubleaveraging.step$mean.mean.theta
          mean.cpt <- doubleaveraging.step$mean.cpt
        } else {
          theta.i <- compute_parameters_simpleaveraging(z.i, z.obs, theta.i, inv.zcov, inv.scaling, gainfactors[step], r.truncation.p2)
        }

        sign.i <- (z.i > z.obs) - (z.i < z.obs)

        # add to the mean
        mean.theta <- mean.theta + theta.i

        # stopping criteria
        if( i > max.iter[step] ){
          stop.iterations <- TRUE
        } else if (i > min.iter[step]) {
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
      #partition.i <- find_startingpoint_single(nodes, numgroups.allowed, sizes.allowed)
      partition.i <- partition
      sign.i_1 <- rep(0, length(unfixed.indexes)) # used to check cross statistics

      # store all stats
      allz <- c()

      # TODO: check whether this is right?
      crossedstats <- rep(FALSE,length(unfixed.indexes))

      # store the estimates collected for all networks simulated after the burn in
      mean.theta <- rep(0,length(unfixed.indexes))

      # Newton-Raphson procedure until generated statistics cross the observed ones
      while(!stop.iterations) {

        # draw chain
        fulltheta.i <- estimates.phase1
        fulltheta.i[unfixed.indexes] <- theta.i
        if(i == 1){
          results.i <- draw_Metropolis_single(fulltheta.i, partition.i, nodes, effects, objects, burnin, 1, 1, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
        } else {
          results.i <- draw_Metropolis_single(fulltheta.i, partition.i, nodes, effects, objects, thining, 1, 1, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
        }
        z.i <- results.i$draws[unfixed.indexes]
        partition.i <- results.i$last.partition

        # store z
        allz <- rbind(allz, z.i)

        # compute new parameters
        if(double.averaging) {
          doubleaveraging.step <- compute_parameters_doubleaveraging(z.i, z.obs[unfixed.indexes], theta.i, mean.mean.z, mean.mean.theta, mean.cpt, inv.zcov, inv.scaling, gainfactors[step], r.truncation.p2)
          theta.i <- double.averaging$theta.i
          mean.mean.z <- double.averaging$mean.mean.z
          mean.mean.theta <- double.averaging$mean.mean.theta
          mean.cpt <- double.averaging$mean.cpt
        } else {
          theta.i <- compute_parameters_simpleaveraging(z.i, z.obs[unfixed.indexes], theta.i, inv.zcov, inv.scaling, gainfactors[step], r.truncation.p2)
        }

        sign.i <- (z.i > z.obs[unfixed.indexes]) - (z.i < z.obs[unfixed.indexes])

        # add to the mean
        mean.theta <- mean.theta + theta.i

        # stopping criteria
        if( i > max.iter[step] ){
          stop.iterations <- TRUE
        } else if (i > min.iter[step]) {
          crossedstats <- crossedstats + (sign.i_1 != sign.i)
          stop.iterations <- all(crossedstats)
        }

        sign.i_1 <- sign.i
        i <- i+1

      }

      mean.theta <- mean.theta / (i-1)
      estimates[unfixed.indexes] <- mean.theta
      lengths.subphases[step] <- i-1

      message(cat("Step",step))
      message(cat("Length of the step",(i-1),"(minimal value:",min.iter[step],"and maximal value:",max.iter[step],")"))
      message(cat("Current estimate",estimates, "\n"))

    }

    message(cat("Difference to estimated statistics after phase 2, step",step))
    message(colMeans(allz) - z.obs)
    message(cat("Estimates after phase 2, step",step))
    message(estimates, "\n")
  }

  return(list(all.estimates = all.estimates,
              final.estimates = estimates,
              lengths.subphases = lengths.subphases))

}



#' Phase 2 wrapper for multiple observation
#'
#' @param partitions observed partitions
#' @param estimates.phase1 vector containing parameter values after phase 1
#' @param inv.zcov inverted covariance matrix
#' @param inv.scaling scaling matrix
#' @param z.obs observed statistics
#' @param presence.tables data frame to indicate which times nodes are present in the partition
#' @param nodes node set (data frame)
#' @param effects effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param objects objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param burnin integer for the number of burn-in steps before sampling
#' @param thining integer for the number of thining steps between sampling
#' @param num.steps number of sub-phases in phase 2
#' @param gainfactors vector of gain factors
#' @param r.truncation.p2 truncation factor 
#' @param min.iter minimum numbers of steps in each subphase
#' @param max.iter maximum numbers of steps in each subphase
#' @param multiplication.iter used to calculate min.iter and max.iter if not specified
#' @param neighborhood vector for the probability of choosing a particular transition in the chain
#' @param fixed.estimates if some parameters are fixed, list with as many elements as effects, these elements equal a fixed value if needed, or NULL if they should be estimated
#' @param numgroups.allowed vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated vector containing the number of groups simulated
#' @param sizes.allowed vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param double.averaging boolean to indicate whether we follow the double-averaging procedure (often leads to better convergence)
#' @param parallel boolean to indicate whether the code should be run in parallel
#' @param cpus number of cpus if parallel = TRUE
#' @return a list
#' @importFrom snowfall sfExport sfLapply
#' @export
run_phase2_multiple <- function(partitions,
                              estimates.phase1,
                              inv.zcov,
                              inv.scaling,
                              z.obs,
                              presence.tables,
                              nodes,
                              effects,
                              objects,
                              burnin,
                              thining,
                              num.steps,
                              gainfactors,
                              r.truncation.p2,
                              min.iter,
                              max.iter,
                              multiplication.iter,
                              neighborhood,
                              fixed.estimates,
                              numgroups.allowed,
                              numgroups.simulated,
                              sizes.allowed,
                              sizes.simulated,
                              double.averaging,
                              parallel = FALSE,
                              cpus = 1) {

  num.effects <- length(effects$names)
  num.nodes <- nrow(nodes)
  num.obs <- ncol(presence.tables)

  estimates <- estimates.phase1
  all.estimates <- c()

  # length of subphases
  if(is.null(min.iter) && is.null(max.iter)){
    min.iter <- rep(0,num.steps)
    max.iter <- rep(0,num.steps)
    for(i in 1:num.steps){
      min.iter[i] <- multiplication.iter * (2.52)^(i-1)
      max.iter[i] <- multiplication.iter * (2.52)^(i-1) + 200
    }
  } else {
    min.iter <- rep(min.iter,num.steps)
    max.iter <- rep(max.iter,num.steps)
  }
  lengths.subphases <- rep(0,num.steps)

  # MAIN STEP: iterate over all substeps of phase 2
  for(step in 1:num.steps) {

    # normal procedure: no fixed estimates
    if(is.null(fixed.estimates)) {

      stop.iterations <- FALSE
      i <- 1
      theta.i <- estimates

      # find a good starting point
      #partitions.i <- find_startingpoint_multiple(presence.tables, nodes, numgroups.allowed, sizes.allowed)
      partitions.i <- partitions
      sign.i_1 <- rep(0, num.effects) # used to check cross statistics

      # store all stats
      allz <- c()

      # TODO: check whether this is right?
      crossedstats <- rep(FALSE,num.effects)

      # store the estimates collected for all networks simulated after the burn in
      mean.theta <- rep(0,num.effects)
      if(double.averaging) {
        mean.mean.z <- rep(0,num.effects)
        mean.mean.theta <- rep(0,num.effects)
        mean.cpt <- 1
      }

      # SUB STEP: until generated statistics cross the observed ones
      while(!stop.iterations) {

        message(theta.i)
        all.estimates <- rbind(all.estimates,matrix(theta.i,nrow=1))

        # draw one element from the chain
        if(parallel){
          
          sfExport("cpus", "theta.i", "partitions.i", "presence.tables", "nodes", "effects", "objects", "burnin", "neighborhood", "numgroups.allowed", "numgroups.simulated", "sizes.allowed", "sizes.simulated")
          res <- sfLapply(1:cpus, fun = function(k) {
            set.seed(k)
            subres <- draw_Metropolis_multiple(theta.i, partitions.i, presence.tables, nodes, effects, objects, burnin, 1, 1, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
            return(subres)
          }
          )
          all.z <- c()
          for(k in 1:cpus) all.z <- rbind(all.z,res[[k]]$draws)
          z.i <- colMeans(all.z)
          partitions.i <- res[[1]]$last.partitions
          
        } else {
          
          results.i <- draw_Metropolis_multiple(theta.i, partitions.i, presence.tables, nodes, effects, objects, burnin, 1, 1, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
          z.i <- results.i$draws
          partitions.i <- results.i$last.partitions
          
        }

        # store z
        allz <- rbind(allz, z.i)

        # compute new parameters
        if(double.averaging) {
          doubleaveraging.step <- compute_parameters_doubleaveraging(z.i, z.obs, theta.i, mean.mean.z, mean.mean.theta, mean.cpt, inv.zcov, inv.scaling, gainfactors[step], r.truncation.p2)
          theta.i <- doubleaveraging.step$theta.i
          mean.mean.z <- doubleaveraging.step$mean.mean.z
          mean.mean.theta <- doubleaveraging.step$mean.mean.theta
          mean.cpt <- doubleaveraging.step$mean.cpt
        } else {
          theta.i <- compute_parameters_simpleaveraging(z.i, z.obs, theta.i, inv.zcov, inv.scaling, gainfactors[step], r.truncation.p2)
        }

        sign.i <- (z.i > z.obs) - (z.i < z.obs)

        # add to the mean
        mean.theta <- mean.theta + theta.i

        # stopping criteria
        if( i > max.iter[step] ){
          stop.iterations <- TRUE
        } else if (i > min.iter[step]) {
          crossedstats <- crossedstats + (sign.i_1 != sign.i)
          stop.iterations <- all(crossedstats)
        }

        sign.i_1 <- sign.i
        i <- i+1

      }

      mean.theta <- mean.theta / (i-1)
      estimates <- mean.theta
      lengths.subphases[step] <- i-1

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
      # partitions.i <- find_startingpoint_multiple(presence.table, nodes, numgroups.allowed, sizes.allowed)
      partitions.i <- partitions
      sign.i_1 <- rep(0, length(unfixed.indexes)) # used to check cross statistics

      # store all stats
      allz <- c()

      # TODO: check whether this is right?
      crossedstats <- rep(FALSE,length(unfixed.indexes))

      # store the estimates collected for all networks simulated after the burn in
      mean.theta <- rep(0,length(unfixed.indexes))

      # Newton-Raphson procedure until generated statistics cross the observed ones
      while(!stop.iterations) {

        # draw chain
        fulltheta.i <- estimates.phase1
        fulltheta.i[unfixed.indexes] <- theta.i
        if(i == 1){
          results.i <- draw_Metropolis_multiple(theta.i, partitions.i, presence.tables, nodes, effects, objects, burnin, 1, 1, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
        } else {
          results.i <- draw_Metropolis_multiple(theta.i, partitions.i, presence.tables, nodes, effects, objects, thining, 1, 1, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
        }
        z.i <- results.i$draws[unfixed.indexes]
        partitions.i <- results.i$last.partitions

        # store z
        allz <- rbind(allz, z.i)

        # compute new parameters
        if(double.averaging) {
          doubleaveraging.step <- compute_parameters_doubleaveraging(z.i, z.obs[unfixed.indexes], theta.i, mean.mean.z, mean.mean.theta, mean.cpt, inv.zcov, inv.scaling, gainfactors[step], r.truncation.p2)
          theta.i <- double.averaging$theta.i
          mean.mean.z <- double.averaging$mean.mean.z
          mean.mean.theta <- double.averaging$mean.mean.theta
          mean.cpt <- double.averaging$mean.cpt
        } else {
          theta.i <- compute_parameters_simpleaveraging(z.i, z.obs[unfixed.indexes], theta.i, inv.zcov, inv.scaling, gainfactors[step], r.truncation.p2)
        }

        sign.i <- (z.i > z.obs[unfixed.indexes]) - (z.i < z.obs[unfixed.indexes])

        # add to the mean
        mean.theta <- mean.theta + theta.i

        # stopping criteria
        if( i > max.iter[step] ){
          stop.iterations <- TRUE
        } else if (i > min.iter[step]) {
          crossedstats <- crossedstats + (sign.i_1 != sign.i)
          stop.iterations <- all(crossedstats)
        }

        sign.i_1 <- sign.i
        i <- i+1

      }

      mean.theta <- mean.theta / (i-1)
      estimates[unfixed.indexes] <- mean.theta
      lengths.subphases[step] <- i-1

    }

    message(cat("Length of the step",step))
    message(cat((i-1),"(minimal value:",min.iter[step],"and maximal value:",max.iter[step],")"), "\n")
    message(cat("Estimated statistics after phase 2, step",step))
    message(colMeans(allz) - z.obs, "\n")
    message(cat("Estimates after phase 2, step",step))
    message(estimates, "\n")
  }

  return(list(all.estimates = all.estimates,
              final.estimates = estimates,
              lengths.subphases = lengths.subphases))

}




# update of parameters in one step of phase 2 (normal procedure)
compute_parameters_simpleaveraging <- function(z.i,
                                      z.obs,
                                      theta.i,
                                      inv.zcov,
                                      inv.scaling,
                                      gainfactor,
                                      r.truncation.p2){

    num.effects <- length(theta.i)

    # compute truncating factor
    r <- 1
    if(r.truncation.p2 > 0){
      diff <- t(z.i - z.obs)
      maxratio <-  max(sqrt((t(diff) %*% inv.zcov %*% diff / num.effects)))
      if(maxratio > r.truncation.p2) {
        r <- r.truncation.p2 / maxratio
      }
    }

    # new theta
    theta.i <- theta.i - gainfactor * r * inv.scaling %*% (z.i - z.obs)

    return(theta.i)
}


# update of parameters in ones tep of phase 2 (double averaging procedure)
compute_parameters_doubleaveraging <- function(z.i,
                                               z.obs,
                                               theta.i,
                                               mean.mean.z,
                                               mean.mean.theta,
                                               mean.cpt,
                                               inv.zcov,
                                               inv.scaling,
                                               gainfactor,
                                               r.truncation.p2){
  num.effects <- length(theta.i)

  # mean.cpt = N, mean.mean.z = average of stats from 1 to N 
  if(mean.cpt == 1) mean.mean.z <- z.i
  if(mean.cpt > 1) mean.mean.z <- (mean.cpt-1) / mean.cpt *  mean.mean.z + z.i / mean.cpt

  # compute truncating factor
  r <- 1
  if(r.truncation.p2 > 0){
    diff <- t(mean.mean.z - z.obs)
    maxratio <-  max(sqrt((t(diff) %*% inv.zcov %*% diff / num.effects)))
    if(maxratio > r.truncation.p2) {
      r <- r.truncation.p2 / maxratio
    }
  }
  
  # mean.mean.theta = average of thetas from 1 to N 
  if(mean.cpt == 1) mean.mean.theta <- theta.i
  if(mean.cpt > 1) mean.mean.theta <- (mean.cpt-1) / mean.cpt * mean.mean.theta + theta.i / mean.cpt
  
  # theta.i (theta_N+1) = [average theta until N] - a_N * N * r * D^-1 * ([average stats until N] - obs stats)  
  theta.i <- mean.mean.theta - gainfactor * mean.cpt * r * inv.scaling %*% (mean.mean.z - z.obs)

  mean.cpt <- mean.cpt + 1

  return(list("theta.i" = theta.i,
              "mean.mean.z" = mean.mean.z,
              "mean.mean.theta" = mean.mean.theta,
              "mean.cpt" = mean.cpt))

}
