######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Functions used to run the phase 2 of the estimation algorithm    ##
## Author: Marion Hoffman                                           ##
######################################################################

#' Phase 2 wrapper for single observation
#'
#' @param estimates.phase1 XXX
#' @param inv.zcov XXX
#' @param inv.scaling XXX
#' @param z.obs XXX
#' @param nodes XXX
#' @param effects XXX
#' @param objects XXX
#' @param burnin XXX
#' @param thining XXX
#' @param num.steps XXX
#' @param gainfactors XXX
#' @param r.truncation.p2 XXX
#' @param min.iter XXX
#' @param max.iter XXX
#' @param multiplication.iter XXX
#' @param neighborhood XXX
#' @param fixed.estimates XXX
#' @param sizes.allowed XXX
#' @param sizes.simulated XXX
#' @param double.averaging XXX
#' @return XXX
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
                       sizes.allowed,
                       sizes.simulated,
                       double.averaging) {

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
      #partition.i <- find_startingpoint_single(nodes,sizes.allowed)
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
        mean.cpt <- 0
      }

      # SUB STEP: until generated statistics cross the observed ones
      while(!stop.iterations) {

        # draw chain
        if(i == 1){
          results.i <- draw_Metropolis_single(theta.i, partition.i, nodes, effects, objects, burnin, 1, 1, neighborhood, sizes.allowed, sizes.simulated)
        } else{
          results.i <- draw_Metropolis_single(theta.i, partition.i, nodes, effects, objects, thining, 1, 1, neighborhood, sizes.allowed, sizes.simulated)
        }
        z.i <- results.i$draws
        partition.i <- results.i$last.partition

        # store z
        allz <- rbind(allz, z.i)

        # compute new parameters
        if(double.averaging) {
          doubleaveraging.step <- compute_parameters_doubleaveraging(z.i, z.obs, theta.i, mean.mean.z, mean.mean.theta, mean.cpt, inv.zcov, inv.scaling, gainfactors[step], r.truncation.p2)
          theta.i <- double.averaging$theta.i
          mean.mean.z <- double.averaging$mean.mean.z
          mean.mean.theta <- double.averaging$mean.mean.theta
          mean.cpt <- double.averaging$mean.cpt
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
      partition.i <- find_startingpoint_single(nodes,sizes.allowed)
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
          results.i <- draw_Metropolis_single(fulltheta.i, partition.i, nodes, effects, objects, burnin, 1, 1, neighborhood, sizes.allowed, sizes.simulated)
        } else {
          results.i <- draw_Metropolis_single(fulltheta.i, partition.i, nodes, effects, objects, thining, 1, 1, neighborhood, sizes.allowed, sizes.simulated)
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

      print(cat("Step",step))
      print(cat("Length of the step",(i-1),"(minimal value:",min.iter[step],"and maximal value:",max.iter[step],")"))
      print(cat("Current estimate",estimates))

    }

    print(cat("Difference to estimated statistics after phase 2, step",step))
    print(colMeans(allz) - z.obs)
    print(cat("Estimates after phase 2, step",step))
    print(estimates)
  }

  return(list(all.estimates = all.estimates,
              final.estimates = estimates,
              lengths.subphases = lengths.subphases))

}

#' Phase 2 wrapper for multiple observation
#'
#' @param partitions XXX
#' @param estimates.phase1 XXX
#' @param inv.zcov XXX
#' @param inv.scaling XXX
#' @param z.obs XXX
#' @param nodes XXX
#' @param effects XXX
#' @param objects XXX
#' @param burnin XXX
#' @param thining XXX
#' @param num.steps XXX
#' @param gainfactors XXX
#' @param r.truncation.p2 XXX
#' @param min.iter XXX
#' @param max.iter XXX
#' @param multiplication.iter XXX
#' @param neighborhood XXX
#' @param fixed.estimates XXX
#' @param sizes.allowed XXX
#' @param sizes.simulated XXX
#' @param double.averaging XXX
#' @return XXX
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
                              sizes.allowed,
                              sizes.simulated,
                              double.averaging) {

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
      #partitions.i <- find_startingpoint_multiple(presence.tables,nodes,sizes.allowed)
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
        mean.cpt <- 0
      }

      # SUB STEP: until generated statistics cross the observed ones
      while(!stop.iterations) {

        print(theta.i)
        all.estimates <- rbind(all.estimates,matrix(theta.i,nrow=1))

        # draw chain
        results.i <- draw_Metropolis_multiple(theta.i, partitions.i, presence.tables, nodes, effects, objects, burnin, 1, 1, neighborhood, sizes.allowed, sizes.simulated)
        #if(i == 1){
        #  results.i <- draw_Metropolis_multiple(theta.i, partitions.i, presence.tables, nodes, effects, objects, burnin, 1, 1, neighborhood, sizes.allowed, sizes.simulated)
        #} else {
        #  results.i <- draw_Metropolis_multiple(theta.i, partitions.i, presence.tables, nodes, effects, objects, thining, 1, 1, neighborhood, sizes.allowed, sizes.simulated)
        #}
        z.i <- results.i$draws
        partition.i <- results.i$last.partitions

        # store z
        allz <- rbind(allz, z.i)

        # compute new parameters
        if(double.averaging) {
          doubleaveraging.step <- compute_parameters_doubleaveraging(z.i, z.obs, theta.i, mean.mean.z, mean.mean.theta, mean.cpt, inv.zcov, inv.scaling, gainfactors[step], r.truncation.p2)
          theta.i <- double.averaging$theta.i
          mean.mean.z <- double.averaging$mean.mean.z
          mean.mean.theta <- double.averaging$mean.mean.theta
          mean.cpt <- double.averaging$mean.cpt
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
      partitions.i <- find_startingpoint_multiple(presence.tables,nodes,sizes.allowed)
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
          results.i <- draw_Metropolis_multiple(theta.i, partitions.i, presence.tables, nodes, effects, objects, burnin, 1, 1, neighborhood, sizes.allowed, sizes.simulated)
        } else {
          results.i <- draw_Metropolis_multiple(theta.i, partitions.i, presence.tables, nodes, effects, objects, thining, 1, 1, neighborhood, sizes.allowed, sizes.simulated)
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

    print(cat("Length of the step",step))
    print(cat((i-1),"(minimal value:",min.iter[step],"and maximal value:",max.iter[step],")"))
    print(cat("Estimated statistics after phase 2, step",step))
    print(colMeans(allz) - z.obs)
    print(cat("Estimates after phase 2, step",step))
    print(estimates)
  }

  return(list(all.estimates = all.estimates,
              final.estimates = estimates,
              lengths.subphases = lengths.subphases))

}


run_phase2_multiple_secondparallel <- function(partitions,
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
                                sizes.allowed,
                                sizes.simulated,
                                double.averaging,
                                parallel,
                                cpus) {

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
      all.estimates <- rbind(all.estimates,matrix(theta.i,nrow=1))

      # find a good starting point
      #partitions.i <- find_startingpoint_multiple(presence.tables,nodes,sizes.allowed)
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
        mean.cpt <- 0
      }

      # SUB STEP: until generated statistics cross the observed ones
      while(!stop.iterations) {

        print(theta.i[1])

        # draw chain
        #if(i == 1){
        #  results.i <- draw_Metropolis_multiple(theta.i, partitions.i, presence.tables, nodes, effects, objects, burnin, 1, 1, neighborhood, sizes.allowed, sizes.simulated)
        #} else {
        #  results.i <- draw_Metropolis_multiple(theta.i, partitions.i, presence.tables, nodes, effects, objects, thining, 1, 1, neighborhood, sizes.allowed, sizes.simulated)
        #}
        results.i <- draw_Metropolis_multiple_secondparallel(theta.i, partitions.i, presence.tables, nodes, effects, objects, burnin, 1, 1, neighborhood, sizes.allowed, sizes.simulated,parallel,cpus)
        z.i <- results.i$draws
        partition.i <- results.i$last.partitions

        # store z
        allz <- rbind(allz, z.i)

        # compute new parameters
        if(double.averaging) {
          doubleaveraging.step <- compute_parameters_doubleaveraging(z.i, z.obs, theta.i, mean.mean.z, mean.mean.theta, mean.cpt, inv.zcov, inv.scaling, gainfactors[step], r.truncation.p2)
          theta.i <- double.averaging$theta.i
          mean.mean.z <- double.averaging$mean.mean.z
          mean.mean.theta <- double.averaging$mean.mean.theta
          mean.cpt <- double.averaging$mean.cpt
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
      partitions.i <- find_startingpoint_multiple(presence.tables,nodes,sizes.allowed)
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
        #if(i == 1){
        #  results.i <- draw_Metropolis_multiple(theta.i, partitions.i, presence.tables, nodes, effects, objects, burnin, 1, 1, neighborhood, sizes.allowed, sizes.simulated)
        #} else {
        #  results.i <- draw_Metropolis_multiple(theta.i, partitions.i, presence.tables, nodes, effects, objects, thining, 1, 1, neighborhood, sizes.allowed, sizes.simulated)
        #}
        results.i <- draw_Metropolis_multiple_secondparallel(theta.i, partitions.i, presence.tables, nodes, effects, objects, burnin, 1, 1, neighborhood, sizes.allowed, sizes.simulated,parallel,cpus)
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

    print(cat("Length of the step",step))
    print(cat((i-1),"(minimal value:",min.iter[step],"and maximal value:",max.iter[step],")"))
    print(cat("Estimated statistics after phase 2, step",step))
    print(colMeans(allz) - z.obs)
    print(cat("Estimates after phase 2, step",step))
    print(estimates)
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
    diff <- t(z.i - z.obs)
    maxratio <-  max(sqrt((t(diff) %*% inv.zcov %*% diff / num.effects)))
    if(maxratio > r.truncation.p2) {
      r <- r.truncation.p2 / maxratio
    }

    # new theta
    theta.i <- theta.i - gainfactor * r * inv.scaling %*% t(z.i - z.obs)

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

  mean.mean.z <- mean.cpt / (mean.cpt+1) *  mean.mean.z + z.i / (mean.cpt+1)

  # compute truncating factor
  r <- 1
  diff <- t(mean.mean.z - z.obs)
  maxratio <-  max(sqrt((t(diff) %*% inv.zcov %*% diff / num.effects)))
  if(maxratio > r.truncation.p2) {
    r <- r.truncation.p2 / maxratio
  }

  # new theta
  theta.i <- theta.i - gainfactor * mean.cpt * r * inv.scaling %*% t(mean.mean.z - z.obs)

  if(mean.cpt == 0) {
    mean.mean.theta <- theta.i
  } else {
    mean.mean.theta <- mean.cpt / (mean.cpt+1) * mean.mean.theta + theta.i / (mean.cpt+1)
  }
  mean.cpt <- mean.cpt + 1

  return(list("theta.i" = theta.i,
              "mean.mean.z" = mean.mean.z,
              "mean.mean.theta" = mean.mean.theta,
              "mean.cpt" = mean.cpt,))

}
