######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Functions used to find appropriate burnin and thining in phase 1 ##
## Author: Marion Hoffman                                           ##
######################################################################

## ------ EXAMINE BURN-IN SEPARATELY ---------

#' Simulate burn in single
#'
#' Function that can be used to find a good length for the burn-in of the Markov chain for a given model and a given set of transitions in the chain (the neighborhood).
#' It draws a chain and calculates the mean statistics for different burn-ins.
#'
#' @param partition A partition (vector)
#' @param theta Initial model parameters
#' @param nodes Node set (data frame)
#' @param effects Effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param objects Objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param num.steps Number of samples wanted
#' @param neighborhood Way of choosing partitions: probability vector (proba actors swap, proba merge/division, proba single actor move)
#' @param numgroups.allowed vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated vector containing the number of groups simulated
#' @param sizes.allowed Vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated Vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @return A list with list the draws, the moving.means and the moving means smoothed
#' @importFrom stats loess
#' @export
simulate_burnin_single <- function(partition, 
                                   theta, 
                                   nodes, 
                                   effects, 
                                   objects, 
                                   num.steps, 
                                   neighborhood, 
                                   numgroups.allowed, 
                                   numgroups.simulated,
                                   sizes.allowed, 
                                   sizes.simulated) 
{
  num.effects <- length(effects$names)

  # print("Neighborhood: ")
  # print(neighborhood)
  
  chain <- draw_Metropolis_single(theta, partition, nodes, effects, objects, 1, 1, num.steps, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, return.all.partitions = FALSE)

  # now we check the evolution of the mean for choosing the burnin
  allmeans <- matrix(0,nrow(chain$draws),num.effects)
  allmeans[1,] <- chain$draws[1,]
  for(d in 2:nrow(chain$draws)) {
    allmeans[d,] <- colMeans(chain$draws[1:d,])
  }

  # we smoothen out the curves (for plots), can be useful when the sample is small
  smoothedmeans <- allmeans
  for(eff in 1:num.effects) {
    lo <- loess(y ~ x, data.frame(x=1:nrow(chain$draws), y=allmeans[,eff]))
    smoothedmeans[!is.na(smoothedmeans[,eff]),eff] <- lo$fitted
  }

  return(list("draws" = chain$draws,
              "moving.means" = allmeans,
              "moving.means.smoothed" = smoothedmeans))

}



#' Grid - search burnin single
#'
#' Function that can be used to find a good length for the burn-in of the Markov chain for a given model and differents sets of transitions in the chain (the neighborhoods).
#' For each neighborhood, it draws a chain and calculates the mean statistics for different burn-ins.
#'
#'
#' @param partition A partition (vector)
#' @param theta Initial model parameters
#' @param nodes Node set (data frame)
#' @param effects Effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param objects Objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param num.steps Number of samples wanted
#' @param neighborhoods List of probability vectors (proba actors swap, proba merge/division, proba single actor move)
#' @param numgroups.allowed = NULL, # vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated vector containing the number of groups simulated
#' @param sizes.allowed Vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated Vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param parallel False, to run different neighborhoods in parallel
#' @param cpus Equal to 1
#' @return all simulations
#' @export
gridsearch_burnin_single <- function(partition, 
                                     theta, 
                                     nodes, 
                                     effects, 
                                     objects, 
                                     num.steps, 
                                     neighborhoods, 
                                     numgroups.allowed, 
                                     numgroups.simulated,
                                     sizes.allowed, 
                                     sizes.simulated, 
                                     parallel = FALSE, 
                                     cpus = 1) {

  # parallel procedure
  if(parallel){

    #distribute neighborhoods in each cpu
    n <- ceiling(length(neighborhoods) / cpus)
    subindexes <- list()
    for(c in 1:cpus){
      start <- (c-1)*n + 1
      end <- c*n
      if(c == cpus) end <- length(neighborhoods)
      subindexes[[c]] <- start:end
    }

    sfExport("partition", "theta", "nodes", "effects", "objects", "num.steps", "neighborhoods", "numgroups.allowed", "numgroups.simulated", "sizes.allowed", "sizes.simulated", "subindexes")
    res <- sfLapply(1:cpus, fun = function(k) {
      subres <- list()
      for(i in 1:length(subindexes[[k]])){
        index <- subindexes[[k]][i]
        subneighborhood <- neighborhoods[[index]]
        subres[[i]] <- simulate_burnin_single(partition, theta, nodes, effects, objects, num.steps, subneighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
      }
      return(subres)
    }
    )

    # append all results
    allsimulations <- list()
    for(c in 1:cpus) allsimulations <- append(allsimulations, res[[c]])

  } else {

    # just go through all neighborhoods one by one
    allsimulations <- list()
    for(i in 1:length(neighborhoods)){
      allsimulations[[i]] <- simulate_burnin_single(partition, theta, nodes, effects, objects, num.steps, neighborhoods[[i]], numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated)
    }

  }

  return(list("all.simulations" = allsimulations))

}


## ------ EXAMINE THINING SEPARATELY ---------

#' Simulate thining single
#'
#' Function that can be used to find a good length for the thining of the Markov chain for a given model and a set of transitions in the chain (the neighborhood).
#' It draws a chain and calculates the autocorrelation of statistics for different thinings.
#'
#' @param partition A partition (vector)
#' @param theta Initial model parameters
#' @param nodes Node set (data frame)
#' @param effects Effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param objects Objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param num.steps Number of samples wanted
#' @param neighborhood Way of choosing partitions: probability vector (proba actors swap, proba merge/division, proba single actor move)
#' @param numgroups.allowed vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated vector containing the number of groups simulated
#' @param sizes.allowed Vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated Vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param burnin number of simulated steps for the burn-in
#' @param max.thining maximal number of simulated steps in the thining
#' @param verbose logical: should intermediate results during the estimation be printed or not? Defaults to FALSE.
#' @return A list
#' @importFrom stats cor loess
#' @export
# SINGLE PARTITION PROCEDURE
simulate_thining_single <- function(partition, 
                                    theta, 
                                    nodes, 
                                    effects, 
                                    objects,
                                    num.steps, 
                                    neighborhood, 
                                    numgroups.allowed,
                                    numgroups.simulated,
                                    sizes.allowed, 
                                    sizes.simulated, 
                                    burnin,
                                    max.thining,
                                    verbose = FALSE)
{
  num.effects <- length(effects$names)

  chain <- draw_Metropolis_single(theta, partition, nodes, effects, objects, burnin, 1, num.steps, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, return.all.partitions = FALSE)

  # first look at autocorrelations find ok thining (either max thining reached or autocor < 0.4)
  ok <- FALSE
  current.draws <- chain$draws
  current.last <- chain$last.partition
  current.thining <- 1

  allautocors <- c()

  while(!ok){

    # we check the auto for the current burnin
    draws <- current.draws[seq(1, nrow(current.draws), current.thining),]
    autocors <- rep(0,num.effects)
    for(e in 1:num.effects){
      autocors[e] <- cor(draws[1:(num.steps-1),e],draws[2:num.steps,e])
    }
    if(current.thining %% 50 == 0){
      if (verbose) {
        cat("thining\n")
        cat(current.thining, "\n\n")
        cat("autocorrelations\n")
        cat(autocors, "\n\n")
      }
    }

    allautocors <- rbind(allautocors,autocors)

    #if(max(autocors,na.rm = TRUE) < 0.4) {
    if(current.thining >= max.thining ) {

      # we can stop
      ok <- TRUE

    } else{

      # we continue the chain
      new.chain <- draw_Metropolis_single(theta, current.last, nodes, effects, objects, 1, 1, num.steps, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, return.all.partitions = FALSE)
      current.draws <- rbind(current.draws,new.chain$draws)
      current.last <- new.chain$last.partition

      current.thining <- current.thining + 1
    }


  }

  # we smoothen out the curves (for plots), can be useful when the sample is small
  smoothedautocors <- allautocors
  for(eff in 1:num.effects) {
    lo <- loess(y ~ x, data.frame(x=1:max.thining, y=allautocors[,eff]))
    smoothedautocors[!is.na(smoothedautocors[,eff]),eff] <- lo$fitted
  }


  return(list("draws" = current.draws,
              "autocorrelations" = allautocors,
              "autocorrelations.smoothed" = smoothedautocors))

}



#' Grid - search thining single
#'
#' Function that can be used to find a good length for the thining of the Markov chain for a given model and differents sets of transitions in the chain (the neighborhoods).
#' For each neighborhood, it draws a chain and calculates the autocorrelation of statistics for different thinings.
#'
#' @param partition A partition (vector)
#' @param theta Initial model parameters
#' @param nodes Node set (data frame)
#' @param effects Effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param objects Objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param num.steps Number of samples wanted
#' @param neighborhoods List of probability vectors (proba actors swap, proba merge/division, proba single actor move)
#' @param numgroups.allowed vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated vector containing the number of groups simulated
#' @param sizes.allowed Vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated Vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param burnin length of the burn-in period
#' @param max.thining maximal value for the thining to be tested
#' @param parallel False, to run different neighborhoods in parallel
#' @param cpus Equal to 1
#' @return all simulations
#' @export
gridsearch_thining_single <- function(partition,
                                      theta, 
                                      nodes, 
                                      effects, 
                                      objects, 
                                      num.steps, 
                                      neighborhoods, 
                                      numgroups.allowed, 
                                      numgroups.simulated,
                                      sizes.allowed, 
                                      sizes.simulated, 
                                      burnin, 
                                      max.thining, 
                                      parallel = FALSE, 
                                      cpus = 1) {

  # parallel procedure
  if(parallel){

    #distribute neighborhoods in each cpu
    n <- ceiling(length(neighborhoods) / cpus)
    subindexes <- list()
    for(c in 1:cpus){
      start <- (c-1)*n + 1
      end <- c*n
      if(c == cpus) end <- length(neighborhoods)
      subindexes[[c]] <- start:end
    }

    sfExport("partition", "theta", "nodes", "effects", "objects", "num.steps", "neighborhoods", "numgroups.allowed", "numgroups.simulated", "sizes.allowed", "sizes.simulated", "burnin", "max.thining", "subindexes")
    res <- sfLapply(1:cpus, fun = function(k) {
      subres <- list()
      for(i in 1:length(subindexes[[k]])){
        index <- subindexes[[k]][i]
        subneighborhood <- neighborhoods[[index]]
        subres[[i]] <- simulate_thining_single(partition, theta, nodes, effects, objects, num.steps, subneighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, burnin, max.thining)
      }
      return(subres)
    }
    )

    # append all results
    allsimulations <- list()
    for(c in 1:cpus) allsimulations <- append(allsimulations, res[[c]])

  } else {

    # just go through all neighborhoods one by one
    allsimulations <- list()
    for(i in 1:length(neighborhoods)){
      allsimulations[[i]] <- simulate_thining_single(partition, theta, nodes, effects, objects, num.steps, neighborhoods[[i]], numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, burnin, max.thining)
    }

  }

  # collect properties for all neighborhoods
  p1 <- c()
  p2 <- c()
  p3 <- c()
  reach0.4 <- c()
  firstthining0.4 <- c()
  reach0.5 <- c()
  firstthining0.5 <- c()
  finalmaxautocor <- c()
  for(i in 1:length(neighborhoods)){
    p1 <- c(p1, neighborhoods[[i]][1])
    p2 <- c(p2, neighborhoods[[i]][2])
    p3 <- c(p3, neighborhoods[[i]][3])

    allmaxautocors <- apply(allsimulations[[i]]$autocorrelations.smoothed, 1, FUN=max)
    finalmaxautocor <- c(finalmaxautocor, allmaxautocors[max.thining])

    reach0.4 <- c(reach0.4, sum(allmaxautocors <= 0.4, na.rm = TRUE))
    reach0.5 <- c(reach0.5, sum(allmaxautocors <= 0.5, na.rm = TRUE))

    if(length(which(allmaxautocors <= 0.4)) > 0) {
      firstthining0.4 <- c(firstthining0.4, min(which(allmaxautocors <= 0.4)))
    } else {
      firstthining0.4 <- c(firstthining0.4, -1)
    }
    if(length(which(allmaxautocors <= 0.5)) > 0) {
      firstthining0.5 <- c(firstthining0.5, min(which(allmaxautocors <= 0.5)))
    } else {
      firstthining0.5 <- c(firstthining0.5, -1)
    }
  }
  results.search <- data.frame(p_swap = p1,
                               p_mergediv = p2,
                               p_single = p3,
                               final_max_autocorr = finalmaxautocor,
                               reach_autocorr0.4 = reach0.4,
                               first_thining_autocorr0.4 = firstthining0.4,
                               reach_autocorr0.5 = reach0.5,
                               first_thining_autocorr0.5 = firstthining0.5)

  return(list("results.search" = results.search,
              "all.simulations" = allsimulations))

}



## ------ EXAMINE BURN-IN AND THINING TOGETHER ---------


#' Simulate burnin thining single
#'
#' Function that simulates the Markov chain for a given model and a set of transitions (the neighborhood), for a single partition.
#' It calculates the autocorrelation of statistics for different thinings and the average statistics for different burn-ins.
#'
#' @param partition Observed partition (vector)
#' @param theta Initial model parameters
#' @param nodes Node set (data frame)
#' @param effects Effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param objects Objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param num.steps Number of samples wanted
#' @param neighborhood Way of choosing partitions: probability vector (proba actors swap, proba merge/division, proba single actor move)
#' @param numgroups.allowed vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated vector containing the number of groups simulated
#' @param sizes.allowed Vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated Vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param max.thining maximal number of simulated steps in the thining
#' @param verbose logical: should intermediate results during the estimation be printed or not? Defaults to FALSE.
#' @return A list
#' @importFrom stats loess
#' @export
simulate_burninthining_single <- function(partition, 
                                          theta,
                                          nodes,
                                          effects,
                                          objects,
                                          num.steps,
                                          neighborhood,
                                          numgroups.allowed,
                                          numgroups.simulated,
                                          sizes.allowed, 
                                          sizes.simulated, 
                                          max.thining,
                                          verbose = FALSE)
{
  num.effects <- length(effects$names)

  chain <- draw_Metropolis_single(theta, partition, nodes, effects, objects, 1, 1, num.steps, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, return.all.partitions = FALSE)

  # first look at autocorrelations find ok thining (either max thining reached or autocor < 0.4)
  ok <- FALSE
  current.draws <- chain$draws
  current.last <- chain$last.partition
  current.thining <- 1

  allautocors <- c()

  while(!ok){

    # we check the auto for the current burnin
    draws <- current.draws[seq(1, nrow(current.draws), current.thining),]
    autocors <- rep(0,num.effects)
    for(e in 1:num.effects){
      autocors[e] <- cor(draws[1:(num.steps-1),e],draws[2:num.steps,e])
    }
    if(current.thining %% 50 == 0){
      if (verbose) {
        cat("thining\n")
        cat(current.thining, "\n\n")
        cat("autocorrelations\n")
        cat(autocors, "\n\n")
      }
    }

    allautocors <- rbind(allautocors,autocors)

    #if(max(autocors,na.rm = TRUE) < 0.4) {
    if(current.thining >= max.thining ) {

      # we can stop
      ok <- TRUE

    } else{

      # we continue the chain
      new.chain <- draw_Metropolis_single(theta, current.last, nodes, effects, objects, 1, 1, num.steps, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, return.all.partitions = FALSE)
      current.draws <- rbind(current.draws,new.chain$draws)
      current.last <- new.chain$last.partition

      current.thining <- current.thining + 1
    }


  }

  # now we check the evolution of the mean for choosing the burnin
  allmeans <- matrix(0,nrow(draws),num.effects)
  allmeans[1,] <- draws[1,]
  for(d in 2:nrow(draws)) {
    allmeans[d,] <- colMeans(draws[1:d,])
  }

  # we smoothen out the curves (for plots), can be useful when the sample is small
  smoothedautocors <- allautocors
  smoothedmeans <- allmeans
  for(eff in 1:num.effects) {
    lo <- loess(y ~ x, data.frame(x=1:max.thining, y=allautocors[,eff]))
    smoothedautocors[!is.na(smoothedautocors[,eff]),eff] <- lo$fitted
    lo <- loess(y ~ x, data.frame(x=1:nrow(draws), y=allmeans[,eff]))
    smoothedmeans[!is.na(smoothedmeans[,eff]),eff] <- lo$fitted
  }


  return(list("draws" = current.draws,
              "autocorrelations" = allautocors,
              "moving.means" = allmeans,
              "autocorrelations.smoothed" = smoothedautocors,
              "moving.means.smoothed" = smoothedmeans))

}


# MULTIPLE PARTITION PROCEDURE


#' Simulate burnin thining multiple
#'
#' Function that simulates the Markov chain for a given model and a set of transitions (the neighborhood), for multiple partitions.
#' It calculates the autocorrelation of statistics for different thinings and the average statistics for different burn-ins.
#'
#' @param partitions Observed partitions
#' @param theta Initial model parameters
#' @param presence.tables to indicate which nodes were present when
#' @param nodes Node set (data frame)
#' @param effects Effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param objects Objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param num.steps Number of samples wanted
#' @param neighborhood Way of choosing partitions: probability vector (proba actors swap, proba merge/division, proba single actor move)
#' @param numgroups.allowed vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated vector containing the number of groups simulated
#' @param sizes.allowed Vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated Vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param max.thining maximal number of simulated steps in the thining
#' @param verbose logical: should intermediate results during the estimation be printed or not? Defaults to FALSE.
#' @return A list
#' @importFrom stats cor loess
#' @export
simulate_burninthining_multiple <- function(partitions, 
                                            presence.tables,
                                            theta,
                                            nodes, 
                                            effects,
                                            objects, 
                                            num.steps,
                                            neighborhood, 
                                            numgroups.allowed,
                                            numgroups.simulated, 
                                            sizes.allowed,
                                            sizes.simulated,
                                            max.thining,
                                            verbose = FALSE)
{
  num.effects <- length(effects$names)

  chain <- draw_Metropolis_multiple(theta, partitions, presence.tables, nodes, effects, objects, 1, 1, num.steps, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, return.all.partitions = FALSE)

  # first look at autocorrelations find ok thining (either max thining reached or autocor < 0.4)
  ok <- FALSE
  current.draws <- chain$draws
  current.last <- chain$last.partitions
  current.thining <- 1

  allautocors <- c()

  while(!ok){

    # we check the auto for the current burnin
    draws <- current.draws[seq(1, nrow(current.draws), current.thining),]
    autocors <- rep(0,num.effects)
    for(e in 1:num.effects){
      autocors[e] <- cor(draws[1:(num.steps-1),e],draws[2:num.steps,e])
    }
    if(current.thining %% 50 == 0){
      if (verbose) {
        cat("thining\n")
        cat(current.thining, "\n\n")
        cat("autocorrelations\n")
        cat(autocors, "\n\n")
      }
    }

    allautocors <- rbind(allautocors,autocors)

    #if(max(autocors,na.rm = TRUE) < 0.4) {
    if(current.thining >= max.thining ) {

      # we can stop
      ok <- TRUE

    } else{

      # we continue the chain
      new.chain <- draw_Metropolis_multiple(theta, current.last, presence.tables, nodes, effects, objects, 1, 1, num.steps, neighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, return.all.partitions = FALSE)
      current.draws <- rbind(current.draws,new.chain$draws)
      current.last <- new.chain$last.partitions

      current.thining <- current.thining + 1
    }


  }

  # now we check the evolution of the mean for choosing the burnin
  allmeans <- matrix(0,nrow(draws),num.effects)
  allmeans[1,] <- draws[1,]
  for(d in 2:nrow(draws)) {
    allmeans[d,] <- colMeans(draws[1:d,])
  }

  # we smoothen out the curves (for plots), can be useful when the sample is small
  smoothedautocors <- allautocors
  smoothedmeans <- allmeans
  for(eff in 1:num.effects) {
    lo <- loess(y ~ x, data.frame(x=1:max.thining, y=allautocors[,eff]))
    smoothedautocors[!is.na(smoothedautocors[,eff]),eff] <- lo$fitted
    lo <- loess(y ~ x, data.frame(x=1:nrow(draws), y=allmeans[,eff]))
    smoothedmeans[!is.na(smoothedmeans[,eff]),eff] <- lo$fitted
  }


  return(list("draws"= current.draws,
              "autocorrelations" = allautocors,
              "moving.means" = allmeans,
              "autocorrelations.smoothed" = smoothedautocors,
              "moving.means.smoothed" = smoothedmeans))

}

# GRID SEARCH FOR SINGLE PARTITION

#' Grid - search burnin thining single
#'
#' Function that simulates the Markov chain for a given model and several sets of transitions (the neighborhoods), for a single partition.
#' For each neighborhood, it calculates the autocorrelation of statistics for different thinings and the average statistics for different burn-ins.
#' Then the best neighborhood can be selected along with good values for burn-in and thining
#'
#' @param partition A partition (vector)
#' @param theta Initial model parameters
#' @param nodes Node set (data frame)
#' @param effects Effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param objects Objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param num.steps Number of samples wanted
#' @param neighborhoods List of probability vectors (proba actors swap, proba merge/division, proba single actor move)
#' @param numgroups.allowed vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated vector containing the number of groups simulated
#' @param sizes.allowed Vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated Vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param max.thining Where to stop adding thining
#' @param parallel False, to run different neighborhoods in parallel
#' @param cpus Equal to 1
#' @return list
#' @export
gridsearch_burninthining_single <- function(partition,
                                            theta,
                                            nodes,
                                            effects, 
                                            objects, 
                                            num.steps, 
                                            neighborhoods,
                                            numgroups.allowed, 
                                            numgroups.simulated,
                                            sizes.allowed, 
                                            sizes.simulated, 
                                            max.thining, 
                                            parallel = FALSE, 
                                            cpus = 1) {

  # parallel procedure
  if(parallel){

    #distribute neighborhoods in each cpu
    n <- ceiling(length(neighborhoods) / cpus)
    subindexes <- list()
    for(c in 1:cpus){
      start <- (c-1)*n + 1
      end <- c*n
      if(c == cpus) end <- length(neighborhoods)
      subindexes[[c]] <- start:end
    }

    sfExport("partition", "theta", "nodes", "effects", "objects", "num.steps", "neighborhoods", "numgroups.allowed", "numgroups.simulated", "sizes.allowed", "sizes.simulated", "max.thining", "subindexes")
    res <- sfLapply(1:cpus, fun = function(k) {
      subres <- list()
      for(i in 1:length(subindexes[[k]])){
        index <- subindexes[[k]][i]
        subneighborhood <- neighborhoods[[index]]
        subres[[i]] <- simulate_burninthining_single(partition, theta, nodes, effects, objects, num.steps, subneighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, max.thining)
      }
      return(subres)
    }
    )

    # append all results
    allsimulations <- list()
    for(c in 1:cpus) allsimulations <- append(allsimulations, res[[c]])

  } else {

    # just go through all neighborhoods one by one
    allsimulations <- list()
    for(i in 1:length(neighborhoods)){
      allsimulations[[i]] <- simulate_burninthining_single(partition, theta, nodes, effects, objects, num.steps, neighborhoods[[i]], numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, max.thining)
    }

  }

  # collect properties for all neighborhoods
  p1 <- c()
  p2 <- c()
  p3 <- c()
  reach0.4 <- c()
  firstthining0.4 <- c()
  reach0.5 <- c()
  firstthining0.5 <- c()
  finalmaxautocor <- c()
  for(i in 1:length(neighborhoods)){
    p1 <- c(p1, neighborhoods[[i]][1])
    p2 <- c(p2, neighborhoods[[i]][2])
    p3 <- c(p3, neighborhoods[[i]][3])

    allmaxautocors <- apply(allsimulations[[i]]$autocorrelations.smoothed, 1, FUN=max)
    finalmaxautocor <- c(finalmaxautocor, allmaxautocors[max.thining])

    reach0.4 <- c(reach0.4, sum(allmaxautocors <= 0.4, na.rm = TRUE))
    reach0.5 <- c(reach0.5, sum(allmaxautocors <= 0.5, na.rm = TRUE))

    if(length(which(allmaxautocors <= 0.4)) > 0) {
      firstthining0.4 <- c(firstthining0.4, min(which(allmaxautocors <= 0.4)))
    } else {
      firstthining0.4 <- c(firstthining0.4, -1)
    }
    if(length(which(allmaxautocors <= 0.5)) > 0) {
      firstthining0.5 <- c(firstthining0.5, min(which(allmaxautocors <= 0.5)))
    } else {
      firstthining0.5 <- c(firstthining0.5, -1)
    }
  }
  results.search <- data.frame(p_swap = p1,
                               p_mergediv = p2,
                               p_single = p3,
                               final_max_autocorr = finalmaxautocor,
                               reach_autocorr0.4 = reach0.4,
                               first_thining_autocorr0.4 = firstthining0.4,
                               reach_autocorr0.5 = reach0.5,
                               first_thining_autocorr0.5 = firstthining0.5)

  return(list("results.search" = results.search,
              "all.simulations" = allsimulations))

}


# GRID SEARCH FOR MULTIPLE PARTITIONS

#' Grid - search burnin thining multiple
#'
#' Function that simulates the Markov chain for a given model and several sets of transitions (the neighborhoods), for multiple partitions.
#' For each neighborhood, it calculates the autocorrelation of statistics for different thinings and the average statistics for different burn-ins.
#' Then the best neighborhood can be selected along with good values for burn-in and thining
#'
#' @param partitions Observed partitions
#' @param presence.tables Presence of nodes
#' @param theta Initial model parameters
#' @param nodes Node set (data frame)
#' @param effects Effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param objects Objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param num.steps Number of samples wanted
#' @param neighborhoods List of probability vectors (proba actors swap, proba merge/division, proba single actor move)
#' @param numgroups.allowed vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated vector containing the number of groups simulated
#' @param sizes.allowed Vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated Vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param max.thining Where to stop adding thining
#' @param parallel False, to run different neighborhoods in parallel
#' @param cpus Equal to 1
#' @return list
#' @importFrom snowfall sfExport sfLapply
#' @export
gridsearch_burninthining_multiple <- function(partitions, 
                                              presence.tables, 
                                              theta, 
                                              nodes, 
                                              effects, 
                                              objects, 
                                              num.steps, 
                                              neighborhoods,
                                              numgroups.allowed,
                                              numgroups.simulated,
                                              sizes.allowed, 
                                              sizes.simulated,
                                              max.thining, 
                                              parallel = FALSE, 
                                              cpus = 1) {

  # parallel procedure
  if(parallel){

    #distribute neighborhoods in each cpu
    n <- ceiling(length(neighborhoods) / cpus)
    subindexes <- list()
    for(c in 1:cpus){
      start <- (c-1)*n + 1
      end <- c*n
      if(c == cpus) end <- length(neighborhoods)
      subindexes[[c]] <- start:end
    }

    sfExport("partitions", "presence.tables", "theta", "nodes", "effects", "objects", "num.steps", "neighborhoods", "numgroups.allowed", "numgroups.simulated", "sizes.allowed", "sizes.simulated", "max.thining", "subindexes")
    res <- sfLapply(1:cpus, fun = function(k) {
      subres <- list()
      for(i in 1:length(subindexes[[k]])){
        index <- subindexes[[k]][i]
        subneighborhood <- neighborhoods[[index]]
        subres[[i]] <- simulate_burninthining_multiple(partitions, presence.tables, theta, nodes, effects, objects, num.steps, subneighborhood, numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, max.thining)
      }
      return(subres)
    }
    )

    # append all results
    allsimulations <- list()
    for(c in 1:cpus) allsimulations <- append(allsimulations, res[[c]])

  } else {

    # just go through all neighborhoods one by one
    allsimulations <- list()
    for(i in 1:length(neighborhoods)){
      # print(neighborhoods[[i]])
      allsimulations[[i]] <- simulate_burninthining_multiple(partitions, presence.tables, theta, nodes, effects, objects, num.steps, neighborhoods[[i]], numgroups.allowed, numgroups.simulated, sizes.allowed, sizes.simulated, max.thining)
    }

  }

  # collect properties for all neighborhoods
  p1 <- c()
  p2 <- c()
  p3 <- c()
  reach0.4 <- c()
  firstthining0.4 <- c()
  reach0.5 <- c()
  firstthining0.5 <- c()
  finalmaxautocor <- c()
  for(i in 1:length(neighborhoods)){
    p1 <- c(p1, neighborhoods[[i]][1])
    p2 <- c(p2, neighborhoods[[i]][2])
    p3 <- c(p3, neighborhoods[[i]][3])

    allmaxautocors <- apply(allsimulations[[i]]$autocorrelations.smoothed, 1, FUN=max)
    finalmaxautocor <- c(finalmaxautocor, allmaxautocors[max.thining])

    reach0.4 <- c(reach0.4, sum(allmaxautocors <= 0.4, na.rm = TRUE))
    reach0.5 <- c(reach0.5, sum(allmaxautocors <= 0.5, na.rm = TRUE))

    if(length(which(allmaxautocors <= 0.4)) > 0) {
      firstthining0.4 <- c(firstthining0.4, min(which(allmaxautocors <= 0.4)))
    } else {
      firstthining0.4 <- c(firstthining0.4, -1)
    }
    if(length(which(allmaxautocors <= 0.5)) > 0) {
      firstthining0.5 <- c(firstthining0.5, min(which(allmaxautocors <= 0.5)))
    } else {
      firstthining0.5 <- c(firstthining0.5, -1)
    }
  }
  results.search <- data.frame(p_swap = p1,
                               p_mergediv = p2,
                               p_single = p3,
                               final_max_autocorr = finalmaxautocor,
                               reach_autocorr0.4 = reach0.4,
                               first_thining_autocorr0.4 = firstthining0.4,
                               reach_autocorr0.5 = reach0.5,
                               first_thining_autocorr0.5 = firstthining0.5)

  return(list("results.search" = results.search,
              "all.simulations" = allsimulations))

}
