######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Function implementing the Metropolis HAstings algorithm to       ##
## through partitions given a certain model specification           ##
## Author: Marion Hoffman                                           ##
######################################################################


## --- FUNCTIONS TO DRAW CHAINS ----

#' Draw Metropolis single
#'
#' Function to sample the model with a Markov chain (single partition procedure).
#'
#'
#' @param theta model parameters
#' @param first.partition, starting partition for the Markov chain
#' @param nodes nodeset (data frame)
#' @param effects effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param objects objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param burnin integer for the number of burn-in steps before sampling
#' @param thining integer for the number of thining steps between sampling
#' @param num.steps number of samples
#' @param neighborhood = c(0.7, 0.3, 0), way of choosing partitions: probability vector (2 actors swap, merge/division, single actor move, single pair move, 2 pairs swap, 2 groups reshuffle)
#' @param numgroups.allowed = NULL, # vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated = NULL, # vector containing the number of groups simulated
#' @param sizes.allowed = NULL,  vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated = NULL, vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param return.all.partitions = FALSE option to return the sampled partitions on top of their statistics (for GOF)
#' @return A list
#' @importFrom stats runif
#' @export
#' @examples
#' # define an arbitrary set of n = 6 nodes with attributes, and an arbitrary covariate matrix
#' n <- 6
#' nodes <- data.frame(label = c("A","B","C","D","E","F"),
#'                     gender = c(1, 1, 2, 1, 2, 2),
#'                     age = c(20, 22, 25, 30, 30, 31))
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
#' # set parameter values for each of these effects
#' parameters <- c(-0.2, 0.2, -0.1, 0.5)
#'
#' \donttest{
#' # generate simulated sample, by setting the desired additional parameters for the
#' # Metropolis sampler and choosing a starting point for the chain (first.partition)
#' nsteps <- 100
#' sample <- draw_Metropolis_single(theta = parameters,
#'                                  first.partition = c(1, 1, 2, 2, 3, 3),
#'                                  nodes = nodes,
#'                                  effects = effects,
#'                                  objects = objects,
#'                                  burnin = 100,
#'                                  thining = 10,
#'                                  num.steps = nsteps,
#'                                  neighborhood = c(0, 1, 0),
#'                                  numgroups.allowed = 1:n,
#'                                  numgroups.simulated = 1:n,
#'                                  sizes.allowed = 1:n,
#'                                  sizes.simulated = 1:n,
#'                                  return.all.partitions = TRUE)
#'
#'
#' # or: simulate an estimated model
#' partition <- c(1, 1, 2, 2, 2, 3) # the partition already defined for the (previous) estimation
#' nsimulations <- 1000
#' simulations <- draw_Metropolis_single(theta = estimation$results$est,
#'                                       first.partition = partition,
#'                                       nodes = nodes,
#'                                       effects = effects,
#'                                       objects = objects,
#'                                       burnin = 100,
#'                                       thining = 20,
#'                                       num.steps = nsimulations,
#'                                       neighborhood = c(0, 1, 0),
#'                                       sizes.allowed = 1:n,
#'                                       sizes.simulated = 1:n,
#'                                       return.all.partitions = TRUE)
#' }
#'
draw_Metropolis_single <- function(theta,
                                   first.partition,
                                   nodes,
                                   effects,
                                   objects,
                                   burnin,
                                   thining,
                                   num.steps,
                                   neighborhood = c(0.7, 0.3, 0),
                                   numgroups.allowed = NULL,
                                   numgroups.simulated = NULL,
                                   sizes.allowed = NULL,
                                   sizes.simulated = NULL,
                                   return.all.partitions = FALSE) {

  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)

  # turn the neighborhood weights into probabilities if needed
  if (sum(neighborhood) != 1) {
    neighborhood <- neighborhood / sum(neighborhood)
  }

  # check whether there are constraints
  constraints <- FALSE
  if (!is.null(numgroups.allowed) || !is.null(sizes.allowed)) {
    constraints <- TRUE
    if (is.null(numgroups.allowed)) numgroups.allowed <- 1:nrow(nodes)
    if (is.null(numgroups.simulated)) numgroups.simulated <- numgroups.allowed
    if (is.null(sizes.allowed)) sizes.allowed <- 1:nrow(nodes)
    if (is.null(sizes.simulated)) sizes.simulated <- sizes.allowed
  }

  # instantiate with the starting network
  current.partition <- first.partition
  current.z <- computeStatistics(current.partition, nodes, effects, objects)
  current.logit <- theta * current.z

  # store the statistics collected for all networks simulated after the burn in
  all.z <- c()

  # store the partitions if needed (for GOF)
  if (return.all.partitions) {
    all.partitions <- c()
  }

  end.walk <- FALSE
  cpt_burnin <- 0
  cpt_thining <- 0
  cpt_steps <- 0

  while (!end.walk) {

    new.step <- draw_step_single(theta,
                                 current.partition,
                                 current.logit,
                                 current.z,
                                 nodes,
                                 effects,
                                 objects,
                                 neighborhood,
                                 numgroups.simulated,
                                 sizes.simulated)

    proba.change <- min(1, new.step$hastings.ratio)
    change.made <- (runif(1, 0, 1) <= proba.change)

    # update partition if needed
    old.partition <- current.partition
    if (change.made) {
      current.partition <- new.step$new.partition
      current.z <- new.step$new.z
      current.logit <- new.step$new.logit
    }

    cpt_burnin <- cpt_burnin + 1
    if (cpt_burnin > burnin) cpt_thining <- cpt_thining + 1


    # store the results if we are out of burnin
    if (cpt_burnin >= burnin && cpt_thining == thining) {

      store <- TRUE
      if (constraints) store <- check_sizes(current.partition, sizes.allowed, numgroups.allowed)

      if (store) {
        all.z <- rbind(all.z, current.z)
        cpt_thining <- 0
        cpt_steps <- cpt_steps + 1
        if (return.all.partitions) { # store all partitions if needed
          all.partitions <- rbind(all.partitions, current.partition)
        }
      } else {
        cpt_thining <- thining - 1
      }

    }

    # stop the walk if number of steps reached
    end.walk <- (cpt_steps >= num.steps)

  }

  # TODO: change this hack
  row.names(all.z) <- rep("", dim(all.z)[1])

  # compute the average statistics and the final network generated
  if (!return.all.partitions) {
    return(list("draws" = all.z,
                "last.partition" = current.partition,
                "all.partitions" = NULL))
  } else {
    return(list("draws" = all.z,
                "last.partition" = current.partition,
                "all.partitions" = all.partitions))
  }

}


# MULTIPLE PARTITIONS PROCEDURE

#' Draw Metropolis multiple
#'
#' Function to sample the model with a Markov chain (single partition procedure).
#'
#' @param theta model parameters
#' @param first.partitions starting partition for the Markov chain
#' @param presence.tables  matrix indicating which actors were present for each observations (mandatory)
#' @param nodes node set (data frame)
#' @param effects effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param objects objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param burnin integer for the number of burn-in steps before sampling
#' @param thining integer for the number of thining steps between sampling
#' @param num.steps number of samples
#' @param neighborhood = c(0.7, 0.3, 0), way of choosing partitions: probability vector (2 actors swap, merge/division, single actor move, single pair move, 2 pairs swap, 2 groups reshuffle)
#' @param numgroups.allowed = NULL, # vector containing the number of groups allowed in the partition (now, it only works with vectors like num_min:num_max)
#' @param numgroups.simulated = NULL, # vector containing the number of groups simulated
#' @param sizes.allowed = NULL,  vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
#' @param sizes.simulated = NULL, vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max)
#' @param return.all.partitions = FALSE, option to return the sampled partitions on top of their statistics (for GOF)
#' @param verbose logical: should intermediate results during the estimation be printed or not? Defaults to FALSE.
#' @return A list
#' @importFrom stats runif
#' @export
#' @examples
#' # define an arbitrary set of n = 6 nodes with attributes, and an arbitrary covariate matrix
#' n <- 6
#' nodes <- data.frame(label = c("A","B","C","D","E","F"),
#'                     gender = c(1, 1, 2, 1, 2, 2),
#'                     age = c(20, 22, 25, 30, 30, 31))
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
#' effects_multiple <- list(names = c("num_groups","same","diff","tie","inertia_1"),
#'                 objects = c("partitions","gender","age","friendship","partitions"),
#'                 objects2 = c("","","","",""))
#' objects_multiple <- list()
#' objects_multiple[[1]] <- list(name = "friendship", object = friendship)
#'
#' # set parameter values for each of these effects
#' parameters <- c(-0.2, 0.2,-0.1, 0.5, 1)
#'
#' # set a starting point for the simulation
#' first.partitions <- matrix(c(1, 1, 2, 2, 2, 3,
#'                              NA, 1, 1, 2, 2, 2,
#'                              1, NA, 2, 3, 3, 1), 6, 3)
#'
#' \donttest{
#' # generate the simulated sample
#' nsteps <- 50
#' sample <- draw_Metropolis_multiple(theta = parameters,
#'                                    first.partitions = first.partitions,
#'                                    nodes = nodes,
#'                                    presence.tables = presence.tables,
#'                                    effects = effects_multiple,
#'                                    objects = objects_multiple,
#'                                    burnin = 100,
#'                                    thining = 100,
#'                                    num.steps = nsteps,
#'                                    neighborhood = c(0, 1, 0),
#'                                    numgroups.allowed = 1:n,
#'                                    numgroups.simulated = 1:n,
#'                                    sizes.allowed = 1:n,
#'                                    sizes.simulated = 1:n,
#'                                    return.all.partitions = TRUE)
#' }
#'
draw_Metropolis_multiple <- function(theta,
                                     first.partitions,
                                     presence.tables,
                                     nodes,
                                     effects,
                                     objects,
                                     burnin,
                                     thining,
                                     num.steps,
                                     neighborhood = c(0.7, 0.3, 0),
                                     numgroups.allowed,
                                     numgroups.simulated,
                                     sizes.allowed,
                                     sizes.simulated,
                                     return.all.partitions = FALSE,
                                     verbose = FALSE) {

  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)
  num.obs <- ncol(presence.tables)

  # turn the neighborhood weights into probabilities if needed
  if (sum(neighborhood) != 1) {
    neighborhood <- neighborhood / sum(neighborhood)
  }

  # check whether there are constraints
  constraints <- FALSE
  if (!is.null(numgroups.allowed) || !is.null(sizes.allowed)) {
    constraints <- TRUE
    if (is.null(numgroups.allowed)) numgroups.allowed <- 1:nrow(nodes)
    if (is.null(numgroups.simulated)) numgroups.simulated <- numgroups.allowed
    if (is.null(sizes.allowed)) sizes.allowed <- 1:nrow(nodes)
    if (is.null(sizes.simulated)) sizes.simulated <- sizes.allowed
  }

  # instantiate with the starting network
  current.partitions <- first.partitions
  current.z.contributions <- computeStatistics_multiple(current.partitions, presence.tables, nodes, effects, objects)
  current.z <- rowSums(current.z.contributions)
  current.logit <- theta * current.z

  # store the statistics collected for all networks simulated after the burn in
  all.z <- c()

  # store the partitions if needed (for GOF)
  if (return.all.partitions) {
    all.partitions <- list()
    cpt_all <- 0
  }

  end.walk <- FALSE
  cpt_burnin <- 0
  cpt_thining <- 0
  cpt_steps <- 0

  while (!end.walk) {

    new.step <- draw_step_multiple(theta,
                                   current.partitions,
                                   current.logit,
                                   current.z.contributions,
                                   current.z,
                                   presence.tables,
                                   nodes,
                                   effects,
                                   objects,
                                   neighborhood,
                                   numgroups.simulated,
                                   sizes.simulated)

    proba.change <- min(1, new.step$hastings.ratio)
    change.made <- (runif(1, 0, 1) <= proba.change)

    move <- new.step$move

    # update partition if needed
    old.partitions <- current.partitions
    if (change.made) {
      current.partitions <- new.step$new.partitions
      current.z.contributions <- new.step$new.z.contributions
      current.z <- new.step$new.z
      current.logit <- new.step$new.logit
    }

    cpt_burnin <- cpt_burnin + 1
    if (cpt_burnin > burnin) cpt_thining <- cpt_thining + 1

    if (cpt_burnin %% 10000 == 0 && verbose) cat(cpt_burnin, "\n")

    # store the results if we are out of burnin
    if (cpt_burnin >= burnin && cpt_thining == thining) {

      store <- TRUE
      if (constraints) { # check all partitions one by one
        for (o in 1:num.obs) store <- store && check_sizes(current.partitions[as.logical(presence.tables[, o]), o], sizes.allowed, numgroups.allowed)
      }

      if (store) {
        all.z <- rbind(all.z, current.z)
        cpt_thining <- 0
        cpt_steps <- cpt_steps + 1
        # store all partitions if needed
        if (return.all.partitions) {
          all.partitions[[cpt_steps]] <- current.partitions
        }
      } else {
        cpt_thining <- thining - 1
      }

    }

    # stop the walk if number of steps reached
    end.walk <- (cpt_steps >= num.steps)

  }

  # TODO: change this hack
  row.names(all.z) <- rep("", dim(all.z)[1])

  # compute the average statistics and the final network generated
  if (!return.all.partitions) {
    return(list("draws" = all.z,
                "last.partitions" = current.partitions,
                "all.partitions" = NULL))
  } else {
    return(list("draws" = all.z,
                "last.partitions" = current.partitions,
                "all.partitions" = all.partitions))
  }

}



# MULTIPLE PARTITIONS PROCEDURE



## --- FUNCTIONS TO DRAW ONE STEP IN THE CHAIN ----


# function to draw next partition and calculate HAstings ratio (one step in the Metropolis algorithm)
draw_step_single <- function(theta,
                             current.partition,
                             current.logit,
                             current.z,
                             nodes,
                             effects,
                             objects,
                             neighborhood,
                             numgroups.simulated,
                             sizes.simulated) {

  # pick a neighborhood
  # 1 = swap two actors,
  # 2 = merge 2 groups or split a group in two,
  # 3 = move one actor,
  ## 4 = move one pair of nodes (FOR NOW REMOVED)
  ## 5 = swap two pairs of nodes (FOR NOW REMOVED)
  ## 6 = reshuffle the members of two groups (FOR NOW REMOVED)
  n_neighborhood <- length(neighborhood)
  move <- sample(1:n_neighborhood, 1, prob = neighborhood)

  current.sizes <- rep(0, n_neighborhood)
  new.sizes <- rep(0, n_neighborhood)

  # calculate current size of neighborhood and pick new partition for the chosen move
  current.size <- compute_size_neighborhood(move, current.partition, numgroups.simulated, sizes.simulated)
  if (current.size$total > 0) {
    new.partition <- sample_new_partition(move, current.partition, current.size, numgroups.simulated, sizes.simulated)
    new.size <- compute_size_neighborhood(move, new.partition, numgroups.simulated, sizes.simulated)
  } else {
    new.partition <- current.partition
    new.size <- current.size
  }
  current.sizes[move] <- current.size$total
  new.sizes[move] <- new.size$total

  # if change, calculate sizes for all potential moves (depending on available neighborhoods)
  if (current.size$total > 0) {
    for (i in 1:n_neighborhood) {
      if (neighborhood[i] > 0 && i != move && reachable(i, current.partition, new.partition)) {
        cs <- compute_size_neighborhood(i, current.partition, numgroups.simulated, sizes.simulated)
        ns <- compute_size_neighborhood(i, new.partition, numgroups.simulated, sizes.simulated)
        current.sizes[i] <- cs$total
        new.sizes[i] <- ns$total
      }
    }
  }

  # compute new statistics only if it changed
  if (current.size$total > 0) {
    new.z <- computeStatistics(new.partition, nodes, effects, objects)
  } else {
    new.z <- current.z
  }

  new.logit <- theta * new.z

  # calculate acceptance ratio if needed
  if (current.size$total > 0) {
    indexes <- current.sizes > 0
    neighborhoods.ratio <- sum(neighborhood[indexes] / new.sizes[indexes]) /
      sum(neighborhood[indexes] / current.sizes[indexes])
  } else {
    neighborhoods.ratio <- 1
  }

  hastings.ratio <- (exp(sum(new.logit) - sum(current.logit))) * neighborhoods.ratio

  return(list("hastings.ratio" = hastings.ratio,
              "new.partition" = new.partition,
              "new.z" = new.z,
              "new.logit" = new.logit))
}




# function to draw next partition and calculate HAstings ratio (one step in the Metropolis algorithm)
draw_step_multiple <- function(theta,
                               current.partitions,
                               current.logit,
                               current.z.contributions,
                               current.z,
                               presence.tables,
                               nodes,
                               effects,
                               objects,
                               neighborhood,
                               numgroups.simulated,
                               sizes.simulated) {

  num.nodes <- nrow(nodes)
  num.obs <- ncol(presence.tables)
  new.partitions <- current.partitions

  rand.o <- sample(1:num.obs, 1)
  nodes.rand.o <- as.logical(presence.tables[, rand.o])
  current.partition <- current.partitions[nodes.rand.o, rand.o]

  # calculate sizes for all potential moves (depending on available neighborhoods)
  current.sizes <- rep(0, length(neighborhood))

  # pick a neighborhood
  # 1 = swap two actors,
  # 2 = merge 2 groups or split a group in two,
  # 3 = move one actor,
  # 4 = move one pair of nodes (FOR NOW REMOVED)
  # 5 = swap two pairs of nodes (FOR NOW REMOVED)
  # 6 = reshuffle the members of two groups (FOR NOW REMOVED)
  # 7 = reproduce a pair from the previous partition (empty for first time point) (FOR NOW REMOVED)
  n_neighborhood <- length(neighborhood)
  move <- sample(1:n_neighborhood, 1, prob = neighborhood)

  # calculate current size of neighborhood and pick new partition for neighborhoods 1 to 6
  for (i in 1:n_neighborhood) {
    if (neighborhood[i] > 0) {
      current.size <- compute_size_neighborhood(i, current.partition,  numgroups.simulated, sizes.simulated)
      current.sizes[i] <- current.size$total
      if (move == i) {
        if (current.size$total > 0) new.partition <- sample_new_partition(i, current.partition, current.size,  numgroups.simulated, sizes.simulated)
        if (current.size$total == 0) new.partition <- current.partition
        new.partitions[, rand.o] <- rep(NA, num.nodes)
        new.partitions[nodes.rand.o, rand.o] <- new.partition
      }
    }
  }
  # # this should be for neighborhood 7
  # if (neighborhood[7] > 0) {
  #   if (rand.o == 1) current.partition2 <- 1:num.nodes
  #   else current.partition2 <- current.partitions[, rand.o - 1]
  #   current.size <- compute_size_neighborhood(7, current.partitions[, rand.o], partition2 = current.partition2,  numgroups.simulated, sizes.simulated)
  #   current.sizes[7] <- current.size$total
  #   if (move == 7) {
  #     if (current.size$total > 0) new.partition <- sample_new_partition(7, current.partitions[, rand.o], current.size, partition2 = current.partition2,  numgroups.simulated, sizes.simulated)
  #     if (current.size$total == 0) new.partition <- current.partition
  #     new.partitions[, rand.o] <- new.partition
  #   }
  # }


  # compute new statistics only if it changed
  if (all(current.partitions == new.partitions, na.rm = TRUE)) {
    new.z.contributions <- current.z.contributions
    new.z <- current.z
  } else {
    recalculated.stats <- recalculate_statistics(new.partitions, rand.o, nodes.rand.o, nodes, effects, objects, current.z.contributions)
    new.z.contributions <- recalculated.stats$new.z.contributions
    new.z <- recalculated.stats$new.z
  }
  new.logit <- theta * new.z

  # chose whether to change or not
  new.sizes <- rep(0, n_neighborhood)
  for (i in 1:n_neighborhood) {
    if (neighborhood[i] > 0) {
      new.size <- compute_size_neighborhood(i, new.partition, numgroups.simulated, sizes.simulated)
      new.sizes[i] <- new.size$total
    }
  }
  # # this should be for neighborhood 7
  # if (neighborhood[7] > 0) {
  #   new.size <- compute_size_neighborhood(7, new.partition, current.partition2, numgroups.simulated, sizes.simulated)
  #   new.sizes[7] <- new.size$total
  # }
  neighborhoods.ratio <- sum(current.sizes * neighborhood) /
    sum(new.sizes * neighborhood)

  hastings.ratio <- (exp(sum(new.logit) - sum(current.logit))) * neighborhoods.ratio

  return(list("hastings.ratio" = hastings.ratio,
              "new.partitions" = new.partitions,
              "new.z.contributions" = new.z.contributions,
              "new.z" = new.z,
              "new.logit" = new.logit,
              "move" = move))
}



# Recalculate statistics for a changed partition in multiple estimation
recalculate_statistics <- function(new.partitions, rand.o, nodes.rand.o, nodes, effects, objects, current.z.contributions) {

  # store new statistics
  new.z.contributions <- current.z.contributions

  # calculate separately each effect
  for (e in 1:length(effects$names)) {

    new.z.contributions[e, rand.o] <- step_recalculate(new.partitions, rand.o, nodes.rand.o, nodes, effects, objects, e)
  }

  new.z <- rowSums(new.z.contributions)

  return(list(new.z.contributions = new.z.contributions,
              new.z = new.z))

}


step_recalculate <- function(new.partitions, rand.o, nodes.rand.o, nodes, effects, objects, e) {

  num.nodes <- nrow(nodes)
  num.obs <- ncol(new.partitions)
  object.name <- effects$objects[e]

  # TIE EFFECT: keep only present nodes
  if (effects$names[e] == "tie") {

    effects.temp <- list(names = "tie", objects = "net.temp")
    for (ob in 1:length(objects)) {
      if (objects[[ob]][[1]] == object.name) {
        net <- objects[[ob]][[2]]
        objects.temp <- list(list(name = "net.temp", object = net[nodes.rand.o, nodes.rand.o]))
      }
    }
    nodes.temp <- nodes[nodes.rand.o, ]
    new.z.contribution <- as.numeric(computeStatistics(new.partitions[nodes.rand.o, rand.o], nodes.temp, effects.temp, objects.temp))

    #  SAME_VAR EFFECT: adapt the varying covariate
  } else if (effects$names[e] == "same_var") {

    effects.temp <- list(names = "same", objects = "var.temp")
    for (ob in 1:length(objects)) {
      if (objects[[ob]][[1]] == object.name) {
        atts <- objects[[ob]][[2]]
      }
    }
    objects.temp <- list()
    nodes.temp <- nodes[nodes.rand.o, ]
    nodes.temp$var.temp <- atts[nodes.rand.o, rand.o]
    new.z.contribution <- as.numeric(computeStatistics(new.partitions[nodes.rand.o, rand.o], nodes.temp, effects.temp, objects.temp))

    #  TIE_VAR EFFECT: adapt the varying network
  } else if (effects$names[e] == "tie_var") {

    effects.temp <- list(names = "tie", objects = "net.temp")
    for (ob in 1:length(objects)) {
      if (objects[[ob]][[1]] == object.name) {
        nets <- objects[[ob]][[2]]
        net <- nets[[rand.o]]
        objects.temp <- list(list(name = "net.temp", object = net[nodes.rand.o, nodes.rand.o]))
      }
    }
    nodes.temp <- nodes[nodes.rand.o, ]
    new.z.contribution <- as.numeric(computeStatistics(new.partitions[nodes.rand.o, rand.o], nodes.temp, effects.temp, objects.temp))

    #  TIE_VAR_X_DIFF EFFECT: adapt the varying network
  } else if (effects$names[e] == "tie_var_X_diff") {

    effects.temp <- list(names = "tie_X_diff", objects = "net.temp", objects2 = effects$objects2[e])
    for (ob in 1:length(objects)) {
      if (objects[[ob]][[1]] == object.name) {
        nets <- objects[[ob]][[2]]
        net <- nets[[rand.o]]
        objects.temp <- list(list(name = "net.temp", object = net[nodes.rand.o, nodes.rand.o]))
      }
    }
    nodes.temp <- nodes[nodes.rand.o, ]
    new.z.contribution <- as.numeric(computeStatistics(new.partitions[nodes.rand.o, rand.o], nodes.temp, effects.temp, objects.temp))


    # INERTIA_1 EFFECT: need to recalculate the partition before and after
  } else if (effects$names[e] == "inertia_1") {

    aff.temp <- as.matrix(table(data.frame(actor = 1:num.nodes, group = new.partitions[, rand.o])))
    adj.temp <- aff.temp %*% t(aff.temp)
    diag(adj.temp) <- 0

    if (rand.o > 1) { # removed  && length(groups[[rand.o]]) > 0
      aff.temp2 <- as.matrix(table(data.frame(actor = 1:num.nodes, group = new.partitions[, rand.o - 1])))
      adj.temp2 <- aff.temp2 %*% t(aff.temp2)
      diag(adj.temp2) <- 0
      adj12 <- adj.temp * adj.temp2
      new.z.contribution <- sum(0.5 * adj12)
    } else {
      new.z.contribution <- 0
    }

    if (rand.o < num.obs) { # removed  && length(groups[[rand.o]]) > 0
      aff.temp2 <- as.matrix(table(data.frame(actor = 1:num.nodes, group = new.partitions[, rand.o + 1])))
      adj.temp2 <- aff.temp2 %*% t(aff.temp2)
      diag(adj.temp2) <- 0
      adj12 <- adj.temp * adj.temp2
      new.z.contribution <- sum(0.5 * adj12)
    }

    # ANY OTHER EFFECT: the effect also exists in the single partition version
  } else {

    effects.temp <- list(names = effects$names[e], objects = effects$objects[e])
    objects.temp <- objects
    nodes.temp <- nodes[nodes.rand.o, ]
    new.z.contribution <- as.numeric(computeStatistics(new.partitions[nodes.rand.o, rand.o], nodes.temp, effects.temp, objects.temp))

  }

  return(new.z.contribution)

}



## --- FUNCTIONS FOR THE PROPOSALS MADE AT EACH STEP ----

## GENERIC FUNCTION TO COMPUTE NEIGHBORHOOD SIZE
compute_size_neighborhood <- function(i,
                                      partition,
                                      numgroups.simulated = NULL,
                                      sizes.simulated = NULL) {#, partition2 = NULL, # for neighborhood 7, now deprecated) {

  # check whether there are constraints
  constraints <- FALSE
  if (!is.null(numgroups.simulated) || !is.null(sizes.simulated)) {
    constraints <- TRUE
  }

  if (i == 1) {
    if (!constraints) return(compute_size_neighborhood_p1(partition))
    else return(compute_size_neighborhood_p1_restricted(partition, numgroups.simulated, sizes.simulated))
  }
  if (i == 2) {
    if (!constraints) return(compute_size_neighborhood_p2(partition))
    else return(compute_size_neighborhood_p2_restricted(partition, numgroups.simulated, sizes.simulated))
  }
  if (i == 3) {
    if (!constraints) return(compute_size_neighborhood_p3(partition))
    else return(compute_size_neighborhood_p3_restricted(partition, numgroups.simulated, sizes.simulated))
  }

  # For now, other neighborhoods removed
  # if (i == 4) {
  #   if (is.null(sizes.simulated)) return(compute_size_neighborhood_p4(partition))
  #   else return(compute_size_neighborhood_p4_restricted(partition, numgroups.simulated, sizes.simulated))
  # }
  # if (i == 5) {
  #   if (is.null(sizes.simulated)) return(compute_size_neighborhood_p5(partition))
  #   else return(compute_size_neighborhood_p5_restricted(partition, numgroups.simulated, sizes.simulated))
  # }
  # if (i == 6) {
  #   if (is.null(sizes.simulated)) return(compute_size_neighborhood_p6(partition))
  #   else return(compute_size_neighborhood_p6_restricted(partition, numgroups.simulated, sizes.simulated))
  # }
  # if (i == 7) {
  #   if (is.null(sizes.simulated)) return(compute_size_neighborhood_p7(partition))
  #   else return(compute_size_neighborhood_p7_restricted(partition, numgroups.simulated, sizes.simulated))
  # }
}

## GENERIC FUNCTION TO COMPUTE NEIGHBORHOOD SIZE
sample_new_partition <- function(i,
                                 current.partition,
                                 size_neighborhood,
                                 numgroups.simulated = NULL,
                                 sizes.simulated = NULL) {#, current.partition2 = NULL, # for neighborhood 7, now deprecated) {

  # check whether there are constraints
  constraints <- FALSE
  if (!is.null(numgroups.simulated) || !is.null(sizes.simulated)) {
    constraints <- TRUE
  }

  if (i == 1) {
    if (!constraints) return(sample_new_partition_p1(current.partition, size_neighborhood))
    else return(sample_new_partition_p1_restricted(current.partition, size_neighborhood, numgroups.simulated, sizes.simulated))
  }
  if (i == 2) {
    if (!constraints) return(sample_new_partition_p2(current.partition, size_neighborhood))
    else return(sample_new_partition_p2_restricted(current.partition, size_neighborhood, numgroups.simulated, sizes.simulated))
  }
  if (i == 3) {
    if (!constraints) return(sample_new_partition_p3(current.partition, size_neighborhood))
    else return(sample_new_partition_p3_restricted(current.partition, size_neighborhood, numgroups.simulated, sizes.simulated))
  }

  # for now, other neighborhoods removed
  # if (i == 4) {
  #   if (is.null(sizes.simulated)) return(sample_new_partition_p4(current.partition, size_neighborhood))
  #   else return(sample_new_partition_p4_restricted(current.partition, size_neighborhood, numgroups.simulated, sizes.simulated))
  # }
  # if (i == 5) {
  #   if (is.null(sizes.simulated)) return(sample_new_partition_p5(current.partition, size_neighborhood))
  #   else return(sample_new_partition_p5_restricted(current.partition, size_neighborhood, numgroups.simulated, sizes.simulated))
  # }
  # if (i == 6) {
  #   if (is.null(sizes.simulated)) return(sample_new_partition_p6(current.partition, size_neighborhood))
  #   else return(sample_new_partition_p6_restricted(current.partition, size_neighborhood, numgroups.simulated, sizes.simulated))
  # }
  # if (i == 7) {
  #   if (is.null(sizes.simulated)) return(sample_new_partition_p7(current.partition, size_neighborhood, current.partition2))
  #   else return(sample_new_partition_p7_restricted(current.partition, size_neighborhood, current.partition2, numgroups.simulated, sizes.simulated))
  # }
}

## GENERIC FUNCTION TO TEST WHETHER A PARTITION IS REACHABLE WITH A GIVEN NEIGHBORHOOD
reachable <- function(i, partition1, partition2) {
  if (i == 1) {
    return(reachable_p1(partition1, partition2))
  } else if (i == 2) {
    return(reachable_p2(partition1, partition2))
  } else if (i == 3) {
    return(reachable_p3(partition1, partition2))
  }

  # for now, other nieghborhoods, removed
  #  else if (i == 4) {
  #   return(reachable_p4(partition1, partition2))
  # } else if (i == 5) {
  #   return(reachable_p5(partition1, partition2))
  # } else if (i == 6) {
  #   return(reachable_p6(partition1, partition2))
  # }
}


## --- SPECIFIC PROPOSALS ----


## NEIGHBORHOOD PI 1: only swaps of two nodes (careful, cannot be used alone)
compute_size_neighborhood_p1 <- function(partition) {

  # find isolates, pairs and groups > 2
  num.groups <- max(partition)
  num.nodes <- length(partition)
  sizes <- table(partition)
  isolates <- as.vector(which(sizes == 1))
  pairs <- as.vector(which(sizes == 2))
  others <- as.vector(which(sizes > 2))

  nums.swaps <- matrix(0, num.nodes, num.nodes)
  for (i in 1:(num.nodes - 1)) {
    for (j in (i + 1):num.nodes) {
      g_i <- partition[i]
      g_j <- partition[j]

      allowed <- TRUE
      if (g_i == g_j) allowed <- FALSE # two members of the same group cannot swap, it wouldn't change anything
      if (g_i %in% isolates && g_j %in% isolates) {
        allowed <- FALSE # two isolates cannot swap, it wouldn't change anything
      }

      if (allowed) nums.swaps[i, j] <- 1
    }
  }

  return(list(isolates = isolates,
              pairs = pairs,
              others = others,
              nums.swaps = nums.swaps,
              num.swaps = sum(nums.swaps),
              total = sum(nums.swaps)))

}

# sample partition neighborhood 1
sample_new_partition_p1 <- function(current.partition, size_neighborhood) {

  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- size_neighborhood
  isolates <- size$isolates
  pairs <- size$pairs
  others <- size$others
  nums.swaps <- size$nums.swaps
  num.swaps <- size$num.swaps
  total <- size$total

  if (total == 0) return(current.partition)

  pick.1 <- sample(total, 1)

  # decide which actors to swap
  new.partition <- current.partition
  swap <- which(nums.swaps == 1)[pick.1]
  j <- ((swap - 1) %/% num.nodes) + 1
  i <- (swap - 1) %% num.nodes + 1
  new.partition[i] <- current.partition[j]
  new.partition[j] <- current.partition[i]
  new.partition <- order_groupids(new.partition)

  return(new.partition)

}

# reachable neighborhood 1
reachable_p1 <- function(partition1, partition2) {

  # they are not neighbors if they have a different number of groups
  if (max(partition1) != max(partition2)) return(FALSE)

  num.nodes <- length(partition1)
  check <- FALSE
  i <- 1
  j <- i + 1

  # try remove all pairs of nodes
  while (!check && i <= num.nodes) {

    newp1 <- order_groupids(partition1[-c(i, j)])
    newp2 <- order_groupids(partition2[-c(i, j)])

    # check if they swapped
    if (all(newp1 == newp2)) {
      othersi1 <- which(partition1 == partition1[i])
      othersi1 <- othersi1[othersi1 != i]
      othersj1 <- which(partition1 == partition1[j])
      othersj1 <- othersj1[othersj1 != j]
      othersi2 <- which(partition2 == partition2[i])
      othersi2 <- othersi2[othersi2 != i]
      othersj2 <- which(partition2 == partition2[j])
      othersj2 <- othersj2[othersj2 != j]
      check <- (setequal(othersi1, othersj2) && setequal(othersi2, othersj1))
      if (check) break
    }

    j <- j + 1
    if (j > num.nodes) {
      i <- i + 1
      j <- i + 1
    }
  }

  return(check)
}



## NEIGHBORHOOD PI 2: only merges and divisions of 2 groups
compute_size_neighborhood_p2 <- function(partition) {

  # calculate the number of neighbor partitions
  num.nodes <- length(partition)
  num.groups <- max(partition)
  num.merges <- num.groups * (num.groups - 1) * 0.5
  nums.divisions <- rep(0, num.groups)
  num.divisions <- 0

  for (k in 1:num.groups) {
    sk <- length(which(partition == k))
    nums.divisions[k] <- 2^(sk - 1) - 1
    num.divisions <- num.divisions + nums.divisions[k]
  }

  return(list(num.merges = num.merges,
              num.divisions = num.divisions,
              nums.divisions = nums.divisions,
              total = num.merges + num.divisions))
}

# sample partition neighborhood 2
sample_new_partition_p2 <- function(current.partition, size_neighborhood) {

  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- size_neighborhood
  num.merges <- size$num.merges
  num.divisions <- size$num.divisions
  nums.divisions <- size$nums.divisions
  total <- size$total

  if (total == 0) return(current.partition)

  # decide between merge or division (or self loop)
  pick.1 <- sample(total, 1)
  new.partition <- current.partition

  # if merge
  if (pick.1 <= num.merges) {

    # pick 2 groups to merge
    old_gs <- sample(num.groups, 2, replace = FALSE)
    # reassign one of the groups and remove useless ids
    new.partition[which(new.partition == old_gs[2])] <- old_gs[1]
    new.partition <- order_groupids(new.partition)

  }else if (pick.1 <= (num.merges + num.divisions)) {

    # pick group to divide
    pick.2 <- sample(num.divisions, 1)
    sums <- unlist(lapply(1:num.groups, FUN = function(x) {
      sum(nums.divisions[1:x])
    }))
    starts <- c(0, sums[1:(num.groups - 1)])
    ends <- sums[1:num.groups]
    old_g <- which(starts < pick.2 & ends >= pick.2)

    # reassign one of the groups and remove useless ids
    new.groups <- sample(c(old_g, num.groups + 1), length(which(new.partition == old_g)), replace = TRUE)
    found <- (length(unique(new.groups)) > 1)
    while (!found) {
      new.groups <- sample(c(old_g, num.groups + 1), length(which(new.partition == old_g)), replace = TRUE)
      found <- (length(unique(new.groups)) > 1)
    }

    new.partition[which(new.partition == old_g)] <- new.groups
    new.partition <- order_groupids(new.partition)
  }

  return(new.partition)

}

# reachable neighborhood 2
reachable_p2 <- function(partition1, partition2) {

  # they are not neighbors if one does not have one more group
  if (abs(max(partition1) - max(partition2)) != 1) return(FALSE)

  num.nodes <- length(partition1)
  num.groups1 <- max(partition1)
  num.groups2 <- max(partition2)
  check <- FALSE
  g1 <- 1
  g2 <- g1 + 1

  # try remove groups (to check divisions) and pairs of groups (to check merges)
  while (!check && g1 <= num.groups1) {

    # first check whether g1 is not divided in partition2
    membersg1 <- which(partition1 == g1)
    if (length(membersg1) == num.nodes) {
      newp1 <- c()
      newp2 <- c()
    } else {
      newp1 <- order_groupids(partition1[-membersg1])
      newp2 <- order_groupids(partition2[-membersg1])
    }

    if (all(newp1 == newp2)) {
      check <- length(unique(partition2[membersg1])) == 2
      if (check) break
    }

    # second check whether g1 and g2 are merged in partition2
    membersg2 <- which(partition1 == g2)
    if (length(c(membersg1, membersg2)) == num.nodes) {
      newp1 <- c()
      newp2 <- c()
    } else {
      newp1 <- order_groupids(partition1[-c(membersg1, membersg2)])
      newp2 <- order_groupids(partition2[-c(membersg1, membersg2)])
    }

    if (all(newp1 == newp2)) {
      check <- length(unique(partition2[c(membersg1, membersg2)])) == 1
      if (check) break
    }

    g2 <- g2 + 1
    if (g2 > num.nodes) {
      g1 <- g1 + 1
      g2 <- g1 + 1
    }
  }

  return(check)
}


## NEIGHBORHOOD PI 3: only swaps of one node
compute_size_neighborhood_p3 <- function(partition) {

  # find isolates, pairs and groups > 2
  num.groups <- max(partition)
  num.nodes <- length(partition)
  sizes <- table(partition)
  isolates <- as.vector(which(sizes == 1))
  pairs <- as.vector(which(sizes == 2))
  others <- as.vector(which(sizes > 2))

  nums.swaps <- rep(0, num.nodes)
  done.pairs <- rep(0, length(pairs))
  done.isolates <- rep(0, length(isolates))
  for (i in 1:num.nodes) {
    g <- partition[i]
    if (g %in% isolates) {
      done.isolates[isolates == g] <- 1
      nums.swaps[i] <- num.groups - sum(done.isolates) # can join any other group, except other isolates before (because already done) and itself
    }
    if (g %in% pairs) {
      nums.swaps[i] <- num.groups - done.pairs[pairs == g] # can join any other group or create its isolate, except if the possibilty of splitting the pair in 2 isolates is already counted
      done.pairs[pairs == g] <- 1
    }
    if (g %in% others) {
      nums.swaps[i] <- num.groups # can join any other group or create its own isolate
    }
  }

  return(list(isolates = isolates,
              pairs = pairs,
              others = others,
              nums.swaps = nums.swaps,
              num.swaps = sum(nums.swaps),
              total = sum(nums.swaps)))

}


# sample partition neighborhood 3
sample_new_partition_p3 <- function(current.partition, size_neighborhood) {

  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- size_neighborhood
  isolates <- size$isolates
  pairs <- size$pairs
  others <- size$others
  nums.swaps <- size$nums.swaps
  num.swaps <- size$num.swaps
  total <- size$total

  if (total == 0) return(current.partition)

  pick.1 <- sample(total, 1)

  # decide which actor to swap
  all.groups <- 1:num.groups
  new.partition <- current.partition
  done.pairs <- rep(0, length(pairs))
  done.isolates <- rep(0, length(isolates))

  for (i in 1:num.nodes) {

    if (i == 1) {
      start <- 0
    } else {
      start <- sum(nums.swaps[1:i - 1])
    }
    if (i == num.nodes) {
      end <- num.swaps
    } else {
      end <- sum(nums.swaps[1:i])
    }

    g <- current.partition[i]
    if (g %in% isolates) done.isolates[isolates == g] <- 1

    if (pick.1 > start && pick.1 <= end) {
      g <- current.partition[i]

      # an isolate is randomly joining another group except isolates before i
      if (g %in% isolates) {
        previousisolates <- isolates[done.isolates == 1]
        tosample <- all.groups[!all.groups %in% previousisolates]
        if (length(tosample) == 1) newg <- tosample
        if (length(tosample) >= 2) newg <- sample(tosample, 1)
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }

      # a member of a pair is randomly joining another group
      # or creating its isolate if it's the first of the pair to appear
      if (g %in% pairs) {
        if (done.pairs[pairs == g] == 0) {
          newg <- sample(all.groups, 1)
          if (newg == g) {
            newg <- num.groups + 1
          }
        } else {
          tosample <- all.groups[all.groups != g]
          if (length(tosample) == 1) newg <- tosample
          if (length(tosample) >= 2) newg <- sample(tosample, 1)
        }
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }

      # a member of a bigger group is randomly joining another group
      # or creating its isolate
      if (g %in% others) {
        newg <- sample(all.groups, 1)
        if (newg == g) {
          newg <- num.groups + 1
        }
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }
    }

    if (g %in% pairs) done.pairs[pairs == g] <- 1
  }

  return(new.partition)

}

# reachable neighborhood 3
reachable_p3 <- function(partition1, partition2) {

  # they are not neighbors if one doesn't have the same number of groups or one more
  if (abs(max(partition1) - max(partition2)) > 1) return(FALSE)

  num.nodes <- length(partition1)
  check <- FALSE
  i <- 1

  # try remove all nodes
  while (!check && i <= num.nodes) {

    newp1 <- order_groupids(partition1[-i])
    newp2 <- order_groupids(partition2[-i])

    # check if it moved
    if (all(newp1 == newp2)) {
      othersi1 <- which(partition1[-i] == partition1[i])
      othersi2 <- which(partition2[-i] == partition2[i])
      check <- !setequal(othersi1, othersi2)
      if (check) break
    }

    i <- i + 1
  }

  return(check)
}


## --- SPECIFIC PROPOSALS RESTRICTED VERSION ----

## NEIGHBORHOOD PI 1 RESTRICTED: only swaps of two nodes (careful, cannot be used alone)
# WARNING!!!!: for now it only works if sizes allowed are an interval ([size_min, size_max])
compute_size_neighborhood_p1_restricted <- function(partition, numgroups.simulated, sizes.simulated) {
  # check if current partition is allowed
  if (!check_sizes(partition, sizes.simulated, numgroups.simulated)) stop("The partition we are in is not allowed.")
  return(compute_size_neighborhood_p1(partition))
}

# sample partition neighborhood 1 restricted
sample_new_partition_p1_restricted <- function(current.partition, size_neighborhood, numgroups.simulated, sizes.simulated) {
  return(sample_new_partition_p1(current.partition, size_neighborhood))
}


## NEIGHBORHOOD PI 2 RESTRICTED: only merges and divisions of 2 groups
compute_size_neighborhood_p2_restricted <- function(partition, numgroups.simulated, sizes.simulated) {

  # check if current partition is allowed
  if (!check_sizes(partition, sizes.simulated, numgroups.simulated)) stop("The partition we are in is not allowed.")

  # find minimum and maximum size allowed
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)

  # calculate the number of neighbor partitions
  num.nodes <- length(partition)
  num.groups <- max(partition)
  sizes <- table(partition)
  num.merges <- 0
  merges <- matrix(0, num.groups, num.groups)
  nums.divisions <- rep(0, num.groups)
  num.divisions <- 0

  # merges
  if (num.groups > 1 && (num.groups - 1) %in% numgroups.simulated) { # we count merges if we have more than one group and having one less group is allowed
    for (g1 in 1:(num.groups - 1)) {
      for (g2 in (g1 + 1):num.groups) {
        if (length(which(partition == g1)) + length(which(partition == g2)) <= smax) { # merges are only possible if the new group size is allowed
          num.merges <- num.merges + 1
          merges[g1, g2] <- 1
        }
      }
    }
  }


  # divisions
  if ((num.groups + 1) %in% numgroups.simulated) { # we count divisions if having one more group is allowed
    for (k in 1:num.groups) {
      sg <- sizes[k]
      if (smin > 1) {
        extras <- sum(unlist(lapply(1:(smin - 1), FUN = function(x) {
          choose(sg, x)
        })))
        if ((sg %% 2) == 0 && sg * 0.5 < smin) extras <- extras - choose(sg, sg * 0.5) * 0.5
      } else {
        extras <- 0
      }
      nums.divisions[k] <- 2^(sg - 1) - 1 - extras
      num.divisions <- num.divisions + nums.divisions[k]
    }
  }


  return(list(num.merges = num.merges,
              merges = merges,
              num.divisions = num.divisions,
              nums.divisions = nums.divisions,
              total = num.merges + num.divisions))
}



# sample partition neighborhood 2 restricted
sample_new_partition_p2_restricted <- function(current.partition, size_neighborhood, numgroups.simulated, sizes.simulated) {

  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- size_neighborhood

  num.merges <- size$num.merges
  merges <- size$merges
  num.divisions <- size$num.divisions
  nums.divisions <- size$nums.divisions
  total <- size$total

  if (total == 0) return(current.partition)

  # sizes allowed
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)

  # decide between merge or division (or self loop)
  pick.1 <- sample(total, 1)
  all.groups <- 1:num.groups
  new.partition <- current.partition

  # if merge
  if (pick.1 <= num.merges) {

    # pick 2 groups to merge
    allindexes <- which(merges == 1)
    p <- allindexes[pick.1]
    indexes <- arrayInd(p, dim(merges))
    old_g1 <- indexes[, 1]
    old_g2 <- indexes[, 2]

    # reassign one of the groups and remove useless ids
    new.partition[which(new.partition == old_g2)] <- old_g1
    new.partition <- order_groupids(new.partition)

  } else if (pick.1 <= (num.merges + num.divisions)) {

    # pick group to divide
    pick.2 <- sample(num.divisions, 1)
    old_g1 <- 1
    cpt2 <- 1
    g <- 1
    found <- FALSE
    while (g <= num.groups && !found) {
      if (pick.2 >= cpt2 && pick.2 <= sum(nums.divisions[1:g])) {
        old_g1 <- g
        found <- TRUE
      }
      cpt2 <- cpt2 + nums.divisions[g]
      g <- g + 1
    }
    old_g2 <- 0
    new_g1 <- old_g1
    new_g2 <- num.groups + 1

    # reassign one of the groups and remove useless ids
    new.groups <- sample(c(old_g1, num.groups + 1), length(which(new.partition == old_g1)), replace = TRUE)
    found <- (length(unique(new.groups)) > 1) && min(table(new.groups)) >= smin
    while (!found) {
      new.groups <- sample(c(old_g1, num.groups + 1), length(which(new.partition == old_g1)), replace = TRUE)
      found <- (length(unique(new.groups)) > 1) && min(table(new.groups)) >= smin
    }
    new.partition[which(new.partition == old_g1)] <- new.groups
    new.partition <- order_groupids(new.partition)
  }

  return(new.partition)

}




## NEIGHBORHOOD PI 3 RESTRICTED: only swaps one node
# WARNING!!!!: for now it only works if sizes allowed are an interval ([size_min, size_max])
compute_size_neighborhood_p3_restricted <- function(partition, numgroups.simulated, sizes.simulated) {

  # check if current partition is allowed
  if (!check_sizes(partition, sizes.simulated, numgroups.simulated)) stop("The partition we are in is not allowed.")

  # find maximum size allowed
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)
  sizes <- table(partition)

  # find isolates, pairs and groups > 2
  num.groups <- max(partition)
  num.nodes <- length(partition)
  isolates <- as.vector(which(sizes == 1))
  pairs <- as.vector(which(sizes == 2))
  others <- as.vector(which(sizes > 2))
  abovemin_groups <- as.vector(which(sizes > smin))
  belowmax_groups <- as.vector(which(sizes < smax))

  nums.swaps <- rep(0, num.nodes)
  done.pairs <- rep(0, length(pairs))
  done.isolates <- rep(0, length(isolates))
  for (i in 1:num.nodes) {
    g <- partition[i]

    if (g %in% isolates && (num.groups - 1) %in% numgroups.simulated) {
      # i can join any other group (<smax), except other isolates done before (because already done) and this one (that would not change anything)
      # since it's an isolate, it means that isolates are allowed in that context so we don't worry about smin (=1)
      done.isolates[isolates == g] <- 1
      nums.swaps[i] <- length(belowmax_groups) - sum(done.isolates)
    } # to move an isolate, a partition with one group less must be allowed

    if (g %in% pairs && g %in% abovemin_groups) {
      # if isolates are allowed, then i can join any other group (smax), and splitting if it's not already counted
      # here we don't worry about smin (=1) because 2 > smin
      gnotbelowmax <- !(g %in% belowmax_groups)
      if ((num.groups + 1) %in% numgroups.simulated && done.pairs[pairs == g] == 0) { # we check whether it's possible to create an isolate
        nums.swaps[i] <- length(belowmax_groups) + gnotbelowmax
      } else {
        nums.swaps[i] <- length(belowmax_groups) + gnotbelowmax - 1
      }
      done.pairs[pairs == g] <- 1
    } # to break a pair, size 2 should not be the minimum allowed

    if (g %in% others && g %in% abovemin_groups) {
      # can join any other group (<smax) or create its own isolate (if it's allowed)
      # we don't care here about checking that other groups are above min size, if they exist it means they are already allowed
      gnotbelowmax <- !(g %in% belowmax_groups)
      if (smin == 1 && (num.groups + 1) %in% numgroups.simulated) { # we check whether it's possible to create an isolate
        nums.swaps[i] <- length(belowmax_groups) + gnotbelowmax
      } else {
        nums.swaps[i] <- length(belowmax_groups) + gnotbelowmax - 1
      }
    } # to break a group, it should not be of minimal size
  }

  return(list(isolates = isolates,
              pairs = pairs,
              others = others,
              abovemin_groups = abovemin_groups,
              belowmax_groups = belowmax_groups,
              nums.swaps = nums.swaps,
              num.swaps = sum(nums.swaps),
              total = sum(nums.swaps)))

}

# sample partition neighborhood 3 restricted
sample_new_partition_p3_restricted <- function(current.partition, size_neighborhood, numgroups.simulated, sizes.simulated) {

  # calculate the number of neighbor partitions
  num.nodes <- length(current.partition)
  num.groups <- max(current.partition)
  size <- size_neighborhood
  isolates <- size$isolates
  pairs <- size$pairs
  others <- size$others
  abovemin_groups <- size$abovemin_groups
  belowmax_groups <- size$belowmax_groups
  nums.swaps <- size$nums.swaps
  num.swaps <- size$num.swaps
  total <- size$total

  if (total == 0) return(current.partition)

  # allowed sizes
  smax <- max(sizes.simulated)
  smin <- min(sizes.simulated)

  pick.1 <- sample(total, 1)

  # decide which actor to swap
  all.groups <- 1:num.groups
  new.partition <- current.partition
  done.pairs <- rep(0, length(pairs))
  done.isolates <- rep(0, length(isolates))

  for (i in 1:num.nodes) {

    if (i == 1) {
      start <- 0
    } else {
      start <- sum(nums.swaps[1:i - 1])
    }
    if (i == num.nodes) {
      end <- num.swaps
    } else {
      end <- sum(nums.swaps[1:i])
    }

    g <- current.partition[i]
    if (g %in% isolates) done.isolates[isolates == g] <- 1

    if (pick.1 > start && pick.1 <= end) {

      # an isolate can join any other group (<smax), except other isolates done before (because already done) and this one (that would not change anything)
      # if we are here, this implies that a partition with one less group is allowed
      if (g %in% isolates) {
        previousisolates <- isolates[done.isolates == 1]
        tosample <- all.groups[(all.groups %in% belowmax_groups) & !(all.groups %in% previousisolates)]
        if (length(tosample) == 1) newg <- tosample
        if (length(tosample) >= 2) newg <- sample(tosample, 1)
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }

      # a member of a pair can join any other group (<smax), and splitting if it's not already counted, if isolates are allowed
      if (g %in% pairs && g %in% abovemin_groups) {
        if ((num.groups + 1) %in% numgroups.simulated && done.pairs[pairs == g] == 0) { # if a split of the pair is possible
          tosample <- all.groups[(all.groups %in% belowmax_groups) | (all.groups == g)]
          if (length(tosample) == 1) newg <- tosample
          if (length(tosample) >= 2) newg <- sample(tosample, 1)
          if (newg == g) {
            newg <- num.groups + 1
          }
        } else {
          tosample <- all.groups[(all.groups %in% belowmax_groups) & (all.groups != g)]
          if (length(tosample) == 1) newg <- tosample
          if (length(tosample) >= 2) newg <- sample(tosample, 1)
        }
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }

      # a member of a bigger group (>smin) can join any other group (<smax) or create its own isolate (if it's allowed)
      if (g %in% others && g %in% abovemin_groups) {
        if (smin == 1 && (num.groups + 1) %in% numgroups.simulated) { # if we can create a isolate
          tosample <- all.groups[(all.groups %in% belowmax_groups) | (all.groups == g)]
          if (length(tosample) == 1) newg <- tosample
          if (length(tosample) >= 2) newg <- sample(tosample, 1)
          if (newg == g) {
            newg <- num.groups + 1
          }
        } else {
          tosample <- all.groups[(all.groups %in% belowmax_groups) & (all.groups != g)]
          if (length(tosample) == 1) newg <- tosample
          if (length(tosample) >= 2) newg <- sample(tosample, 1)
        }
        new.partition[i] <- newg
        new.partition <- order_groupids(new.partition)
      }

    }

    if (g %in% pairs) done.pairs[pairs == g] <- 1
  }

  return(new.partition)

}


### Deprecated: unused neighborhoods
### Can be reintroduced, but they don't work with constraints on number of groups (only on sizes)

# ## NEIGHBORHOOD PI 4: swaps a pair of nodes
#
# compute_size_neighborhood_p4 <- function(partition) {
#
#   # find isolates, pairs and groups > 2
#   num.groups <- max(partition)
#   num.nodes <- length(partition)
#   sizes <- table(partition)
#   isolates <- as.vector(which(sizes == 1))
#   pairs <- as.vector(which(sizes == 2))
#   others <- as.vector(which(sizes > 2))
#
#   nums.swaps <- rep(0, max(partition))
#
#   for (g in 1:num.groups) {
#
#     #if (g %in% isolates) {}# nothing happens, there is no pair to swap here
#
#     if (g %in% pairs) {
#       nums.swaps[g] <- length(pairs[pairs > g]) + length(isolates) + length(others)
#     }
#     if (g %in% others) {
#       sizeg <- sum(partition == g)
#       nums.swaps[g] <- (length(pairs) + length(others) + length(isolates)) * sizeg * (sizeg - 1) / 2
#     }
#   }
#
#   return(list(isolates = isolates,
#               pairs = pairs,
#               others = others,
#               nums.swaps = nums.swaps,
#               num.swaps = sum(nums.swaps),
#               total = sum(nums.swaps)))
#
# }
#
#
# sample_new_partition_p4 <- function(current.partition, size_neighborhood) {
#
#   # calculate the number of neighbor partitions
#   num.nodes <- length(current.partition)
#   num.groups <- max(current.partition)
#   size <- size_neighborhood
#   isolates <- size$isolates
#   pairs <- size$pairs
#   others <- size$others
#   nums.swaps <- size$nums.swaps
#   num.swaps <- size$num.swaps
#   total <- size$total
#
#   if (total == 0) return(current.partition)
#
#   pick.1 <- sample(total, 1)
#
#   # decide which group to pick a pair from
#   all.groups <- 1:num.groups
#   new.partition <- current.partition
#   #done.pairs <- rep(0, length(pairs))
#
#   sums <- unlist(lapply(1:num.groups, FUN = function(x) {sum(nums.swaps[1:x])}))
#   starts <- c(0, sums[1:(num.groups - 1)])
#   ends <- sums[1:num.groups]
#   old_g <- which(starts < pick.1 & ends >= pick.1)
#
#
#   # a pair can join other available groups, unless it's a pair that has already been "joined"
#   if (old_g %in% pairs) {
#     members <- which(current.partition == old_g)
#     tosample <- all.groups[!(all.groups %in% pairs & all.groups <= old_g)]
#     if (length(tosample) == 1) newg <- tosample else newg <- sample(tosample, 1)
#     new.partition[members] <- newg
#     new.partition <- order_groupids(new.partition)
#   }
#
#   # a pair in a bigger group can join any other group or create its own isolated pair
#   if (old_g %in% others) {
#     members <- which(current.partition == old_g)
#     sizeg <- length(members)
#     pairg <- sample(members, 2)
#     tosample <- all.groups
#     if (length(tosample) == 1) newg <- tosample else newg <- sample(tosample, 1)
#     if (newg == old_g) newg <- num.groups + 1
#     new.partition[pairg] <- newg
#     new.partition <- order_groupids(new.partition)
#   }
#
#   return(new.partition)
#
# }
#
#
# reachable_p4 <- function(partition1, partition2) {
#
#   # they are not neighbors if one doesn't have the same number of groups or one more
#   if (abs(max(partition1) - max(partition2)) > 1) return(FALSE)
#
#   num.nodes <- length(partition1)
#   check <- FALSE
#   i <- 1
#   j <- i + 1
#
#   # try remove all pairs of nodes
#   while (!check && i <= num.nodes) {
#
#     if (partition1[i] == partition1[j] && partition2[i] == partition2[j]) {
#       newp1 <- order_groupids(partition1[-c(i, j)])
#       newp2 <- order_groupids(partition2[-c(i, j)])
#
#       # check if the pair moved together
#       if (all(newp1 == newp2)) {
#         others1 <- which(partition1[-c(i, j)] == partition1[i])
#         others2 <- which(partition2[-c(i, j)] == partition2[i])
#         check <- !setequal(others1, others2)
#         if (check) break
#       }
#     }
#
#     j <- j + 1
#     if (j > num.nodes) {
#       i <- i + 1
#       j <- i + 1
#     }
#   }
#
#   return(check)
# }
#
#
# ## NEIGHBORHOOD PI 5: swap two pairs of nodes
#
# compute_size_neighborhood_p5 <- function(partition) {
#
#   # find isolates, pairs and groups > 2
#   num.groups <- max(partition)
#   num.nodes <- length(partition)
#   sizes <- table(partition)
#   isolates <- as.vector(which(sizes == 1))
#   pairs <- as.vector(which(sizes == 2))
#   others <- as.vector(which(sizes > 2))
#
#   aff <- as.matrix(table(data.frame(actor = 1:num.nodes, group = partition)))
#   net <- aff %*% t(aff)
#   net[upper.tri(net)] <- 0
#   diag(net) <- 0
#   paired.actors <- which(net == 1)
#   num.paired.actors <- sum(net)
#
#   groups.paired.actors <- rep(0, num.paired.actors)
#   pairs.paired.actors <- rep(0, num.paired.actors)
#
#   impossible.swaps <- 0
#   for (g in 1:num.groups) {
#     members <- which(partition == g)
#
#     if (g %in% pairs) {
#       i <- members[1]
#       j <- members[2]
#       p <- (i - 1 * num.nodes + j
#       pairs.paired.actors[which(paired.actors == p)] <- 1
#       groups.paired.actors[which(paired.actors == p)] <- g
#     }
#
#     if (g %in% others) {
#       for (i in 1:(length(members)-1)) {
#         for (j in (i + 1):length(members)) {
#           p <- (members[i]-1 * num.nodes + members[j]
#           groups.paired.actors[which(paired.actors == p)] <- g
#         }
#       }
#       num.pairs.g <- length(members)*(length(members)-1) * 0.5
#       impossible.swaps <- impossible.swaps + num.pairs.g * (num.pairs.g - 1) * 0.5
#     }
#
#   }
#
#   num.swaps <- num.paired.actors * (num.paired.actors - 1) / 2 - impossible.swaps - sum(pairs.paired.actors) * (sum(pairs.paired.actors)-1) / 2
#
#   return(list(isolates = isolates,
#               pairs = pairs,
#               others = others,
#               paired.actors = paired.actors,
#               groups.paired.actors = groups.paired.actors,
#               pairs.paired.actors = pairs.paired.actors,
#               num.swaps = num.swaps,
#               total = num.swaps))
#
# }
#
#
# sample_new_partition_p5 <- function(current.partition, size_neighborhood) {
#
#   # calculate the number of neighbor partitions
#   num.nodes <- length(current.partition)
#   num.groups <- max(current.partition)
#   size <- size_neighborhood
#   isolates <- size$isolates
#   pairs <- size$pairs
#   others <- size$others
#   paired.actors <- size$paired.actors
#   groups.paired.actors <- size$groups.paired.actors
#   pairs.paired.actors <- size$pairs.paired.actors
#   num.swaps <- size$num.swaps
#   total <- size$total
#
#   if (total == 0) return(current.partition)
#
#   new.partition <- current.partition
#
#   # decide which actors to swap
#   if (num.swaps > 0) {
#     found <- FALSE
#     while (!found) {
#       pick.1 <- sample(paired.actors, 1)
#       pick.2 <- sample(paired.actors, 1)
#       if (groups.paired.actors[which(paired.actors == pick.1)] != groups.paired.actors[which(paired.actors == pick.2)] &&
#          !(pairs.paired.actors[which(paired.actors == pick.1)] && pairs.paired.actors[which(paired.actors == pick.2)])) {
#         found <- TRUE
#       }
#     }
#
#     j1 <- ((pick.1 - 1) %/% num.nodes) + 1
#     i1 <- (pick.1 - 1) %% num.nodes + 1
#     j2 <- ((pick.2 - 1) %/% num.nodes) + 1
#     i2 <- (pick.2 - 1) %% num.nodes + 1
#     new.partition[i1] <- current.partition[i2]
#     new.partition[j1] <- current.partition[j2]
#     new.partition[i2] <- current.partition[i1]
#     new.partition[j2] <- current.partition[j1]
#     new.partition <- order_groupids(new.partition)
#   }
#
#
#   return(new.partition)
#
# }
#
#
# reachable_p5 <- function(partition1, partition2) {
#
#   # they are not neighbors if they have a different number of groups
#   if (max(partition1) != max(partition2)) return(FALSE)
#
#   num.nodes <- length(partition1)
#   check <- FALSE
#   i <- 1
#   j <- i + 1
#   k <- 1
#   l <- k + 1
#
#   # try remove all pairs of pairs
#   while (!check && i <= num.nodes) {
#
#     if (partition1[i] == partition1[j] && partition2[i] == partition2[j] &&
#        partition1[k] == partition1[l] && partition2[k] == partition2[l]) {
#       newp1 <- order_groupids(partition1[-c(i, j, k, l)])
#       newp2 <- order_groupids(partition2[-c(i, j, k, l)])
#
#       # check if the pairs were swapped
#       if (all(newp1 == newp2)) {
#         othersij1 <- which(partition1[-c(i, j)] == partition1[i])
#         othersij2 <- which(partition2[-c(i, j)] == partition2[i])
#         otherskl1 <- which(partition1[-c(k, l)] == partition1[k])
#         otherskl2 <- which(partition2[-c(k, l)] == partition2[k])
#         check <- (setequal(othersij1, otherskl2) && setequal(othersij2, otherskl1))
#         if (check) break
#       }
#     }
#
#     l <- l + 1
#     if (l > num.nodes) {
#       k <- k + 1
#       l <- k + 1
#       if (k > num.nodes) {
#         j <- j + 1
#         if (j > num.nodes) {
#           i <- i + 1
#           j <- i + 1
#         }
#       }
#     }
#   }
#
#   return(check)
# }
#
#
# ## NEIGHBORHOOD PI 6: pick two groups, re attribute nodes within them (careful, cannot be used alone)
#
# compute_size_neighborhood_p6 <- function(partition) {
#
#   # find isolates, pairs and groups > 2
#   num.groups <- max(partition)
#   num.nodes <- length(partition)
#   sizes <- table(partition)
#   isolates <- as.vector(which(sizes == 1))
#   pairs <- as.vector(which(sizes == 2))
#   others <- as.vector(which(sizes > 2))
#
#   group.pairs <- matrix(0, num.groups, num.groups)
#   nums.swaps <- matrix(0, num.groups, num.groups)
#
#   if (num.groups > 1) {
#
#     for (g1 in 1:(num.groups - 1)) {
#       for (g2 in (g1 + 1):num.groups) {
#         size.1 <- sum(partition == g1)
#         size.2 <- sum(partition == g2)
#
#         if (size.1 == 1 && size.2 == 1) {
#           nums.swaps[g1, g2] <- 0
#         } else {
#           k <- size.1
#           n <- size.1 + size.2
#           nums.swaps[g1, g2] <- factorial(n) / (factorial(k * factorial(n - k))
#           if (size.1 == size.2) nums.swaps[g1, g2] <- nums.swaps[g1, g2]  * 0.5
#           nums.swaps[g1, g2] <- nums.swaps[g1, g2] -1
#         }
#       }
#     }
#   }
#
#   return(list(isolates = isolates,
#               pairs = pairs,
#               others = others,
#               nums.swaps = nums.swaps,
#               num.swaps = sum(nums.swaps),
#               total = sum(nums.swaps)))
#
# }
#
#
# sample_new_partition_p6 <- function(current.partition, size_neighborhood) {
#
#   # calculate the number of neighbor partitions
#   num.nodes <- length(current.partition)
#   num.groups <- max(current.partition)
#   size <- size_neighborhood
#   isolates <- size$isolates
#   pairs <- size$pairs
#   others <- size$others
#   nums.swaps <- size$nums.swaps
#   num.swaps <- size$num.swaps
#   total <- size$total
#
#   if (total == 0) return(current.partition)
#
#   # decide which actors to swap
#   new.partition <- current.partition
#   pick <- sample(1:num.swaps, 1)
#   start <- 0
#   end <- 0
#
#   if (num.groups > 1) {
#     for (g1 in 1:(num.groups - 1)) {
#       for (g2 in (g1 + 1):num.groups) {
#
#         end <- end + nums.swaps[g1, g2]
#
#         if (start < pick && pick <= end) {
#           found <- FALSE
#           while (!found) {
#             tosample <- c(which(current.partition == g1), which(current.partition == g2))
#             new.g1 <- sample(tosample, sum(current.partition == g1))
#             new.g2 <- tosample[-which(tosample %in% new.g1)]
#             if (!(setequal(new.g1, current.partition[g1]) && setequal(new.g2, current.partition[g2])) &&
#                 !(setequal(new.g1, current.partition[g2]) && setequal(new.g2, current.partition[g1]))) found <- TRUE
#           }
#           new.partition[new.g1] <- g1
#           new.partition[new.g2] <- g2
#           new.partition <- order_groupids(new.partition)
#         }
#
#         start <- start + nums.swaps[g1, g2]
#       }
#     }
#   }
#
#   return(new.partition)
#
# }
#
#
# reachable_p6 <- function(partition1, partition2) {
#
#   # they are not neighbors if they have a different number of groups
#   if (max(partition1) != max(partition2)) return(FALSE)
#
#   num.nodes <- length(partition1)
#   num.groups1 <- max(partition1)
#   num.groups2 <- max(partition2)
#   check <- FALSE
#   g1 <- 1
#   g2 <- g1 + 1
#
#   # try remove groups (to check divisions) and pairs of groups (to check merges)
#   while (!check && g <= num.groups1) {
#
#     # first check whether g1 is not divided in partition2
#     membersg1 <- which(partition1 == g1)
#     newp1 <- order_groupids(partition1[-membersg1])
#     newp2 <- order_groupids(partition2[-membersg1])
#
#     if (all(newp1 == newp2)) {
#       check <- length(unique(partition2[membersg1])) == 2
#       if (check) break
#     }
#
#     # second check whether g1 and g2 are merged in partition2
#     membersg2 <- which(partition1 == g2)
#     newp1 <- order_groupids(partition1[-c(membersg1, membersg2)])
#     newp2 <- order_groupids(partition2[-c(membersg1, membersg2)])
#
#     if (all(newp1 == newp2)) {
#       check <- length(unique(partition2[c(membersg1, membersg2)])) == 1
#       if (check) break
#     }
#
#     g2 <- g2 + 1
#     if (g2 > num.nodes) {
#       g1 <- g1 + 1
#       g2 <- g2 + 1
#     }
#   }
#
#   return(check)
# }
#
# ## NEIGHBORHOOD PI 7: take two nodes together in partition 2, bring them together (single move)
# ## otherwise separate then
# ## careful need isolates
#
# compute_size_neighborhood_p7 <- function(partition, partition2) {
#
#   num.nodes <- length(partition)
#   present <- !is.na(partition)
#   present2 <- !is.na(partition2)
#
#   num.groups <- max(partition, na.rm = TRUE)
#   sizes <- table(partition)
#   isolates <- as.vector(which(sizes == 1))
#   pairs <- as.vector(which(sizes == 2))
#
#   affiliation <- as.matrix(table(data.frame(actor = 1:num.nodes, group = partition)))
#   affiliation2 <- as.matrix(table(data.frame(actor = 1:num.nodes, group = partition2)))
#   adjacency <- affiliation %*% t(affiliation)
#   adjacency2 <- affiliation2 %*% t(affiliation2)
#   diag(adjacency) <- 0
#   diag(adjacency2) <- 0
#
#   nums.swaps <- matrix(0, num.nodes, num.nodes)
#   for (i in 1:num.nodes) {
#     for (j in 1:num.nodes) {
#       ad <- adjacency[i, j]
#       ad2 <- adjacency2[i, j]
#       if (!is.na(ad) && !is.na(ad2) && ad2 == 1) {
#         if (ad == 0) {  # if separated join if possible (two isolates will only join once)
#           if (partition[i] %in% isolates && partition[j] %in% isolates && j < i) {
#             nums.swaps[i, j] <- 0
#           } else {
#             nums.swaps[i, j] <- 1
#           }
#         }
#         if (ad == 1) {  # if together separate if possible (a pair will only separate once)
#           if (partition[i] %in% pairs && j < i) {
#             nums.swaps[i, j] <- num.groups - 1
#           } else {
#             nums.swaps[i, j] <- num.groups
#           }
#
#         }
#       }
#     }
#   }
#
#   return(list(adjacency = adjacency,
#               pairs = pairs,
#               nums.swaps = nums.swaps,
#               num.swaps = sum(nums.swaps),
#               total = sum(nums.swaps)))
# }
#
#
# sample_new_partition_p7 <- function(current.partition, current.partition2, size_neighborhood) {
#
#   num.nodes <- length(current.partition)
#   num.groups <- max(current.partition, na.rm = TRUE)
#   adjacency <- size_neighborhood$adjacency
#   pairs <- size_neighborhood$pairs
#   nums.swaps <- size_neighborhood$nums.swaps
#   num.swaps <- size_neighborhood$num.swaps
#   total <- size$total
#
#   if (total == 0) return(current.partition1)
#
#   present <- !(is.na(current.partition))
#   new.partition <- current.partition
#
#   probas_pairs <- as.vector(nums.swaps) / num.swaps
#   pick <- sample(1:(num.nodes * num.nodes), 1, prob = probas_pairs)
#   i <- (pick - 1) %% num.nodes + 1
#   j <- ((pick - 1) %/% num.nodes) + 1
#
#   if (adjacency[i, j] == 0) {
#     new.partition[i] <- current.partition[j]
#   } else {
#     if (current.partition[i] %in% pairs && j < i) {
#       tosample <- (1:num.groups)[(1:num.groups) != current.partition[i]]
#       if (length(tosample) == 1) pickgroup <- tosample
#       else pickgroup <- sample(tosample, 1)
#       if (pickgroup == current.partition[i]) new.partition[i] <- num.groups + 1
#       else new.partition[i] <- pickgroup
#     }else {
#       tosample <- 1:num.groups
#       if (length(tosample) == 1) pickgroup <- tosample
#       else pickgroup <- sample(tosample, 1)
#       if (pickgroup == current.partition[i]) new.partition[i] <- num.groups + 1
#       else new.partition[i] <- pickgroup
#     }
#   }
#   new.partition[present] <- order_groupids(new.partition[present])
#
#   return(new.partition)
#
# }
#
#
#
# ## NEIGHBORHOOD PI 4 RESTRICTED: swaps a pair of nodes
# # WARNING!!!!: for now it only works if sizes allowed are an interval ([size_min, size_max])
# compute_size_neighborhood_p4_restricted <- function(partition, sizes.simulated) {
#
#   # check if current partition is allowed
#   if (!check_sizes(partition, sizes.simulated)) {
#     print(partition)
#     #stop("The partition we are in is not allowed.")
#   }
#
#   # find maximum size allowed
#   smax <- max(sizes.simulated)
#   smin <- min(sizes.simulated)
#
#   # find isolates, pairs and groups > 2
#   num.groups <- max(partition)
#   num.nodes <- length(partition)
#   sizes <- table(partition)
#   isolates <- as.vector(which(sizes == 1))
#   pairs <- as.vector(which(sizes == 2))
#   others <- as.vector(which(sizes > 2))
#   abovemin_groups <- as.vector(which(sizes > (smin + 1)))
#   belowmax_groups <- as.vector(which(sizes < (smax - 1)))
#
#   nums.swaps <- rep(0, max(partition))
#   done.pairs <- rep(0, length(pairs))
#   for (g in 1:num.groups) {
#
#     if (g %in% isolates) {
#       # nothing happens, there is no pair to swap here
#     }
#     if (g %in% pairs) {
#       # then the pair can move to a group <smax - 1, unless it's another pair and the combination of the two pairs is already counted
#       # is 2 <smax - 1 we have to remove the pair from the count of possible other groups
#       if (2 < smax - 1) nums.swaps[g] <- length(belowmax_groups) - 1 - sum(done.pairs) + done.pairs[pairs == g]
#       if (2 >= smax - 1) nums.swaps[g] <- length(belowmax_groups) # then it's only isolates
#       done.pairs[pairs == g] <- 1
#     }
#     if (g %in% others && g %in% abovemin_groups) {
#       # then any pair inside the group can be separated and added to an available group that is <smax - 1, unless the group is not >smin + 1
#       sizeg <- sum(partition == g)
#       gbelowmax <- (g %in% belowmax_groups)
#       # if it can create its own isolated pair
#       if (gbelowmax && 2 %in% sizes.simulated) nums.swaps[g] <- length(belowmax_groups) * sizeg * (sizeg - 1) / 2
#       if (!gbelowmax && 2 %in% sizes.simulated) nums.swaps[g] <- (length(belowmax_groups) + 1) * sizeg * (sizeg - 1) / 2
#       # if not
#       if (gbelowmax && !(2 %in% sizes.simulated)) nums.swaps[g] <- (length(belowmax_groups)-1) * sizeg * (sizeg - 1) / 2
#       if (!gbelowmax && !(2 %in% sizes.simulated)) nums.swaps[g] <- length(belowmax_groups) * sizeg * (sizeg - 1) / 2
#     } # we can never break a group that is already of minimal size
#   }
#
#   return(list(isolates = isolates,
#               pairs = pairs,
#               others = others,
#               abovemin_groups = abovemin_groups,
#               belowmax_groups = belowmax_groups,
#               nums.swaps = nums.swaps,
#               num.swaps = sum(nums.swaps),
#               total = sum(nums.swaps)))
#
# }
#
# sample_new_partition_p4_restricted <- function(current.partition, size_neighborhood, sizes.simulated) {
#
#   # calculate the number of neighbor partitions
#   num.nodes <- length(current.partition)
#   num.groups <- max(current.partition)
#   size <- size_neighborhood
#   isolates <- size$isolates
#   pairs <- size$pairs
#   others <- size$others
#   abovemin_groups <- size$abovemin_groups
#   belowmax_groups <- size$belowmax_groups
#   nums.swaps <- size$nums.swaps
#   num.swaps <- size$num.swaps
#   total <- size$total
#
#   # allowed sizes
#   smax <- max(sizes.simulated)
#   smin <- min(sizes.simulated)
#
#   pick.1 <- sample(total, 1)
#
#   # decide which group to pick a pair from
#   all.groups <- 1:num.groups
#   new.partition <- current.partition
#   done.pairs <- rep(0, length(pairs))
#
#   for (g in 1:num.groups) {
#
#     if (g == 1) {start <- 0} else {start <- sum(nums.swaps[1:g - 1])}
#     if (g == num.groups) {end <- num.swaps} else {end <- sum(nums.swaps[1:g])}
#
#     if (pick.1 > start && pick.1 <= end) {
#
#       # a pair can join other available groups (<smax - 1), unless it's a pair that has already been "joined"
#       if (g %in% pairs) {
#         members <- which(current.partition == g)
#         previous.pairs <- pairs[done.pairs == 1]
#         if (2 < smax - 1) {
#           tosample <- all.groups[(all.groups %in% belowmax_groups) & (all.groups != g) & !(all.groups %in% previous.pairs)]
#           if (length(tosample) == 1) newg <- tosample
#           if (length(tosample) >= 2) newg <- sample(tosample, 1)
#         }else {
#           tosample <- all.groups[(all.groups %in% belowmax_groups)]
#           if (length(tosample) == 1) newg <- tosample
#           if (length(tosample) >= 2) newg <- sample(tosample, 1)
#         }
#         new.partition[members] <- newg
#         new.partition <- order_groupids(new.partition)
#       }
#
#       # a pair in a bigger group (>smin + 1) can join any other group (<smax - 1) or create its isolate
#       if (g %in% others && g %in% abovemin_groups) {
#         members <- which(current.partition == g)
#         sizeg <- length(members)
#         pairg <- sample(members, 2)
#         if (2 %in% sizes.simulated) tosample <- all.groups[(all.groups %in% belowmax_groups) || (all.groups == g)]
#         else tosample <- all.groups[(all.groups %in% belowmax_groups) & (all.groups != g)]
#         if (length(tosample) == 1) newg <- tosample
#         if (length(tosample) >= 2) newg <- sample(tosample, 1)
#         if (newg == g) newg <- num.groups + 1
#         new.partition[pairg] <- newg
#         new.partition <- order_groupids(new.partition)
#       }
#
#     }
#
#     if (g %in% pairs) done.pairs[pairs == g] <- 1
#   }
#
#   return(new.partition)
#
# }
#
#
# ## NEIGHBORHOOD PI 5 RESTRICTED: only swaps of two pairs of nodes (careful, cannot be used alone)
# compute_size_neighborhood_p5_restricted <- function(partition, sizes.simulated) {
#   # check if current partition is allowed
#   if (!check_sizes(partition, sizes.simulated)) stop("The partition we are in is not allowed.")
#   return(compute_size_neighborhood_p5(partition))
# }
#
# sample_new_partition_p5_restricted <- function(current.partition, size_neighborhood, sizes.simulated) {
#   return(sample_new_partition_p5(current.partition, size_neighborhood))
# }
#
#
# ## NEIGHBORHOOD PI 6 RESTRICTED: exchanges of two groups (careful, cannot be used alone)
# compute_size_neighborhood_p6_restricted <- function(partition, sizes.simulated) {
#   # check if current partition is allowed
#   if (!check_sizes(partition, sizes.simulated)) stop("The partition we are in is not allowed.")
#   return(compute_size_neighborhood_p6(partition))
# }
#
# sample_new_partition_p6_restricted <- function(current.partition, size_neighborhood, sizes.simulated) {
#   return(sample_new_partition_p6(current.partition, size_neighborhood))
# }
#
#
# ## NEIGHBORHOOD PI 7 RESTRICTED: take two nodes together in partition 2, bring them together (single move)
# ## otherwise separate then
# ## careful need isolates
# compute_size_neighborhood_p7_restricted <- function(partition, partition2, sizes.simulated) {
#
#   # find maximum size allowed
#   smax <- max(sizes.simulated)
#   smin <- min(sizes.simulated)
#
#   num.nodes <- length(partition)
#   present <- !is.na(partition)
#   present2 <- !is.na(partition2)
#
#   num.groups <- max(partition, na.rm = TRUE)
#   sizes <- table(partition)
#   isolates <- as.vector(which(sizes == 1))
#   pairs <- as.vector(which(sizes == 2))
#   abovemin_groups <- as.vector(which(sizes > smin))
#   belowmax_groups <- as.vector(which(sizes < smax))
#   max_groups <- as.vector(which(sizes == smax))
#
#   affiliation <- as.matrix(table(data.frame(actor = 1:num.nodes, group = partition)))
#   affiliation2 <- as.matrix(table(data.frame(actor = 1:num.nodes, group = partition2)))
#   adjacency <- affiliation %*% t(affiliation)
#   adjacency2 <- affiliation2 %*% t(affiliation2)
#   diag(adjacency) <- 0
#   diag(adjacency2) <- 0
#
#   nums.swaps <- matrix(0, num.nodes, num.nodes)
#   for (i in 1:num.nodes) {
#     for (j in 1:num.nodes) {
#       ad <- adjacency[i, j]
#       ad2 <- adjacency2[i, j]
#       if (!is.na(ad) && !is.na(ad2) && ad2 == 1) {
#         if (ad == 0) { # if separated join if possible (two isolates will only join once)
#           if (partition[i] %in% isolates && partition[j] %in% isolates && j < i) {
#             nums.swaps[i, j] <- 0
#           } else {
#             nums.swaps[i, j] <- (partition[i] %in% c(abovemin_groups, isolates)) & (partition[j] %in% belowmax_groups)
#           }
#         }
#         if (ad == 1) { # if together separate if possible (a pair will only separate once)
#           if (partition[i] %in% abovemin_groups) {
#             if (partition[i] %in% pairs && j < i) {
#               nums.swaps[i, j] <- num.groups - length(max_groups) + (partition[i] %in% max_groups) - !(1 %in% sizes.simulated) - 1
#             } else {
#               nums.swaps[i, j] <- num.groups - length(max_groups) + (partition[i] %in% max_groups) - !(1 %in% sizes.simulated)
#             }
#           }
#         }
#       }
#     }
#   }
#
#   return(list(abovemin_groups = abovemin_groups,
#               belowmax_groups = belowmax_groups,
#               adjacency = adjacency,
#               pairs = pairs,
#               nums.swaps = nums.swaps,
#               num.swaps = sum(nums.swaps),
#               total = sum(nums.swaps)))
# }
#
# sample_new_partition_p7_restricted <- function(current.partition, size_neighborhood, current.partition2, sizes.simulated) {
#
#   # find maximum size allowed
#   smax <- max(sizes.simulated)
#   smin <- min(sizes.simulated)
#
#   num.nodes <- length(current.partition)
#   num.groups <- max(current.partition, na.rm = TRUE)
#   abovemin_groups <- size_neighborhood$abovemin_groups
#   belowmax_groups <- size_neighborhood$belowmax_groups
#   adjacency <- size_neighborhood$adjacency
#   pairs <- size_neighborhood$pairs
#   nums.swaps <- size_neighborhood$nums.swaps
#   num.swaps <- size_neighborhood$num.swaps
#   total <- size$total
#
#   if (total == 0) return(current.partition)
#
#   present <- !(is.na(current.partition))
#   new.partition <- current.partition
#
#   probas_pairs <- as.vector(nums.swaps) / num.swaps
#   pick <- sample(1:(num.nodes * num.nodes), 1, prob = probas_pairs)
#   i <- (pick - 1) %% num.nodes + 1
#   j <- ((pick - 1) %/% num.nodes) + 1
#
#   if (adjacency[i, j] == 0) {
#     new.partition[i] <- current.partition[j]
#   }
#   if (adjacency[i, j] == 1) {
#     if (current.partition[i] %in% pairs && j < i) {
#       tosample <- belowmax_groups[belowmax_groups != current.partition[i]]
#     } else {
#       if (1 %in% sizes.simulated) tosample <- unique(c(belowmax_groups, current.partition[i]))
#       else tosample <- belowmax_groups[belowmax_groups != current.partition[i]]
#     }
#     if (length(tosample) == 1) pickgroup <- tosample
#     else pickgroup <- sample(tosample, 1)
#     if (pickgroup == current.partition[i]) new.partition[i] <- num.groups + 1
#     else new.partition[i] <- pickgroup
#   }
#   new.partition[present] <- order_groupids(new.partition[present])
#
#   return(new.partition)
#
# }
#
#



##############################################################################

### Deprecated 2: neighborhoods for the multiple partition estimation

# ## NEIGHBORHOOD PI 7: pick two partitions and exchange then (only for multiple partition estimation)
# # since we force one change, and it's reversible (bijection), we don't have to care about the size of the neighborhood
# # CAREFUL, only works when singletons are allowed
# compute_size_neighborhood_p7 <- function(partition1, partition2) {
#
#   num.nodes <- length(partition1)
#   present1 <- !is.na(partition1)
#   present2 <- !is.na(partition2)
#
#   # calculate the number of ways to reconstruct p1 after having done the exchange with p2
#   p1withoutabsentp2 <- partition1[present1 & present2]
#   max1 <- length(unique(p1withoutabsentp2))
#   nodesunknown1 <- sum(present1 & !present2)
#   if (nodesunknown1 == 0) {
#     possibilities1 <- 1
#   } else {
#     possibilities1 <- 0
#     for (k in 1:nodesunknown1) {
#       possibilities1 <- possibilities1 + Stirling(nodesunknown1, k)*(max1 + 1)^k
#     }
#   }
#
#   # same the other way around
#   p2withoutabsentp1 <- partition2[present2 & present1]
#   max2 <- length(unique(p2withoutabsentp1))
#   nodesunknown2 <- sum(present2 & !present1)
#   if (nodesunknown2 == 0) {
#     possibilities2 <- 1
#   } else {
#     possibilities2 <- 0
#     for (k in 1:nodesunknown2) {
#       possibilities2 <- possibilities2 + Stirling(nodesunknown2, k)*(max2 + 1)^k
#     }
#   }
#
#   return(possibilities1 * possibilities2)
# }
#
# sample_new_partition_p7 <- function(current.partition1, current.partition2) {
#
#   num.nodes <- length(current.partition1)
#   present1 <- which(!is.na(current.partition1))
#   present2 <- which(!is.na(current.partition2))
#
#   new.partition1 <- rep(0, num.nodes)
#   new.partition2 <- rep(0, num.nodes)
#
#   # fill in the first with the second, and create isolates for people not present in the second
#   new.partition1[present1] <- current.partition2[present1]
#   max1 <- max(new.partition1, na.rm = TRUE)
#   singletons1 <- which(is.na(new.partition1))
#   if (length(singletons1) > 0) { # assign randomly the ones who are not present in the other one
#     for (s in singletons1) {
#       newgroup <- sample(1:(max1 + 1), 1)
#       new.partition1[s] <- newgroup
#       if (newgroup == (max1 + 1)) max1 <- max1 + 1
#     }
#   }
#   new.partition1[present1] <- order_groupids(new.partition1[present1])
#   new.partition1[new.partition1 == 0] <- NA
#
#   # same the other way around
#   new.partition2[present2] <- current.partition1[present2]
#   max2 <- max(new.partition2, na.rm = TRUE)
#   singletons2 <- which(is.na(new.partition2))
#   if (length(singletons2) > 0) { # assign randomly the ones who are not present in the other one
#     for (s in singletons2) {
#       newgroup <- sample(1:(max2 + 1), 1)
#       new.partition2[s] <- newgroup
#       if (newgroup == (max2 + 1)) max2 <- max2 + 1
#     }
#   }
#   new.partition2[present2] <- order_groupids(new.partition2[present2])
#   new.partition2[new.partition2 == 0] <- NA
#
#   return(list(new.partition1 = new.partition1,
#               new.partition2 = new.partition2))
#
# }
#
#
# ## NEIGHBORHOOD PI 8: pick another partition, a group in it, and reproduce it in the current partition, or reorganizes it
# # CAREFUL, only works when singletons are allowed
# compute_size_neighborhood_p8 <- function(partition1, partition2) {
#
#   num.nodes <- length(partition1)
#   present1 <- !is.na(partition1)
#   present2 <- !is.na(partition2)
#
#   # count number of groups from the second partition
#   groups.toreorganize <- unique(partition2[present1 & present2])
#   nums.reorganizations <- rep(0, length(groups.toreorganize))
#
#   # count for each of them number of moves
#   for (g in 1:length(groups.toreorganize)) {
#     group <- groups.toreorganize[g]
#     members <- which(partition2 == group)
#     max1 <- length(unique(partition1[present1 & partition1 != group]))
#     for (k in 1:length(members)) {
#       nums.reorganizations[g] <- nums.reorganizations[g] + Stirling(length(members), k)*(max1 + 1)^k
#     }
#     nums.reorganizations[g] <- nums.reorganizations[g] + 1 # when we just copy the group
#   }
#
#   return(list(groups.toreorganize = groups.toreorganize,
#               nums.reorganizations = nums.reorganizations,
#               num.reorganizations = sum(nums.reorganizations)))
# }
#
# sample_new_partition_p8 <- function(current.partition1, current.partition2, size_neighborhood) {
#
#   groups.toreorganize <- size_neighborhood$groups.toreorganize
#   nums.reorganizations <- size_neighborhood$nums.reorganizations
#   num.reorganizations <- size_neighborhood$num.reorganizations
#
#   num.nodes <- length(current.partition1)
#   present1 <- !is.na(current.partition1)
#   present2 <- !is.na(current.partition2)
#
#   new.partition <- current.partition1
#
#   # pick
#   probas_groups <- nums.reorganizations / num.reorganizations
#   pick <- sample(c(0, groups.toreorganize), 1, prob = c(1 /  num.reorganizations, probas_groups))
#
#   #randomize group
#   if (pick > 0) {
#     group <- pick
#     members <- which(current.partition2 == group && present1)
#     max1 <- length(unique(current.partition1[present1 & current.partition1 != group]))
#     for (a in members) {
#       newgroup <- sample(1:(max1 + 1), 1)
#       new.partition[a] <- newgroup
#       if (newgroup == (max1 + 1)) max1 <- max1 + 1
#     }
#     new.partition[present1] <- order_groupids(new.partition[present1])
#   }
#
#   # copy group
#   if (pick ==0) {
#     # reproduce the group
#     new.partition[members] <- max(new.partition, na.rm = TRUE) + 1
#     #reorder
#     new.partition[present1] <- order_groupids(new.partition[present1])
#   }
#
#   return(new.partition)
#
# }
#
# #sample_new_partition_p8 <- function(current.partition1, current.partition2) {
# #
# #  num.nodes <- length(current.partition1)
# #  present1 <- !is.na(current.partition1)
# #  present2 <- !is.na(current.partition2)
# #
# #  new.partition <- current.partition1
# #
# #  # pick a group to reproduce
# #  topick <- unique(current.partition2[present1 & present2])
# #  g <- sample(topick, 1)
# #  members <- which(current.partition2 == g && present1)
# #
# #  # reproduce the group
# #  new.partition[members] <- max(new.partition, na.rm = TRUE) + 1
# #
# #  #reorder
# #  new.partition[present1] <- order_groupids(new.partition[present1])
# #
# #  return(new.partition)
# #
# #}
#
# ## NEIGHBORHOOD PI 8: pick another partition, a group in it, and reproduce it in the current partition, or reorganizes it
# # CAREFUL, only works when singletons are allowed
# compute_size_neighborhood_p8_restricted <- function(partition1, partition2, sizes.simulated) {
#
#   num.nodes <- length(partition1)
#   present1 <- !is.na(partition1)
#   present2 <- !is.na(partition2)
#   smax <- max(sizes.simulated)
#
#   # count number of groups from the second partition
#   groups.toreorganize <- unique(partition2[present1 & present2])
#   combinations.pergroup <- list()
#   nums.partitions.pergroup <- list()
#   nums.reorganizations.pergroup <- list()
#   nums.reorganizations <- rep(0, length(groups.toreorganize))
#
#   # count for each of them number of moves
#   for (g in 1:length(groups.toreorganize)) {
#     group <- groups.toreorganize[g]
#     members <- which(partition2 == group)
#     members <- which(1:num.nodes %in% members & present1)
#
#     # if the group is not present, we will create it
#     if (length(unique(partition1[members])) == 1 && sum(partition1[-members] == (partition1[members])[1], na.rm = TRUE) == 0) {
#
#       nums.reorganizations.pergroup[[g]] <- NULL
#       nums.reorganizations[g] <- 1
#
#       # if the group is present, we will reshuffle
#     } else {
#       nmembers <- length(members)
#       sizes <- table(partition1[-members])
#
#       # find ways to divide the group
#       combs <- parts(nmembers)
#       combinations.pergroup[[g]] <- combs
#
#       nums.p <- rep(1, ncol(combs))
#       nums.r <- rep(1, ncol(combs))
#       for (c in 1:ncol(combs)) {
#
#         # find number partitions of each
#         ncpt <- nmembers
#         for (s in 1:nrow(combs)) {
#           ncs <- sum(combs[, c] == s)
#           nums.p[c] <- nums.p[c] * (choose(ncpt, s))^ncs / factorial(ncs)
#           ncpt <- ncpt - ncs * s
#         }
#
#         # find ways to reorganize
#         for (s in 1:nrow(combs)) {
#           ncs <- sum(combs[, c] == s)
#           ncpt <- sum(sizes <= (smax - s))
#           for (k in 1:ncs) {
#             nums.r[c] <- nums.r[c] * ncpt
#             ncpt <- ncpt - 1
#           }
#         }
#       }
#
#       nums.partitions.pergroup[[g]] <- nums.p
#       nums.reorganizations.pergroup[[g]] <- nums.r
#       nums.reorganizations[g] <- nums.p * nums.r
#     }
#
#   }
#
#   return(list(groups.toreorganize = groups.toreorganize,
#               combinations.pergroup = combinations.pergroup,
#               nums.partitions.pergroup = nums.partitions.pergroup,
#               nums.reorganizations.pergroup = nums.reorganizations.pergroup,
#               nums.reorganizations = nums.reorganizations,
#               num.reorganizations = sum(nums.reorganizations)))
# }
#
# sample_new_partition_p8_restricted <- function(current.partition1, current.partition2, size_neighborhood, sizes.simulated) {
#
#   groups.toreorganize <- size_neighborhood$groups.toreorganize
#   combinations.pergroup <- size_neighborhood$combinations.pergroup
#   nums.partitions.pergroup <- size_neighborhood$nums.partitions.pergroup
#   nums.reorganizations.pergroup <- size_neighborhood$nums.reorganizations.pergroup
#   nums.reorganizations <- size_neighborhood$nums.reorganizations
#   num.reorganizations <- size_neighborhood$num.reorganizations
#
#   num.nodes <- length(current.partition1)
#   present1 <- !is.na(current.partition1)
#   present2 <- !is.na(current.partition2)
#
#   new.partition <- current.partition1
#
#   # pick
#   probas_groups <- nums.reorganizations / num.reorganizations
#   group <- sample(groups.toreorganize, 1, prob = probas_groups)
#   members <- which(current.partition2 == group)
#   members <- which(1:num.nodes %in% members & present1)
#
#   # if the group is not present, we will create it
#   if (length(unique(current.partition1[members])) == 1 && sum(current.partition1[-members] == (current.partition1[members])[1], na.rm = TRUE) == 0) {
#
#     new.partition[members] <- max(new.partition, na.rm = TRUE) + 1
#     new.partition[present1] <- order_groupids(new.partition[present1])
#
#     # if the group is present, we will reshuffle
#   } else {
#
#     nmembers <- length(members)
#     gleft <- unique(current.partition1[-members])
#     gleft <- gleft[!is.na(gleft)]
#     sizes <- apply(gleft, FUN = function(x) {return(sum(current.partition1[-members] == x))})
#
#     g <- which(groups.toreorganize == group)
#     combs <- combinations.pergroup[[g]]
#
#     # sample a combination
#     probas_combs <- nums.partitions.pergroup[[g]] * nums.reorganizations.pergroup[[g]] / nums.reorganizations[g]
#     c <- sample(1:ncol(combs), 1, prob = probas_combs)
#
#     # create a combination
#     ntaken <- rep(FALSE, nmembers)
#     gtaken <- rep(FALSE, length(gleft))
#     for (s in nrow(combs)) {
#       ncs <- sum(combs[, c] == s)
#       for (k in 1:ncs) {
#         ns <- choose(members[!ntaken], s)
#         topick <- gleft[sizes <= (smax - s) & !gtaken]
#         newg <- sample(c(0, topick), 1)
#         if (newg == 0) {
#           new.partition[ns] <- max(new.partition, na.rm = TRUE) + 1
#         } else {
#           new.partition[ns] <- newg
#           gtaken[gleft == newg] <- TRUE
#         }
#         ntaken[members == ns] <- TRUE
#       }
#     }
#     newnodes <- sample(1:nmembers, nmembers)
#
#     new.partition[present1] <- order_groupids(new.partition[present1])
#
#   }
#
#   return(new.partition)
#
# }
