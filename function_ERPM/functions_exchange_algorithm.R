######################################################################
## Simulation and estimation of Exponential Random Partition Models ## 
## Functions used to run the exhange algorithm estimation (Bayesian)##
## Author: Marion Hoffman                                           ##
######################################################################


# Wrapper function for the exchange algorithm (bayesian estimation) for multiple partitions 
# (without any hidden network)
draw_exchangealgorithm_multiple <- function(partitions, 
                                            z.obs,
                                            presence.tables, 
                                            nodes, 
                                            objects, 
                                            effects, 
                                            mean.priors,
                                            sd.priors, 
                                            start.chains,
                                            burnin.1,
                                            thining.1,
                                            num.chains,
                                            length.chains,
                                            mini.steps.2,
                                            burnin.2, 
                                            neighborhood.partition,
                                            neighborhood.augmentation,
                                            sizes.allowed,
                                            sizes.simulated,
                                            parallel,
                                            cpus) {
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)
  num.obs <- ncol(presence.tables)
  
  # create the different chains to simulate, parallelize if possible
  if(parallel && cpus >= num.chains){
  
    min_cpu_per_chain <- cpus %/% num.chains
    cpus_per_chain <- rep(min_cpu_per_chain,num.chains)  
    remaining_cpus <- cpus %% num.chains
    cpus_per_chain[1:remaining_cpus] <- cpus_per_chain[1:remaining_cpus]+1
    chains_per_cpu <- c()
    for(c in 1:num.chains) chains_per_cpu <- c(chains_per_cpu,rep(c,cpus_per_chain[c]))
    
    sfExport("partitions", "z.obs","presence.tables","nodes", "objects","effects", "mean.priors","sd.priors", "start.chains","burnin.1","thining.1","length.chains","mini.steps.2","burnin.2", "neighborhood.partition","neighborhood.augmentation","sizes.allowed","sizes.simulated","chains_per_cpu")
    all.chains <- sfLapply(1:cpus, fun = function(k) {
      set.seed(k)
      c <- chains_per_cpu[k]
      chain <- exchangealgorithm_multiple(partitions, 
                                          z.obs,
                                          presence.tables, 
                                          nodes, 
                                          objects, 
                                          effects, 
                                          mean.priors,
                                          sd.priors, 
                                          start.chains[[c]],
                                          burnin.1,
                                          thining.1,
                                          length.chains,
                                          mini.steps.2,
                                          burnin.2, 
                                          neighborhood.partition,
                                          neighborhood.augmentation,
                                          sizes.allowed,
                                          sizes.simulated)
      return(chain)
    }
    )
    
  }else{
   
    all.chains <- list()
    for(c in 1:num.chains){
      set.seed(c)
      all.chains[[c]] <- exchangealgorithm_multiple(partitions, 
                                                    z.obs,
                                                    presence.tables, 
                                                    nodes, 
                                                    objects, 
                                                    effects, 
                                                    mean.priors,
                                                    sd.priors, 
                                                    start.chains[[c]],
                                                    burnin.1,
                                                    thining.1,
                                                    length.chains,
                                                    mini.steps.2,
                                                    burnin.2, 
                                                    neighborhood.partition,
                                                    neighborhood.augmentation,
                                                    sizes.allowed,
                                                    sizes.simulated)
    }
    
  }
  
  # Posterior means and standard deviations
  all.draws <- do.call(rbind,all.chains)
  post.mean <- apply(all.draws, 2, mean)
  post.sd <- apply(all.draws, 2, sd)
  
  
  return( list("post.mean" = post.mean, 
               "post.sd" = post.sd,
               "all.chains" = all.chains) )
  
}


# Core function for the exchange algorithm (bayesian estimation) for multiple partitions 
# (without any hidden network)
exchangealgorithm_multiple <- function(partitions, 
                                       z.obs,
                                       presence.tables, 
                                       nodes, 
                                       objects, 
                                       effects, 
                                       mean.priors,
                                       sd.priors, 
                                       start.chain,
                                       burnin.1,
                                       thining.1,
                                       length.chains,
                                       mini.steps.2,
                                       burnin.2, 
                                       neighborhood.partition,
                                       neighborhood.augmentation,
                                       sizes.allowed,
                                       sizes.simulated) {
  
  num.effects <- length(effects$names)
  num.nodes <- nrow(nodes)
  num.obs <- ncol(presence.tables)

  # starts of the chain
  augmented.partitions <- partitions
  current.theta <- start.chain
  new.theta <- current.theta
  
  # store results
  all.theta <- c()
  
  # go through the chain and sample
  cpt <- 0
  cpt_thining <- 0
  end.walk <- F
  
  while(!end.walk){
    
    # draw new theta augmented, and new partitions
    p <- sample(num.effects,1)
    new.theta[p] <- rnorm(1, mean=current.theta[p], sd=neighborhood.augmentation[p])  
    partition.draw <- draw_Metropolis_multiple(new.theta, 
                                         first.partitions = partitions, 
                                         presence.tables, 
                                         nodes, 
                                         effects, 
                                         objects, 
                                         burnin = burnin.2, 
                                         thining = 1, 
                                         num.steps = 1, 
                                         mini.steps = mini.steps.2, 
                                         neighborhood = neighborhood.partition, 
                                         sizes.allowed, 
                                         sizes.simulated,
                                         return.all.partitions = F) 
    augmented.partitions <- partition.draw$last.partitions
    z.augmented <- computeStatistics_multiple(augmented.partitions, presence.tables, nodes, effects, objects)
    
    # compute acceptance ratio
    ratio.probas <- prod(pnorm(new.theta, mean=mean.priors, sd=sd.priors) /
                         pnorm(current.theta, mean=mean.priors, sd=sd.priors))
    ratio.obs <- exp(sum(new.theta*z.obs) - 
                     sum(current.theta*z.obs)) 
    ratio.augmented <- exp(sum(current.theta*z.augmented) - 
                           sum(new.theta*z.augmented)) 
    
    acceptance.ratio <- ratio.obs * ratio.probas * ratio.augmented
    
    #proba.current.theta <- prod(pnorm(current.theta, mean=mean.priors, sd=sd.priors))
    #proba.new.theta <- prod(pnorm(new.theta, mean=mean.priors, sd=sd.priors))
    #denom.current.theta.obs <- exp(sum(current.theta*z.obs))
    #denom.current.theta.augmented <- exp(sum(current.theta*z.augmented))
    #denom.new.theta.obs <- exp(sum(new.theta*z.obs))
    #denom.new.theta.augmented <- exp(sum(new.theta*z.augmented))
    
    #acceptance.ratio <- (denom.new.theta.obs/denom.current.theta.obs) * 
    #  (proba.new.theta/proba.current.theta) * 
    #  (denom.current.theta.augmented/denom.new.theta.augmented)  
    
    # make the update to new theta if needed
    proba.change <- min(1,acceptance.ratio)
    if(is.nan(proba.change)) proba.change <- 0
    change.made <- (runif(1,0,1) <= proba.change)
    
    if(change.made){
      current.theta <- new.theta
    }
    
    # move through the chain
    cpt <- cpt + 1
    if(cpt > burnin.1) cpt_thining <- cpt_thining + 1
    
    # store the results if we are out of burnin
    if(cpt >= burnin.1 && cpt_thining == thining.1) {
      all.theta <- rbind(all.theta, current.theta)
      cpt_thining <- 0
    }
      
    # stop the walk if number of steps reached 
    end.walk <- (cpt >= (burnin.1 + thining.1 * length.chains))
      
  }

  return(all.theta)
}