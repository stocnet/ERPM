run_phase3 <- function(estimates.phase2, 
                       z.obs, 
                       nodes, 
                       effects, 
                       objects, 
                       burnin, 
                       thining,
                       mini.steps, 
                       length.p3, 
                       neighborhood,
                       sizes.allowed,
                       sizes.simulated,
                       fixed.estimates) {
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)
  # TODO: check this
  # length.p3 <- 1000
  
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
  
  # simulate a large sample with the estimates found in phase 2 
  results.phase3 <- draw_Metropolis(estimates.phase2, first.partition, nodes, effects, objects, burnin, thining, length.p3, mini.steps, neighborhood, sizes.allowed, sizes.simulated)
  z.phase3 <- results.phase3$draws
  
  # calculate covariance and scaling
  inverted_matrices <- calculate_inverted_covariance_and_scaling(estimates.phase2, 
                                                        z.obs, 
                                                        nodes, 
                                                        effects, 
                                                        objects, 
                                                        length.phase = length.p3, 
                                                        z.phase = z.phase3,
                                                        fixed.estimates)
  
  # hack for size constraints
  if(!is.null(sizes.allowed)){
    length.p3 <- nrow(z.phase3)
    print("new length of phase 3")
    print(length.p3)
  }
  
  # mean statistics
  finalmean <- rep(0,num.effects)
  for(e in 1:num.effects) {
    finalmean[e] <- Reduce('+',z.phase3[,e]) / length.p3
  }
   
  # standard deviations
  finalsd <- rep(0, num.effects)
  for(e in 1:num.effects) {
    statse <- rep(0,length.p3)
    for(i in 1:length.p3) {
      statse[i] <- z.phase3[i,e]
    }
    finalsd[e] <- sd(statse)
  }
  finalse <- finalsd / sqrt(length.p3)
  
  # convergence ratios
  finalconvratios <- (finalmean - z.obs) / finalsd

  print("Estimated statistics after phase 3")
  print(finalmean)
  print("Estimates after phase 3")
  print(estimates.phase2)
  
  return(list("means" = finalmean, 
              "standard.deviations" = finalsd, 
              "standard.errors" = finalse, 
              "convergence.ratios" = finalconvratios,
              "inv.zcov" = inverted_matrices$inv.zcov,
              "inv.scaling" = inverted_matrices$inv.scaling))
  
}
