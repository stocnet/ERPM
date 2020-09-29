######################################################################
## Simulation and estimation of Exponential Random Partition Models ## 
## Functions used to generate partitions bases on permuting actors  ##
## (that keeps the same number of groups)                           ##
## Author: Marion Hoffman                                           ##
######################################################################

# permutations witht he same number of groups but sized are not constrained
generate_permutations <- function(partition, n_sample, thining){
  
  p <- partition[!is.na(partition)]
  num.nodes <- length(p)
  all.groups <- 1:max(p)
  all.actors <- 1:num.nodes
  
  all_permutations <- c()
  
  cpt <- 0
  cpt_thining <- 0
  all_generated <- F
  
  temp_p <- p
  
  while(!all_generated){
    
    cpt <- cpt + 1
    cpt_thining <- cpt_thining + 1
    
    rand_actor <- sample(all.actors,1)
    old_g <- temp_p[rand_actor]
    new_g <- sample(all.groups[-old_g],1)
    temp_p[rand_actor] <- new_g
    
    if(cpt_thining == thining) {
      p_toadd <- partition
      p_toadd[!is.na(partition)] <- temp_p
      all_permutations <- rbind(all_permutations, p_toadd)
      cpt_thining <- 0
    }
    
    if(length(all_permutations) > 0 && nrow(all_permutations) == n_sample){
      all_generated <- T
    }
  }
  
  return(all_permutations)
}


# permutations witht he same number of groups and sizes ARE constrained
generate_strict_permutations <- function(partition, n_sample, thining){
  
  p <- partition[!is.na(partition)]
  num.nodes <- length(p)
  all.groups <- 1:max(p)
  all.actors <- 1:num.nodes
  
  all_permutations <- c()
  
  cpt <- 0
  cpt_thining <- 0
  all_generated <- F
  
  temp_p <- p
  
  while(!all_generated){
    
    cpt <- cpt + 1
    cpt_thining <- cpt_thining + 1
    
    rand_actor_1 <- sample(all.actors,1)
    old_g_1 <- temp_p[rand_actor_1]
    members <- which(temp_p == old_g_1)
    rand_actor_2 <- sample(all.actors[-members],1)
    old_g_2 <- temp_p[rand_actor_2]
    temp_p[rand_actor_1] <- old_g_2
    temp_p[rand_actor_2] <- old_g_1
    
    if(cpt_thining == thining) {
      p_toadd <- partition
      p_toadd[!is.na(partition)] <- temp_p
      all_permutations <- rbind(all_permutations, p_toadd)
      cpt_thining <- 0
    }
    
    if(length(all_permutations) > 0 && nrow(all_permutations) == n_sample){
      all_generated <- T
    }
  }
  
  return(all_permutations)
}


# permutations witht he same number of groups and sizes ARE constrained with multiple permutations
generate_strict_permutations_multiple <- function(partitions, n_sample, thining){
  
  num.nodes <- nrow(partitions)
  num.obs <- ncol(partitions)
  
  all_permutations <- list()
  
  temp_p <- partitions
  
  for(i in 1:n_sample){

    for(o in 1:num.obs){
      temp_p[,o] <- generate_strict_permutations(temp_p[,o], 1, thining)
    }
    all_permutations[[i]] <- temp_p
  }
  
  return(all_permutations)
  
}
