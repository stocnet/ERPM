######################################################################
## Simulation and estimation of Exponential Random Partition Models ## 
## Functions used to generate partitions bases on permuting actors  ##
## (that keeps the same number of groups)                           ##
## Author: Marion Hoffman                                           ##
######################################################################

# permutations with the same number of groups but sized are not constrained
generate_permutations <- function(partition, n_sample, thining){
  
  p <- partition[!is.na(partition)]
  num.nodes <- length(p)
  all.groups <- 1:max(p)
  all.actors <- 1:num.nodes
  
  all_permutations <- list()
  
  cpt <- 0
  cpt_thining <- 0
  cpt_store <- 1
  all_generated <- F
  
  temp_p <- p
  
  while(!all_generated){
    
    cpt <- cpt + 1
    cpt_thining <- cpt_thining + 1
    
    found_actor <- F
    while(!found_actor) {
      rand_actor <- sample(all.actors,1)
      old_g <- temp_p[rand_actor]
      if(sum(temp_p == old_g) > 1) found_actor <- T
    }
    new_g <- sample(all.groups[-old_g],1)
    temp_p[rand_actor] <- new_g
    
    if(cpt_thining == thining) {
      p_toadd <- partition
      p_toadd[!is.na(partition)] <- temp_p
      all_permutations[[cpt_store]] <- p_toadd
      cpt_store <- cpt_store + 1
      cpt_thining <- 0
    }
    
    if(cpt_store > n_sample){
      all_generated <- T
    }
  }
  
  return(all_permutations)
}


# permutations with the same number of groups and sizes ARE constrained
generate_strict_permutations <- function(partition, n_sample, thining){
  
  p <- partition[!is.na(partition)]
  num.nodes <- length(p)
  all.groups <- 1:max(p)
  all.actors <- 1:num.nodes
  
  all_permutations <- list()
  
  cpt <- 0
  cpt_thining <- 0
  cpt_store <- 1
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
      all_permutations[[cpt_store]] <- p_toadd
      cpt_thining <- 0
      cpt_store <- cpt_store + 1
    }
    
    if(cpt_store > n_sample){
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
      temp_p[,o] <- generate_strict_permutations(temp_p[,o], 1, thining)[[1]]
    }
    all_permutations[[i]] <- temp_p
  }
  
  return(all_permutations)
  
}

# generate permutations of multiple partitions while keeping the permutations of acotrs constant over all observations
generate_strict_permutations_multiple2 <- function(partitions, n_sample){
  
  num.nodes <- nrow(partitions)
  num.obs <- ncol(partitions)
  
  all_permutations <- list()
  
  temp_p <- partitions
  
  for(i in 1:n_sample){
    
    newactors <- sample(1:num.nodes)
    
    for(o in 1:num.obs){
      temp_p[,o] <- partitions[newactors,o]
    }
    
    all_permutations[[i]] <- temp_p
  }
  
  return(all_permutations)
  
}



# permutations with the same number of groups and sizes ARE constrained, 
# and the number of groups of a certain type is also fixed
generate_strict_permutations_withtypes <- function(partition, types, n_sample, thining){
  
  p <- partition[!is.na(partition)]
  num.nodes <- length(p)
  all.groups <- 1:max(p)
  all.actors <- 1:num.nodes
  
  num.types <- max(types)
  all.types <- 1:max(types)
  
  all_permutations <- list()
  
  cpt <- 0
  cpt_thining <- 0
  cpt_store <- 1
  all_generated <- F
  
  temp_p <- p
  temp_types <- types
  
  while(!all_generated){
    
    cpt <- cpt + 1
    cpt_thining <- cpt_thining + 1
    
    if(cpt%%2 == 0){
      
      rand_actor_1 <- sample(all.actors,1)
      old_g_1 <- temp_p[rand_actor_1]
      members <- which(temp_p == old_g_1)
      rand_actor_2 <- sample(all.actors[-members],1)
      old_g_2 <- temp_p[rand_actor_2]
      temp_p[rand_actor_1] <- old_g_2
      temp_p[rand_actor_2] <- old_g_1
      
    } else {
      
      rand_group_1 <- sample(all.groups,1)
      old_t_1 <- temp_types[rand_group_1]
      sames <- which(temp_types == old_t_1)
      rand_group_2 <- sample(all.groups[-sames],1)
      old_t_2 <- temp_types[rand_group_2]
      temp_types[rand_group_1] <- old_t_2
      temp_types[rand_group_2] <- old_t_1
      
    }
    
    
    if(cpt_thining == thining) {
      p_toadd <- partition
      p_toadd[!is.na(partition)] <- temp_p
      all_permutations[[cpt_store]] <- list()
      all_permutations[[cpt_store]]$partition <- p_toadd
      all_permutations[[cpt_store]]$types <- temp_types
      cpt_thining <- 0
      cpt_store <- cpt_store + 1
    }
    
    if(cpt_store > n_sample){
      all_generated <- T
    }
  }
  
  return(all_permutations)
}

