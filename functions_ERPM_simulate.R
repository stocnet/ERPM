simulate_ERPM <- function(theta,nodes,effects,objects,partition,burnin,thining,num.steps,mini.steps,neighborhood,sizes.allowed,sizes.simulated) {
    
   
    num.nodes <- nrow(nodes)
    num.effects <- length(effects$names)

    allpartitions <- c()
    
    # simulate a large sample with the estimates found in phase 2  
    first.partition <- partition
    
    # instantiate with the starting network
    current.partition <- first.partition
    current.z <- computeStatistics(current.partition, nodes, effects, objects)
    current.logit <- theta * current.z
    
    # store the statistics collected for all networks simulated after the burn in
    all.z <-c()
    
    end.walk <- FALSE
    cpt <- 0
    cpt2 <- 0
    cpt_thining <- 0
    
    while(!end.walk){
      
      #print("partition")
      #print(current.partition)
      
      ## IF NEIGHBORHOOD IS: MERGES OR SPLITS OF 2 (Pi_1)
      if(neighborhood == 1 && is.null(sizes.allowed)) {
        new.partition <- sample_new_partition_p1(current.partition, mini.steps)
      } else if(neighborhood == 2 && is.null(sizes.allowed)){
        new.partition <- sample_new_partition_p2(current.partition, mini.steps)
      }
      if(neighborhood == 1 && !is.null(sizes.allowed)) {
        new.partition <- sample_new_partition_p1_restricted(current.partition, mini.steps, sizes.simulated)
      } else if(neighborhood == 2 && !is.null(sizes.allowed)){
        new.partition <- sample_new_partition_p2_restricted(current.partition, mini.steps, sizes.simulated)
      }
      
      #print("new")
      #print(new.partition)
      
      # compute new statistics only if it changed
      if(!all(current.partition == new.partition)) {
        new.z <- computeStatistics(new.partition, nodes, effects, objects)
        #new.z <- computeChangeStatistics(current.z, current.partition, new.partition, old_g1, old_g2, new_g1, new_g2, nodes, effects, objects)
      } else {
        new.z <- current.z
      }
      new.logit <- theta * new.z
      
      #print("stats")
      #print(new.z)
      
      # chose whether to change or not
      if(mini.steps == "normalized") {
        
        if(neighborhood == 1 && is.null(sizes.allowed)) {
          current.size <- compute_size_neighborhood_p1(current.partition)
          new.size <- compute_size_neighborhood_p1(new.partition)
          neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
        } else if(neighborhood == 2 && is.null(sizes.allowed)) {
          current.size <- compute_size_neighborhood_p2(current.partition)
          new.size <- compute_size_neighborhood_p2(new.partition)
          neighborhoods.ratio <- (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
        }
        if(neighborhood == 1 && !is.null(sizes.allowed)) {
          current.size <- compute_size_neighborhood_p1_restricted(current.partition, sizes.simulated)
          new.size <- compute_size_neighborhood_p1_restricted(new.partition, sizes.simulated)
          neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
        } else if(neighborhood == 2 && !is.null(sizes.allowed)) {
          current.size <- compute_size_neighborhood_p2_restricted(current.partition, sizes.simulated)
          new.size <- compute_size_neighborhood_p2_restricted(new.partition, sizes.simulated)
          neighborhoods.ratio <- (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
        }
        hastings.ratio <- (exp(sum(new.logit) - sum(current.logit))) * neighborhoods.ratio
        
      } else if(mini.steps == "selfloops"){
        
        hastings.ratio <- exp(sum(new.logit) - sum(current.logit))
        
      }
      proba.change <- min(1,hastings.ratio)
      change.made <- (runif(1,0,1) <= proba.change)
      
      #print("hastings proba.change and change.made")
      #print(hastings.ratio)
      #print(proba.change)
      #print(change.made)
      
      old.partition <- current.partition
      if(change.made) {
        current.partition <- new.partition
        current.z <- new.z
        current.logit <- new.logit
      }
      
      cpt <- cpt + 1
      if(cpt > burnin) cpt_thining <- cpt_thining + 1
      
      
      # storing results if all sizes are allowed
      if(is.null(sizes.allowed)){
        # store the results if we are out of burnin
        if(cpt >= burnin && cpt_thining == thining) {
          all.z <- rbind(all.z,current.z)
          allpartitions <- rbind(allpartitions,current.partition)
          cpt_thining <- 0
        }
        
        # stop the walk if number of steps reached 
        end.walk <- (cpt >= (burnin+thining*num.steps))
      }
      
      # storing results if sizes are constrained
      if(!is.null(sizes.allowed)){
        
        #print(check_sizes(current.partition,sizes.allowed))
        #print(current.partition)
        
        # if the partition is sampled with the right sizes
        # we keep the swapping strategy
        if(check_sizes(current.partition,sizes.allowed)){
          #neighborhood <- 2
          
          # if wrong sizes, we continue searching, and if we are out of tolerated simulated partitions, we go back
          # when out of allowed partitions we keep the divide/merge strategy
        } else {
          #neighborhood <- 2
          if(!check_sizes(current.partition,sizes.simulated)){
            current.partition <- old.partition
            if(check_sizes(current.partition,sizes.allowed)) neighborhood <- 1
          } 
        }
        
        # store the results if we are out of burnin
        if(cpt >= burnin && cpt_thining == thining) {
          cpt_thining <- 0
          if(check_sizes(current.partition,sizes.allowed)){
            all.z <- rbind(all.z,current.z)
            allpartitions <- rbind(allpartitions,current.partition)
            cpt2 <- cpt2 + 1
          } 
        }
        
        # stop the walk if number of steps reached 
        end.walk <- (cpt2 >= num.steps)
      }
      
      
    }
    
   return(list("partitions" = allpartitions,
               "stats" = all.z))  
  
}


simulate_barplots <- function(thetamin, thetamax, nodes,effects,objects,burnin,thining,num.steps,mini.steps) {
  
  thetamin <- -1
  thetamax <- 1
  burnin <- 100
  thining <- 100
  num.steps <- 500
  effects <- list(names = c("num_groups","product_sizes"), objects = "partition")
  neighborhood <- 1
  mini.steps <- "normalized"
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)

  classcounts <- readRDS("classes_counts_10nodes.rds")
  allclasses <- classcounts$classes
  alldistances <- readRDS("randmatrix_10nodes.rds")
  nclass <- dim(allclasses)[1]
  
  allthetas <- seq(-10,10,0.25)
  allmeans <- rep(0, length(allthetas))
  allsizes <- matrix(0, nrow = length(allthetas), ncol = num.nodes)
  allclass.counts <- matrix(0, nrow = length(allthetas), ncol = nclass)
  
  allthetas <- cbind(allp1,allp2)  
  allmeans <- rep(0, dim(allthetas)[1])
  allsizes <- matrix(0, nrow = dim(allthetas)[1], ncol = num.nodes)
  allclass.counts <- matrix(0, nrow = dim(allthetas)[1], ncol = nclass)
  
  
  #for(t in 1:length(allthetas)) {
  for(t in 1:dim(allthetas)[1]) {
    
    #theta <- allthetas[t]
    theta <- allthetas[t,]
    print(theta)
    
    # simulate a large sample with the estimates found in phase 2  
    first.partition <- 1 + rbinom(num.nodes, as.integer(num.nodes/2), 0.5)
    first.partition <- order_groupids(first.partition)
    
    # instantiate with the starting network
    current.partition <- first.partition
    current.z <- computeStatistics(current.partition, nodes, effects, objects)
    current.logit <- theta * current.z
    
    # store the distribution of the groups
    all.ngs <- c()
    all.z <- c()
    class.counts <- rep(0,nclass)
    
    end.walk <- FALSE
    cpt <- 0
    cpt_thining <- 0
    
    while(!end.walk){
      
      ## IF NEIGHBORHOOD IS: MERGES OR SPLITS OF 2 (Pi_1)
      if(neighborhood == 1) {
        new.partition <- sample_new_partition_p1(current.partition, mini.steps)
      } else if(neighborhood == 2){
        new.partition <- sample_new_partition_p2(current.partition, mini.steps)
      }
      
      # compute new statistics only if it changed
      if(!all(current.partition == new.partition)) {
        new.z <- computeStatistics(new.partition, nodes, effects, objects)
        #new.z <- computeChangeStatistics(current.z, current.partition, new.partition, old_g1, old_g2, new_g1, new_g2, nodes, effects, objects)
      } else {
        new.z <- current.z
      }
      new.logit <- theta * new.z
      
      # chose whether to change or not
      if(mini.steps == "normalized") {
        
        if(neighborhood == 1) {
          current.size <- compute_size_neighborhood_p1(current.partition)
          new.size <- compute_size_neighborhood_p1(new.partition)
          neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
        } else if(neighborhood == 2) {
          current.size <- compute_size_neighborhood_p2(current.partition)
          new.size <- compute_size_neighborhood_p2(new.partition)
          neighborhoods.ratio <- (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
        }
        hastings.ratio <- (exp(sum(new.logit) - sum(current.logit))) * neighborhoods.ratio
        
      } else if(mini.steps == "selfloops"){
        
        hastings.ratio <- exp(sum(new.logit) - sum(current.logit))
        
      }
      proba.change <- min(1,hastings.ratio)
      change.made <- (runif(1,0,1) <= proba.change)
      
      if(change.made) {
        current.partition <- new.partition
        current.z <- new.z
        current.logit <- new.logit
      }
      
      cpt <- cpt + 1
      if(cpt > burnin) cpt_thining <- cpt_thining + 1
      
      # store the results of the currently simulated network if
      # we are out of burnin
      if(cpt >= burnin && cpt_thining == thining) {
        all.z <- c(all.z, current.z[2])
        #all.z <- c(all.z, current.z)
        ngs <- rep(0,num.nodes)
        for(g in 1:max(current.partition)) {
          ngs[sum(current.partition == g)] <- ngs[sum(current.partition == g)] + 1 
        }
        all.ngs <- rbind(all.ngs, ngs)
        for(c in 1:nclass){
          if(all(allclasses[c,] == ngs)) {
            class <- c
          }
        }
        class.counts[class] <-class.counts[class] + 1
        cpt_thining <- 0
      }
      
      # stop the walk if number of steps reached 
      end.walk <- (cpt >= (num.steps * thining) + burnin)
      
    }
    
    # compute the average statistics and the final network generated
    allmeans[t] <- mean(all.z) 
    allsizes[t,] <- colMeans(all.ngs) 
    allclass.counts[t,] <- class.counts
  }
  
  # plot outcome
  df <- data.frame(parameter=allthetas[,2],average_statistic=allmeans)
  ggplot(df, aes(x=parameter, y=average_statistic)) + geom_point()
  df <- data.frame(parameter=allthetas,average_statistic=allmeans)
  ggplot(df, aes(x=parameter, y=average_statistic)) + geom_point()
  plot(allthetas, allmeans)
  plot(allthetas[,2], allmeans)
  
  # plot all sizes
  df <- data.frame(parameter=allthetas, group_count=as.vector(allsizes), 
                   group=rep(paste0("size", 1:10), each=length(allthetas)))
  ggplot(data = df, aes(x=parameter, y=group_count)) + geom_line(aes(colour=group), size = 1.5)
  df <- data.frame(parameter=allthetas[,2], group_count=as.vector(allsizes), 
                   variable=rep(paste0("size", 1:10), each=dim(allthetas)[1]))
  ggplot(data = df, aes(x=parameter, y= group_count)) + geom_line(aes(colour=variable), size = 1.5)
  
  # plot sizes
  df2 <- data.frame(group_size=1:10, count=allsizes[67,])
  ggplot(data=df2, aes(x=group_size, y=count)) +
    geom_bar(stat="identity")
  
  # plot classes
  class.counts <- allclass.counts[24,]
  mode <- which(class.counts == max(class.counts))
  distances <- 1-alldistances[mode,]
  orderedindexes <- order(distances)
  labels <- c()
  for(c in 1:nclass){
    ngs <- allclasses[c,]
    cpt <- n
    name <- ""
    while(cpt>0){
      if(ngs[cpt]>0){
        name <- paste(name,cpt,sep=",")
        ngs[cpt] <- ngs[cpt] - 1
      } else {
        cpt <- cpt-1
      }
    }
    name <- substr(name,2,nchar(name))
    labels <- rbind(labels,name[])
  }
  ordlabels <- labels[orderedindexes]
  orddistances <- distances[orderedindexes]
  ordcounts <- class.counts[orderedindexes]
  df <- data.frame(class = ordlabels, proportion=ordcounts)
  df$class <- factor(df$class, levels = df$class)
  ggplot(data=df, aes(x=class, y=proportion)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  df <- data.frame(averagedistance = orddistances, proportion=ordcounts)
  ggplot(data=df, aes(x=averagedistance, y=proportion)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


simulate_barplots_homophily <- function(thetamin, thetamax, nodes,effects,objects,burnin,thining,num.steps,mini.steps) {
  
  nodes <- data.frame(label = paste("Actor",1:20),
                      att = c(rep(0,10),rep(1,10)))
  
  # Estimation to have p1, with n=20 and 8 groups
  
  partition <- c(1,1,1,2,2,2,3,3,3,3,4,4,4,4,5,5,6,6,7,7,8,8)
  effects <- list( names = c("num_ties"),
                   objects = c("partition"))
  objects <- list()
  multiplicationfactor <- 30
  gainfactor <- 0.1 
  mini.steps <- "normalized"
  burnin <- 100
  thining <- 100
  length.p1 <- 500
  min.iter.p2 <- 200
  max.iter.p2 <- 300
  num.steps.p2 <- 10
  length.p3 <- 500
  neighborhood <- 1

  results.erpm <- estimate_ERPM(partition, nodes, objects, effects, 
                                  startingestimates=0, multiplicationfactor, gainfactor, mini.steps, burnin, thining,
                                  length.p1, min.iter.p2, max.iter.p2, num.steps.p2, length.p3, neighborhood, fixed.estimates=NULL)
    
  p1 <- results.erpm$est[1]
  
  
  # Simulation to have the group sizes distribution
  burnin <- 100
  thining <- 100
  num.steps <- 2000
  effects <- list(names = c("num_groups"), objects = "partition")
  neighborhood <- 1
  mini.steps <- "normalized"
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)

  allthetas <- p1
  allsizes <- rep(0, num.nodes)
  
  for(t in 1:length(allthetas)) {
    
    theta <- allthetas[t]
    print(theta)
    
    # simulate a large sample with the estimates found in phase 2  
    first.partition <- 1 + rbinom(num.nodes, as.integer(num.nodes/2), 0.5)
    first.partition <- order_groupids(first.partition)
    
    # instantiate with the starting network
    current.partition <- first.partition
    current.z <- computeStatistics(current.partition, nodes, effects, objects)
    current.logit <- theta * current.z
    
    # store the distribution of the groups
    all.ngs <- c()
    all.z <- c()
    class.counts <- rep(0,nclass)
    
    end.walk <- FALSE
    cpt <- 0
    cpt_thining <- 0
    
    while(!end.walk){
      
      ## IF NEIGHBORHOOD IS: MERGES OR SPLITS OF 2 (Pi_1)
      if(neighborhood == 1) {
        new.partition <- sample_new_partition_p1(current.partition, mini.steps)
      } else if(neighborhood == 2){
        new.partition <- sample_new_partition_p2(current.partition, mini.steps)
      }
      
      # compute new statistics only if it changed
      if(!all(current.partition == new.partition)) {
        new.z <- computeStatistics(new.partition, nodes, effects, objects)
        #new.z <- computeChangeStatistics(current.z, current.partition, new.partition, old_g1, old_g2, new_g1, new_g2, nodes, effects, objects)
      } else {
        new.z <- current.z
      }
      new.logit <- theta * new.z
      
      # chose whether to change or not
      if(mini.steps == "normalized") {
        
        if(neighborhood == 1) {
          current.size <- compute_size_neighborhood_p1(current.partition)
          new.size <- compute_size_neighborhood_p1(new.partition)
          neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
        } else if(neighborhood == 2) {
          current.size <- compute_size_neighborhood_p2(current.partition)
          new.size <- compute_size_neighborhood_p2(new.partition)
          neighborhoods.ratio <- (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
        }
        hastings.ratio <- (exp(sum(new.logit) - sum(current.logit))) * neighborhoods.ratio
        
      } else if(mini.steps == "selfloops"){
        
        hastings.ratio <- exp(sum(new.logit) - sum(current.logit))
        
      }
      proba.change <- min(1,hastings.ratio)
      change.made <- (runif(1,0,1) <= proba.change)
      
      if(change.made) {
        current.partition <- new.partition
        current.z <- new.z
        current.logit <- new.logit
      }
      
      cpt <- cpt + 1
      if(cpt > burnin) cpt_thining <- cpt_thining + 1
      
      # store the results of the currently simulated network if
      # we are out of burnin
      if(cpt >= burnin && cpt_thining == thining) {
        all.z <- c(all.z, current.z)
        ngs <- rep(0,num.nodes)
        for(g in 1:max(current.partition)) {
          ngs[sum(current.partition == g)] <- ngs[sum(current.partition == g)] + 1 
        }
        all.ngs <- rbind(all.ngs, ngs)
        cpt_thining <- 0
      }
      
      # stop the walk if number of steps reached 
      end.walk <- (cpt >= (num.steps * thining) + burnin)
      
    }
    
    # compute the average statistics and the final network generated
    allsizes <- colMeans(all.ngs) 
  }
  
  df <- data.frame(group_size=1:20, count=allsizes)
  ggplot(data=df, aes(x=group_size, y=count)) +
    geom_bar(stat="identity")
  
  
  # simulations with the homophily parameter
  allp2 <- seq(-2,2,0.05)
  allp1 <- rep(p1,length(allp2))
  
  thetamin <- -1
  thetamax <- 1
  burnin <- 100
  thining <- 100
  num.steps <- 2000
  effects <- list(names = c("num_ties","same"), objects = c("partition","att"))
  neighborhood <- 1
  mini.steps <- "normalized"
  
  num.nodes <- nrow(nodes)
  num.effects <- length(effects$names)

  allthetas <- cbind(allp1,allp2)  
  allmeans <- rep(0, dim(allthetas)[1])
  allsizes <- matrix(0, nrow = dim(allthetas)[1], ncol = num.nodes)
  
  
  #for(t in 1:length(allthetas)) {
  for(t in 1:dim(allthetas)[1]) {
    
    #theta <- allthetas[t]
    theta <- allthetas[t,]
    print(theta)
    
    # simulate a large sample with the estimates found in phase 2  
    first.partition <- 1 + rbinom(num.nodes, as.integer(num.nodes/2), 0.5)
    first.partition <- order_groupids(first.partition)
    
    # instantiate with the starting network
    current.partition <- first.partition
    current.z <- computeStatistics(current.partition, nodes, effects, objects)
    current.logit <- theta * current.z
    
    # store the distribution of the groups
    all.ngs <- c()
    all.z <- c()
    class.counts <- rep(0,nclass)
    
    end.walk <- FALSE
    cpt <- 0
    cpt_thining <- 0
    
    while(!end.walk){
      
      ## IF NEIGHBORHOOD IS: MERGES OR SPLITS OF 2 (Pi_1)
      if(neighborhood == 1) {
        new.partition <- sample_new_partition_p1(current.partition, mini.steps)
      } else if(neighborhood == 2){
        new.partition <- sample_new_partition_p2(current.partition, mini.steps)
      }
      
      # compute new statistics only if it changed
      if(!all(current.partition == new.partition)) {
        new.z <- computeStatistics(new.partition, nodes, effects, objects)
        #new.z <- computeChangeStatistics(current.z, current.partition, new.partition, old_g1, old_g2, new_g1, new_g2, nodes, effects, objects)
      } else {
        new.z <- current.z
      }
      new.logit <- theta * new.z
      
      # chose whether to change or not
      if(mini.steps == "normalized") {
        
        if(neighborhood == 1) {
          current.size <- compute_size_neighborhood_p1(current.partition)
          new.size <- compute_size_neighborhood_p1(new.partition)
          neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
        } else if(neighborhood == 2) {
          current.size <- compute_size_neighborhood_p2(current.partition)
          new.size <- compute_size_neighborhood_p2(new.partition)
          neighborhoods.ratio <- (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
        }
        hastings.ratio <- (exp(sum(new.logit) - sum(current.logit))) * neighborhoods.ratio
        
      } else if(mini.steps == "selfloops"){
        
        hastings.ratio <- exp(sum(new.logit) - sum(current.logit))
        
      }
      proba.change <- min(1,hastings.ratio)
      change.made <- (runif(1,0,1) <= proba.change)
      
      if(change.made) {
        current.partition <- new.partition
        current.z <- new.z
        current.logit <- new.logit
      }
      
      cpt <- cpt + 1
      if(cpt > burnin) cpt_thining <- cpt_thining + 1
      
      # store the results of the currently simulated network if
      # we are out of burnin
      if(cpt >= burnin && cpt_thining == thining) {
        all.z <- c(all.z, current.z[2])
        #all.z <- c(all.z, current.z)
        ngs <- rep(0,num.nodes)
        for(g in 1:max(current.partition)) {
          ngs[sum(current.partition == g)] <- ngs[sum(current.partition == g)] + 1 
        }
        all.ngs <- rbind(all.ngs, ngs)
        cpt_thining <- 0
      }
      
      # stop the walk if number of steps reached 
      end.walk <- (cpt >= (num.steps * thining) + burnin)
      
    }
    
    # compute the average statistics and the final network generated
    allmeans[t] <- mean(all.z) 
    allsizes[t,] <- colMeans(all.ngs) 
  }
  
  # plot outcome
  df <- data.frame(parameter=allthetas[,2],average_statistic=allmeans)
  ggplot(df, aes(x=parameter, y=average_statistic)) + geom_point()
  df <- data.frame(parameter=allthetas,average_statistic=allmeans)
  ggplot(df, aes(x=parameter, y=average_statistic)) + geom_point()
  plot(allthetas, allmeans)
  plot(allthetas[,2], allmeans)
  
  # plot all sizes
  df <- data.frame(parameter=allthetas, group_count=as.vector(allsizes), 
                   group=rep(paste0("size", 1:10), each=length(allthetas)))
  ggplot(data = df, aes(x=parameter, y=group_count)) + geom_line(aes(colour=group), size = 1.5)
  df <- data.frame(parameter=allthetas[,2], group_count=as.vector(allsizes), 
                   variable=rep(paste0("size", 1:10), each=dim(allthetas)[1]))
  ggplot(data = df, aes(x=parameter, y= group_count)) + geom_line(aes(colour=variable), size = 1.5)
  
  # plot sizes
  df2 <- data.frame(group_size=1:10, count=allsizes[81,])
  ggplot(data=df2, aes(x=group_size, y=count)) +
    geom_bar(stat="identity")
  
  # plot classes
  class.counts <- allclass.counts[24,]
  mode <- which(class.counts == max(class.counts))
  distances <- 1-alldistances[mode,]
  orderedindexes <- order(distances)
  labels <- c()
  for(c in 1:nclass){
    ngs <- allclasses[c,]
    cpt <- n
    name <- ""
    while(cpt>0){
      if(ngs[cpt]>0){
        name <- paste(name,cpt,sep=",")
        ngs[cpt] <- ngs[cpt] - 1
      } else {
        cpt <- cpt-1
      }
    }
    name <- substr(name,2,nchar(name))
    labels <- rbind(labels,name[])
  }
  ordlabels <- labels[orderedindexes]
  orddistances <- distances[orderedindexes]
  ordcounts <- class.counts[orderedindexes]
  df <- data.frame(class = ordlabels, proportion=ordcounts)
  df$class <- factor(df$class, levels = df$class)
  ggplot(data=df, aes(x=class, y=proportion)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  df <- data.frame(averagedistance = orddistances, proportion=ordcounts)
  ggplot(data=df, aes(x=averagedistance, y=proportion)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


find_fixedparameter <- function(est, nodes,effects,objects){
  
  nodes <- data.frame(label = paste("Actor",1:10))
  
  allp2 <- seq(-1,1,0.1)
  allp1 <- rep(0,length(allp2))
  
  partition <- c(1,1,1,2,2,2,3,3,3,3)
  effects <- list( names = c("num_ties","sum_log_factorials"),
                   objects = c("partition","partition"))
  objects <- list()
  multiplicationfactor <- 30
  gainfactor <- 0.1 
  mini.steps <- "normalized"
  burnin <- 200
  thining <- 100
  length.p1 <- 200
  min.iter.p2 <- 200
  max.iter.p2 <- 300
  num.steps.p2 <- 6
  length.p3 <- 300
  neighborhood <- 1
  
  for(p in 1:length(allp2)) {
    
    print(allp2[p])
    startingestimates <- c(0,allp2[p])
    fixed.estimates <- list()
    fixed.estimates[[1]] <- NULL
    fixed.estimates[[2]] <- allp2[p]
    
    results.erpm <- estimate_ERPM(partition, nodes, objects, effects, 
                                  startingestimates, multiplicationfactor, gainfactor, mini.steps, burnin, thining,
                                  length.p1, min.iter.p2, max.iter.p2, num.steps.p2, length.p3, neighborhood, fixed.estimates,
                                  sizes.allowed = NULL,sizes.simulated = NULL,double.averaging = F)
    
    allp1[p] <- results.erpm$est[1]
                                          
  }

  lo <- predict(loess(allp1~allp2),allp2)
  plot(allp2,allp1)
  lines(allp2,predict(loess(allp1~allp2),allp2), col="red", lwd=2)
  allp1_smoothed <- lo
  
  # store the stats and distribution of the groups
  all.sums <- rep(0,length(allp2))
  all.sizes <- matrix(0,nrow=length(allp2),ncol=10)
  num.steps <- 1000
  thining <- 200
  
  for(p in 1:length(allp2)) {
    
    theta <- c(allp1[p],allp2[p])
    
    # simulate a large sample with the estimates found in phase 2  
    first.partition <- 1 + rbinom(num.nodes, as.integer(num.nodes/2), 0.5)
    first.partition <- order_groupids(first.partition)
    
    # instantiate with the starting network
    current.partition <- first.partition
    current.z <- computeStatistics(current.partition, nodes, effects, objects)
    current.logit <- theta * current.z
    
    # store the distribution of the groups
    all.ngs <- c()
    all.z <- c()
    
    end.walk <- FALSE
    cpt <- 0
    cpt_thining <- 0
    
    while(!end.walk){
      
      ## IF NEIGHBORHOOD IS: MERGES OR SPLITS OF 2 (Pi_1)
      if(neighborhood == 1) {
        new.partition <- sample_new_partition_p1(current.partition, mini.steps)
      } else if(neighborhood == 2){
        new.partition <- sample_new_partition_p2(current.partition, mini.steps)
      }
      
      # compute new statistics only if it changed
      if(!all(current.partition == new.partition)) {
        new.z <- computeStatistics(new.partition, nodes, effects, objects)
        #new.z <- computeChangeStatistics(current.z, current.partition, new.partition, old_g1, old_g2, new_g1, new_g2, nodes, effects, objects)
      } else {
        new.z <- current.z
      }
      new.logit <- theta * new.z
      
      # chose whether to change or not
      if(mini.steps == "normalized") {
        
        if(neighborhood == 1) {
          current.size <- compute_size_neighborhood_p1(current.partition)
          new.size <- compute_size_neighborhood_p1(new.partition)
          neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
        } else if(neighborhood == 2) {
          current.size <- compute_size_neighborhood_p2(current.partition)
          new.size <- compute_size_neighborhood_p2(new.partition)
          neighborhoods.ratio <- (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
        }
        hastings.ratio <- (exp(sum(new.logit) - sum(current.logit))) * neighborhoods.ratio
        
      } else if(mini.steps == "selfloops"){
        
        hastings.ratio <- exp(sum(new.logit) - sum(current.logit))
        
      }
      proba.change <- min(1,hastings.ratio)
      change.made <- (runif(1,0,1) <= proba.change)
      
      if(change.made) {
        current.partition <- new.partition
        current.z <- new.z
        current.logit <- new.logit
      }
      
      cpt <- cpt + 1
      if(cpt > burnin) cpt_thining <- cpt_thining + 1
      
      # store the results of the currently simulated network if
      # we are out of burnin
      if(cpt >= burnin && cpt_thining == thining) {
        ngs <- rep(0,num.nodes)
        for(g in 1:max(current.partition)) {
          ngs[sum(current.partition == g)] <- ngs[sum(current.partition == g)] + 1 
        }
        all.ngs <- rbind(all.ngs, ngs)
        all.z <- rbind(all.z,current.z[2])
        cpt_thining <- 0
      }
      
      # stop the walk if number of steps reached 
      end.walk <- (cpt >= (num.steps * thining) + burnin)
      
    }
    
    # compute the average statistics and the final network generated
    all.sizes[p,] <- colMeans(all.ngs) 
    all.sums[p] <- mean(all.z)
  }
  
  df <- data.frame(x=firstallp$p2, val=as.vector(firstallsizes), 
                   variable=rep(paste0("size", 1:10), each=length(firstallp$p2)))
  ggplot(data = df, aes(x=x, y=val)) + geom_line(aes(colour=variable), size = 1.5)
  
  df2 <- data.frame(group_size=1:10, proportion=firstallsizes[26,])
  ggplot(data=df2, aes(x=group_size, y=proportion)) +
    geom_bar(stat="identity")
  
  df3 <- data.frame(sim=1:1000, z=all.z)
  ggplot(data=df3, aes(x=z, y=sim)) +
    geom_bar(stat="identity")
  
}


draw_distribution <- function() {
  
  allp2 <- seq(0.8,5,0.1)
  allp1 <- rep(0,length(allp2))
  
  partition <- c(1,1,1,2,2,2,3,3,3,3)
  effects <- list( names = c("num_groups","sizes_squared"),
                   objects = c("partition","partition"))
  objects <- list()
  multiplicationfactor <- 30
  gainfactor <- 0.1 
  mini.steps <- "normalized"
  burnin <- 100
  thining <- 100
  length.p1 <- 500
  min.iter.p2 <- 100
  max.iter.p2 <- 200
  num.steps.p2 <- 10
  length.p3 <- 500
  neighborhood <- 1
  
  mine <- minestimates[e]
  maxe <- maxestimates[e]
  allparameters <- seq(mine, maxe, by = 0.5)
  allmeans <- rep(0,length(allparameters))
  
  effects2 <- list(names = effects$names[e], objects = effects$objects[e])
  
  for(p in 1:length(allparameters)){
    estimates2 <- allparameters[p]
    chain <- draw_Metropolis(estimates2, first.partition, nodes, effects, objects, burnin, thining, num.steps, mini.steps, neighborhood)
    draws <- chain$draws
    allmeans[p] <- mean(draws[,1])
    
    print(paste("parameter",estimates2))
    print(paste("mean statistic", mean(draws[,1])))
  }
  
  plots[[e]] <- list(x = allparameters,
                     y = allmeans)
  
}


simulate_classes <- function(){
  
  thetas <- seq(-5,5,0.1)
  burnin <- 100
  thining <- 200
  num.steps <- 1000
  
  num.nodes <- 10
  num.effects <- 1
  effects2 <- list(names="num_groups",objects="partition")
  
  # simulate a large sample with the estimates found in phase 2  
  first.partition <- 1 + rbinom(num.nodes, as.integer(num.nodes/2), 0.5)
  first.partition <- order_groupids(first.partition)
  
  # instantiate with the starting network
  current.partition <- first.partition
  current.z <- computeStatistics(current.partition, nodes, effects2, objects)
  current.logit <- theta * current.z
  
  # store the distribution of the classes
  allclasses <- classescountsn10$classes
  alldistances <- randmatn10
  nclasses <- dim(allclasses)[1]
  class.counts <- rep(0,nclasses)
  
  end.walk <- FALSE
  cpt <- 0
  cpt_thining <- 0
  
  while(!end.walk){
    
    ## IF NEIGHBORHOOD IS: MERGES OR SPLITS OF 2 (Pi_1)
    if(neighborhood == 1) {
      new.partition <- sample_new_partition_p1(current.partition, mini.steps)
    } else if(neighborhood == 2){
      new.partition <- sample_new_partition_p2(current.partition, mini.steps)
    }
    
    # compute new statistics only if it changed
    if(!all(current.partition == new.partition)) {
      new.z <- computeStatistics(new.partition, nodes, effects2, objects)
      #new.z <- computeChangeStatistics(current.z, current.partition, new.partition, old_g1, old_g2, new_g1, new_g2, nodes, effects, objects)
    } else {
      new.z <- current.z
    }
    new.logit <- theta * new.z
    
    # chose whether to change or not
    if(mini.steps == "normalized") {
      
      if(neighborhood == 1) {
        current.size <- compute_size_neighborhood_p1(current.partition)
        new.size <- compute_size_neighborhood_p1(new.partition)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else if(neighborhood == 2) {
        current.size <- compute_size_neighborhood_p2(current.partition)
        new.size <- compute_size_neighborhood_p2(new.partition)
        neighborhoods.ratio <- (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
      }
      hastings.ratio <- (exp(sum(new.logit) - sum(current.logit))) * neighborhoods.ratio
      
    } else if(mini.steps == "selfloops"){
      
      hastings.ratio <- exp(sum(new.logit) - sum(current.logit))
      
    }
    proba.change <- min(1,hastings.ratio)
    change.made <- (runif(1,0,1) <= proba.change)
    
    if(change.made) {
      current.partition <- new.partition
      current.z <- new.z
      current.logit <- new.logit
    }
    
    cpt <- cpt + 1
    if(cpt > burnin) cpt_thining <- cpt_thining + 1
    
    # store the results of the currently simulated network if
    # we are out of burnin
    if(cpt >= burnin && cpt_thining == thining) {
      ngs <- rep(0,num.nodes)
      for(g in 1:max(current.partition)) {
        ngs[sum(current.partition == g)] <- ngs[sum(current.partition == g)] + 1 
      }
      for(c in 1:dim(allclasses)[1]){
        if(all(allclasses[c,] == ngs)) {
          class <- c
        }
      }
      class.counts[class] <-class.counts[class] + 1
      cpt_thining <- 0
    }
    
    # stop the walk if number of steps reached 
    end.walk <- (cpt >= (num.steps * thining) + burnin)
    
  }
  
  # plot to check bimodality
  df2 <- data.frame(group_size=1:nclasses, proportion=class.counts)
  ggplot(data=df2, aes(x=group_size, y=proportion)) +
    geom_bar(stat="identity")
  
  # create labels
  labels <- c()
  for(c in 1:nclasses){
    ngs <- allclasses[c,]
    cpt <- n
    name <- ""
    while(cpt>0){
      if(ngs[cpt]>0){
        name <- paste(name,cpt,sep=",")
        ngs[cpt] <- ngs[cpt] - 1
      } else {
        cpt <- cpt-1
      }
    }
    name <- substr(name,2,nchar(name))
    labels <- rbind(labels,name[])
  }
  
  # plot around the mode
  mode <- 1
  mode <- which(class.counts == max(class.counts))
  distances <- 1-alldistances[mode,]
  orderedindexes <- order(distances)
  ordlabels <- labels[orderedindexes]
  orddistances <- distances[orderedindexes]
  ordcounts <- class.counts[orderedindexes]
  df <- data.frame(class = ordlabels, proportion=ordcounts)
  ggplot(data=df, aes(x=class, y=proportion)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  df <- data.frame(averagedistance = orddistances, proportion=ordcounts)
  ggplot(data=df, aes(x=averagedistance, y=proportion)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
}


plot_parameter_map <- function(){
  
  thetas1 <- seq(-7,7,0.5)
  thetas2 <- seq(-10,10,0.5)
  burnin <- 100
  thining <- 200
  num.steps <- 1000
  
  num.nodes <- 10
  num.effects <- 1
  effects <- list(names=c("sizes_squared","isolates"),
                  objects=c("partition","partition"))
  
  # simulate a large sample with the estimates found in phase 2  
  first.partition <- 1 + rbinom(num.nodes, as.integer(num.nodes/2), 0.5)
  first.partition <- order_groupids(first.partition)
  
  # instantiate with the starting network
  current.partition <- first.partition
  current.z <- computeStatistics(current.partition, nodes, effects, objects)
  current.logit <- theta * current.z
  
  # distribution of the classes
  allclasses <- classescountsn10$classes
  alldistances <- randmatn10
  nclasses <- dim(allclasses)[1]
  
  # store mode classes
  allcounts <- list()
  cpt <- 1
 
  for(t1 in 1:length(thetas1)) {
    
    for(t2 in 1:length(thetas2)){
      
      theta <- c(thetas1[t1],thetas2[t2])
      class.counts <- rep(0,nclasses)
      print(theta)
      
      end.walk <- FALSE
      cpt <- 0
      cpt_thining <- 0
      
      while(!end.walk){
        
        ## IF NEIGHBORHOOD IS: MERGES OR SPLITS OF 2 (Pi_1)
        if(neighborhood == 1) {
          new.partition <- sample_new_partition_p1(current.partition, mini.steps)
        } else if(neighborhood == 2){
          new.partition <- sample_new_partition_p2(current.partition, mini.steps)
        }
        
        # compute new statistics only if it changed
        if(!all(current.partition == new.partition)) {
          new.z <- computeStatistics(new.partition, nodes, effects, objects)
          #new.z <- computeChangeStatistics(current.z, current.partition, new.partition, old_g1, old_g2, new_g1, new_g2, nodes, effects, objects)
        } else {
          new.z <- current.z
        }
        new.logit <- theta * new.z
        
        # chose whether to change or not
        if(mini.steps == "normalized") {
          
          if(neighborhood == 1) {
            current.size <- compute_size_neighborhood_p1(current.partition)
            new.size <- compute_size_neighborhood_p1(new.partition)
            neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
          } else if(neighborhood == 2) {
            current.size <- compute_size_neighborhood_p2(current.partition)
            new.size <- compute_size_neighborhood_p2(new.partition)
            neighborhoods.ratio <- (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
          }
          hastings.ratio <- (exp(sum(new.logit) - sum(current.logit))) * neighborhoods.ratio
          
        } else if(mini.steps == "selfloops"){
          
          hastings.ratio <- exp(sum(new.logit) - sum(current.logit))
          
        }
        proba.change <- min(1,hastings.ratio)
        change.made <- (runif(1,0,1) <= proba.change)
        
        if(change.made) {
          current.partition <- new.partition
          current.z <- new.z
          current.logit <- new.logit
        }
        
        cpt <- cpt + 1
        if(cpt > burnin) cpt_thining <- cpt_thining + 1
        
        # store the results of the currently simulated network if
        # we are out of burnin
        if(cpt >= burnin && cpt_thining == thining) {
          ngs <- rep(0,num.nodes)
          for(g in 1:max(current.partition)) {
            ngs[sum(current.partition == g)] <- ngs[sum(current.partition == g)] + 1 
          }
          for(c in 1:dim(allclasses)[1]){
            if(all(allclasses[c,] == ngs)) {
              class <- c
            }
          }
          class.counts[class] <-class.counts[class] + 1
          cpt_thining <- 0
        }
        
        # stop the walk if number of steps reached 
        end.walk <- (cpt >= (num.steps * thining) + burnin)
        
      }
      
      # find modes: more than 10% presence
      allcounts[[cpt]] <- list(parameters = theta,
                              counts = class.counts)
      cpt <- cpt + 1
    }
    
  }
  
}


plot_outcome_space <- function(){
  
  classcounts_n10 <- readRDS("classes_counts_10nodes.rds")
  allclasses <- classcounts_n10$classes
  
  nclass <- dim(allclasses)[1]
  n <- dim(allclasses)[2]
  
  # store mean statistics for s1 and s2 for each class
  s0 <- rep(0,nclass)
  s1 <- rep(0,nclass)
  s2 <- rep(0,nclass)
  s3 <- rep(0,nclass)
  s4 <- rep(0,nclass)
  s5 <- rep(0,nclass)
  s6 <- rep(0,nclass)
  s7 <- rep(0,nclass)
  s8 <- rep(0,nclass)
  s9 <- rep(0,nclass)
  s10 <- rep(0,nclass)
  s11 <- rep(0,nclass)
  
  for(c in 1:nclass){
    
    ngs <- allclasses[c,]
    
    # s0: number of groups
    for(g in 1:n){
      s0[c] <- s0[c] + ngs[g]
    }
    
    # s1: sum of squared sizes
    for(g in 1:n){
      s1[c] <- s1[c] + ngs[g]*g^2
    }
    
    # s2: sum of squared sizes normalized
    for(g in 1:n){
      s2[c] <- s2[c] + ngs[g]*g^2
    }
    s2[c] <- s2[c] / sum(ngs)
    
    # s3: number of isolates
    s3[c] <- ngs[1]
    
    # s4
    for(g in 1:n){
      s4[c] <- s4[c] + ngs[g]*g*log(g^2)
    }
    
    # s5 range
    s5[c] <- (max(which(ngs>0)) - min(which(ngs>0)) ) 
    
    # s6 num ties
    for(g in 1:n){
      s6[c] <- s6[c] + ngs[g]*g*(g-1)/2
    }
    
    # s7 num triangles
    for(g in 1:n){
      if(g == 3) s7[c] <- s7[c] + 1
      if(g > 3) s7[c] <- s7[c] + ngs[g]*dim(combn(g,3))[2]
    }
    
    # s8 num groups>2
    for(g in 1:n){
      if(g > 2) s8[c] <- s8[c] + ngs[g]
    }
    
    # s9 sum of factorials
    s9[c] <- ngs[1]
    for(g in 2:n){
      s9[c] <- s9[c] + ngs[g]*factorial(g-1)
    }
    
    # s10 log product sizes
    s10[c] <- 1
    for(g in 1:n){
      s10[c] <- s10[c]*g^(ngs[g])
    }
    s10[c] <- log(s10[c])
    
    # s11 sum log factorials
    s11[c] <- 0
    for(g in 3:n){
      s11[c] <- s11[c] + log(factorial(g-1))*(ngs[g])
    }
    
  }
  
  # plot outcome space
  df <- data.frame(number_groups=s0,squared_sizes=s1,squared_sizes_norm=s2)
  ggplot(df) + 
      geom_point( aes(x=number_groups, y=squared_sizes), color = "black") + 
      geom_point( aes(x=number_groups, y=squared_sizes_norm), color = "red") 
  
  # plot dispersion of classes
  mode <- 1
  distances <- 1-alldistances[mode,]
  orderedindexes <- order(distances)
  labels <- c()
  for(c in 1:nclass){
    ngs <- allclasses[c,]
    cpt <- n
    name <- ""
    while(cpt>0){
      if(ngs[cpt]>0){
        name <- paste(name,cpt,sep=",")
        ngs[cpt] <- ngs[cpt] - 1
      } else {
        cpt <- cpt-1
      }
    }
    name <- substr(name,2,nchar(name))
    labels <- rbind(labels,name[])
  }
  ordlabels <- labels[orderedindexes]
  orddistances <- distances[orderedindexes]
  ordstats <- s24[orderedindexes]
  df <- data.frame(class = ordlabels, s=ordstats)
  df$class <- factor(df$class, levels = df$class)
  ggplot(data=df, aes(x=class, y=s)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


estimate60nodes <- function(){
  nodes <- data.frame(label = paste("Actor",1:60),
                      att = sample(c(0,1), replace=TRUE, size=60))
  num.nodes <- 60
  partition <- c(rep(1,2),
                 rep(2,3),
                 rep(3,4),rep(4,4),rep(5,4),rep(6,4),rep(7,4),
                 rep(8,5),rep(9,5),rep(10,5),rep(11,5),rep(12,5),rep(13,5),rep(14,5))
  effects <- list( names = c("num_groups"),
                   objects = c("partition"))
  effects <- list( names = c("num_groups","sizes_squared"),
                   objects = c("partition","partition"))
  objects <- list()
  multiplicationfactor <- 30
  gainfactor <- 0.1 
  mini.steps <- "normalized"
  burnin <- 200
  thining <- 50
  length.p1 <- 1000
  min.iter.p2 <- 1000
  max.iter.p2 <- 1200
  num.steps.p2 <- 2
  length.p3 <- 1500
  neighborhood <- 2
  startingestimates <- c(-4.863319685,-0.009775098)
  fix <- list()
  fix[[1]] <- NULL
  fix[[2]] <- -100
  results.erpm <- estimate_ERPM(partition, nodes, objects, effects, 
                                startingestimates, multiplicationfactor, gainfactor, mini.steps, burnin, thining,
                                length.p1, min.iter.p2, max.iter.p2, num.steps.p2, length.p3, neighborhood, 
                                fixed.estimates=NULL, sizes.allowed = 2:5, sizes.simulated = 1:5)
    
}

barplots60nodes <- function(thetamin, thetamax, nodes,effects,objects,burnin,thining,num.steps,mini.steps) {
  
  effects <- list( names = c("num_groups","sizes_squared"),
                   objects = c("partition","partition"))
  theta <- c(-4.863319685,-0.009775098)
  sizes.allowed <- 2:5
  sizes.simulated <- 1:5
  
  # simulate a large sample with the estimates found in phase 2  
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
  
  # instantiate with the starting network
  current.partition <- first.partition
  current.z <- computeStatistics(current.partition, nodes, effects, objects)
  current.logit <- theta * current.z
  
  # store the distribution of the groups
  all.ngs <- c()
  #all.z <- c()
  
  end.walk <- FALSE
  cpt <- 0
  cpt2 <- 0
  cpt_thining <- 0
  num.steps <- 1000
  
  while(!end.walk){
    
    ## IF NEIGHBORHOOD IS: MERGES OR SPLITS OF 2 (Pi_1)
    if(neighborhood == 1 && is.null(sizes.allowed)) {
      new.partition <- sample_new_partition_p1(current.partition, mini.steps)
    } else if(neighborhood == 2 && is.null(sizes.allowed)){
      new.partition <- sample_new_partition_p2(current.partition, mini.steps)
    }
    if(neighborhood == 1 && !is.null(sizes.allowed)) {
      new.partition <- sample_new_partition_p1_restricted(current.partition, mini.steps, sizes.simulated)
    } else if(neighborhood == 2 && !is.null(sizes.allowed)){
      new.partition <- sample_new_partition_p2_restricted(current.partition, mini.steps, sizes.simulated)
    }
    
    # compute new statistics only if it changed
    if(!all(current.partition == new.partition)) {
      new.z <- computeStatistics(new.partition, nodes, effects, objects)
      #new.z <- computeChangeStatistics(current.z, current.partition, new.partition, old_g1, old_g2, new_g1, new_g2, nodes, effects, objects)
    } else {
      new.z <- current.z
    }
    new.logit <- theta * new.z
    
    # chose whether to change or not
    if(mini.steps == "normalized") {
      
      if(neighborhood == 1 && is.null(sizes.allowed)) {
        current.size <- compute_size_neighborhood_p1(current.partition)
        new.size <- compute_size_neighborhood_p1(new.partition)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else if(neighborhood == 2 && is.null(sizes.allowed)) {
        current.size <- compute_size_neighborhood_p2(current.partition)
        new.size <- compute_size_neighborhood_p2(new.partition)
        neighborhoods.ratio <- (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
      }
      if(neighborhood == 1 && !is.null(sizes.allowed)) {
        current.size <- compute_size_neighborhood_p1_restricted(current.partition, sizes.simulated)
        new.size <- compute_size_neighborhood_p1_restricted(new.partition, sizes.simulated)
        neighborhoods.ratio <- current.size$num.swaps / new.size$num.swaps 
      } else if(neighborhood == 2 && !is.null(sizes.allowed)) {
        current.size <- compute_size_neighborhood_p2_restricted(current.partition, sizes.simulated)
        new.size <- compute_size_neighborhood_p2_restricted(new.partition, sizes.simulated)
        neighborhoods.ratio <- (new.size$num.merges + new.size$num.divisions)/ (current.size$num.merges + current.size$num.divisions)
      }
      hastings.ratio <- (exp(sum(new.logit) - sum(current.logit))) * neighborhoods.ratio
      
    } else if(mini.steps == "selfloops"){
      
      hastings.ratio <- exp(sum(new.logit) - sum(current.logit))
      
    }
    proba.change <- min(1,hastings.ratio)
    change.made <- (runif(1,0,1) <= proba.change)
    
    if(change.made) {
      current.partition <- new.partition
      current.z <- new.z
      current.logit <- new.logit
    }
    
    cpt <- cpt + 1
    if(cpt > burnin) cpt_thining <- cpt_thining + 1
    
    # storing results if all sizes are allowed
    if(is.null(sizes.allowed)){
      # store the results if we are out of burnin
      if(cpt >= burnin && cpt_thining == thining) {
        ngs <- rep(0,num.nodes)
        for(g in 1:max(current.partition)) {
          ngs[sum(current.partition == g)] <- ngs[sum(current.partition == g)] + 1 
        }
        all.ngs <- rbind(all.ngs, ngs)
        cpt_thining <- 0
      }
      
      # stop the walk if number of steps reached 
      print(cpt)
      end.walk <- (cpt >= (burnin+thining*num.steps))
      print(end.walk)
    }
    
    # storing results if sizes are constrained
    if(!is.null(sizes.allowed)){
      
      #print(check_sizes(current.partition,sizes.allowed))
      #print(current.partition)
      
      # if the partition is sampled with the right sizes
      # we keep the swapping strategy
      if(check_sizes(current.partition,sizes.allowed)){
        #neighborhood <- 2
        
        # if wrong sizes, we continue searching, and if we are out of tolerated simulated partitions, we go back
        # when out of allowed partitions we keep the divide/merge strategy
      } else {
        #neighborhood <- 2
        if(!check_sizes(current.partition,sizes.simulated)){
          current.partition <- old.partition
          if(check_sizes(current.partition,sizes.allowed)) neighborhood <- 1
        } 
      }
      
      # store the results if we are out of burnin
      if(cpt >= burnin && cpt_thining == thining) {
        cpt_thining <- 0
        if(check_sizes(current.partition,sizes.allowed)){
          ngs <- rep(0,num.nodes)
          for(g in 1:max(current.partition)) {
            ngs[sum(current.partition == g)] <- ngs[sum(current.partition == g)] + 1 
          }
          all.ngs <- rbind(all.ngs, ngs)
          cpt2 <- cpt2 + 1
        } 
      }
      
      # stop the walk if number of steps reached 
      end.walk <- (cpt2 >= num.steps)
    }

  }
  
  # plot sizes
  df2 <- data.frame(group_size=1:6, count=colMeans(all.ngs[,1:6]))
  ggplot(data=df2, aes(x=group_size, y=count)) +
    geom_bar(stat="identity")
  
  
}

testburninthining <- function(){
  participants2017 <- read.csv(file="PolyHack/participants-data2017.csv")
  networks2017 <- readRDS(file="PolyHack/networks2017.rds")
  nodes <- data.frame(label = participants2017$key,
                      age = participants2017$age_imputed,
                      language = participants2017$language_imputed,
                      level = participants2017$current.degree,
                      major = participants2017$category)
  net_acquaintances2017 <- networks2017$known_before_imputed
  net_samelanguage2017 <- networks2017$same_language_imputed
  net_samemajor2017 <- networks2017$same_major
  net_colocation2017_total <- networks2017$
    num.nodes <- 60
  partition <- participants2017$Team
  effects <- list( names = c("num_groups","sizes_squared","diff","same","same","same","tie"),
                   objects = c("partition","partition","age","language","level","major","net_acquaintances2017"))
  objects <- list()
  objects[[1]] <- list(name = "net_acquaintances2017", object = net_acquaintances2017)
  startingestimates <- c(-4.325,-0.078,0.017,-0.294,0.328,-0.307,6.575)
  
  burninthining_n2 <- find_burninthining(partition, 
                                      theta = startingestimates, 
                                      nodes, 
                                      effects,
                                      objects,
                                      num.steps = 200,
                                      mini.steps = "normalized",
                                      neighborhood = 2, 
                                      sizes.allowed = 2:5, 
                                      sizes.simulated = 1:5,
                                      max.thining = 1000)
  burninthining_n3 <- find_burninthining(partition, 
                                         theta = startingestimates, 
                                         nodes, 
                                         effects,
                                         objects,
                                         num.steps = 200,
                                         mini.steps = "normalized",
                                         neighborhood = 3, 
                                         sizes.allowed = 2:5, 
                                         sizes.simulated = 1:5,
                                         max.thining = 1000)
  burninthining_n3_0.8 <- find_burninthining(partition, 
                                         theta = startingestimates, 
                                         nodes, 
                                         effects,
                                         objects,
                                         num.steps = 200,
                                         mini.steps = "normalized",
                                         neighborhood = 3, 
                                         sizes.allowed = 2:5, 
                                         sizes.simulated = 1:5,
                                         max.thining = 1000)
  burninthining_n3_0.7 <- find_burninthining(partition, 
                                             theta = startingestimates, 
                                             nodes, 
                                             effects,
                                             objects,
                                             num.steps = 200,
                                             mini.steps = "normalized",
                                             neighborhood = c(0.6,0.3,0.1), 
                                             sizes.allowed = 2:5, 
                                             sizes.simulated = 1:5,
                                             max.thining = 150)
  
  smoothedautocor <- burninthining_n3_0.7$autocorrelations
  for(eff in 1:7) {
    lo <- loess(y ~ x, data.frame(x=1:150,y=burninthining_n3_0.7$autocorrelations[,eff]))
    smoothedautocor[!is.na(smoothedautocor[,eff]),eff] <- lo$fitted
  }
  ggplot(data.frame(autocorr = as.vector(smoothedautocor),
                    thining = rep(1:150,7),
                    effect = c(rep("1",150),rep("2",150),rep("3",150),rep("4",150),rep("5",150),rep("6",150),rep("7",150)))) +
    geom_point(aes(x=thining,y=autocorr,color=effect)) +
    scale_color_discrete()
  
  ggplot(data.frame(autocorr = as.vector(burninthining_n3_0.7$moving.means),
                    thining = rep(1:200,7),
                    effect = c(rep("1",200),rep("2",200),rep("3",200),rep("4",200),rep("5",200),rep("6",200),rep("7",200)))) +
    geom_point(aes(x=thining,y=autocorr,color=effect)) +
    scale_color_discrete()
  
  neighborhoods <- list(c(0.6,0.3,0.1),c(0.5,0.2,0.3))#,c(0.3,0.1,0.6),c(0.1,0.3,0.6))
  gridsearch <- gridsearch_burninthining_single(partition, # observed partition
                                              theta = startingestimates, # initial model parameters
                                              nodes, # nodeset (data frame)
                                              effects, # effects/sufficient statistics (list with a vector "names", and a vector "objects")
                                              objects, # objects used for statistics calculation (list with a vector "name", and a vector "object")
                                              num.steps = 100, # number of samples wanted in phase 1
                                              mini.steps = "normalized", # type of transition, either "normalized", either "self-loops" (take "normalized")
                                              neighborhoods, # list of probability vectors (proba actors swap, proba merge/division, proba single actor move)
                                              sizes.allowed = 2:5, # vector of group sizes allowed in sampling (now, it only works for vectors like size_min:size_max)
                                              sizes.simulated = 1:5, # vector of group sizes allowed in the Markov chain but not necessraily sampled (now, it only works for vectors like size_min:size_max) 
                                              max.thining = 50, # where to stop adding thining
                                              parallel = F, # to run different neighborhoods in parallel
                                              cpus = 1)
    
}

estimate2017 <- function(){
  participants2017 <- read.csv(file="PolyHack/participants-data2017.csv")
  networks2017 <- readRDS(file="PolyHack/networks2017.rds")
  nodes <- data.frame(label = participants2017$key,
                      age = participants2017$age_imputed,
                      language = participants2017$language_imputed,
                      level = participants2017$current.degree,
                      major = participants2017$category)
  net_acquaintances2017 <- networks2017$known_before_imputed
  net_samelanguage2017 <- networks2017$same_language_imputed
  net_samemajor2017 <- networks2017$same_major
  net_colocation2017_total <- networks2017$
  num.nodes <- 60
  partition <- participants2017$Team
  effects <- list( names = c("num_groups","sizes_squared","diff","same","same","tie"),
                   objects = c("partition","partition","age","language","level","net_acquaintances2017"))
  effects <- list( names = c("num_groups","sizes_squared","diff","same","same","same","tie"),
                   objects = c("partition","partition","age","language","level","major","net_acquaintances2017"))
  objects <- list()
  objects[[1]] <- list(name = "net_acquaintances2017", object = net_acquaintances2017)
  multiplicationfactor <- 30
  gainfactor <- 0.1 
  mini.steps <- "normalized"
  burnin <- 500
  thining <- 1000
  length.p1 <- 500
  min.iter.p2 <- 200
  max.iter.p2 <- 300
  num.steps.p2 <- 20
  length.p3 <- 2000
  neighborhood <- 2
  startingestimates <- c(-5,0,0,0,0,0)
  startingestimates <- c(-4.325,-0.078,0.017,-0.294,0.328,-0.307,6.575)
  fix <- list()
  fix[[1]] <- NULL
  fix[[2]] <- -100
  results.erpm <- estimate_ERPM(partition, nodes, objects, effects, 
                                startingestimates, multiplicationfactor, gainfactor, mini.steps, burnin, thining,
                                length.p1, min.iter.p2, max.iter.p2, num.steps.p2, length.p3, neighborhood, 
                                fixed.estimates=NULL, sizes.allowed = 2:5, sizes.simulated = 1:5, double.averaging = F)
  
  library(gmp)
  estimate_logL()
  
}



estimate2018 <- function(){
  participants2018 <- read.csv(file="PolyHack/participants-data2018.csv")
  networks2018 <- readRDS(file="PolyHack/networks2018.rds")
  nodes <- data.frame(label = participants2018$key,
                      age = participants2018$age_imputed,
                      language = participants2018$language_imputed,
                      level = participants2018$current.degree,
                      major = participants2018$category)
  net_acquaintances2018 <- networks2018$known_before_imputed

  num.nodes <- 58
  partition <- participants2018$Team
  effects <- list( names = c("num_groups","sizes_squared","diff","same","same","same","tie"),
                   objects = c("partition","partition","age","language","level","major","net_acquaintances2018"))
  effects <- list( names = c("num_groups","diff","same","same","number_attributes","tie"),
                   objects = c("partition","age","language","level","major","net_acquaintances2018"))
  objects <- list()
  objects[[1]] <- list(name = "net_acquaintances2018", object = net_acquaintances2018)
  
  multiplicationfactor <- 30
  gainfactor <- 0.5
  mini.steps <- "normalized"
  burnin <- 500
  thining <- 10
  length.p1 <- 1000
  min.iter.p2 <- 200
  max.iter.p2 <- 300
  num.steps.p2 <- 20
  length.p3 <- 3000
  neighborhood <- 2
  startingestimates <- c(-1.0526, -0.0289, 0.5362, 0.04225, 0.0962, 2.4916)
  
  fix <- list()
  fix[[1]] <- NULL
  fix[[2]] <- -100
  
  results.erpm <- estimate_ERPM(partition, nodes, objects, effects, 
                                startingestimates, multiplicationfactor, gainfactor, mini.steps, burnin, thining,
                                length.p1, min.iter.p2, max.iter.p2, num.steps.p2, length.p3, neighborhood, 
                                fixed.estimates=NULL, sizes.allowed = 3:5, sizes.simulated = 1:5, double.averaging = F)
  

  startingestimates <- c(-1.0722,-0.0295,0.5350,0.0355,-0.0809,2.5003)
  results.phase3 <- estimate_ERPM_p3(partition, 
                               nodes, 
                               objects, 
                               effects, 
                               startingestimates, 
                               mini.steps = "normalized", 
                               burnin = 500, 
                               thining = 50,
                               length.p3 = 5000,
                               neighborhood = 2,
                               fixed.estimates = NULL,
                               sizes.allowed = 3:5,
                               sizes.simulated = 1:5)
  
  library(gmp)
  estimate_logL()
  
}




GOF2017 <- function(){
  
  participants2017 <- read.csv(file="PolyHack/participants-data2017.csv")
  networks2017 <- readRDS(file="PolyHack/networks2017.rds")
  nodes <- data.frame(label = participants2017$key,
                      age = participants2017$age_imputed,
                      language = participants2017$language_imputed,
                      level = participants2017$current.degree,
                      major = participants2017$category)
  net_acquaintances2017 <- networks2017$known_before_imputed
  net_samelanguage2017 <- networks2017$same_language_imputed
  net_samemajor2017 <- networks2017$same_major
  
  num.nodes <- 60
  partition <- participants2017$Team
  objects <- list()
  objects[[1]] <- list(name = "net_acquaintances2017", object = net_acquaintances2017)
  
  distlan <- as.matrix(dist(as.numeric(nodes$language)))
  distlan <- distlan == 0
  diag(distlan) <- 0
  distlev <- as.matrix(dist(as.numeric(nodes$level)))
  distlev <- distlev == 0
  diag(distlev) <- 0
  distmaj <- as.matrix(dist(as.numeric(nodes$major)))
  distmaj <- distmaj == 0
  diag(distmaj) <- 0
  
  # Model 1
  effects <- list( names = c("num_groups","sizes_squared","diff","same","same","same","tie"),
                   objects = c("partition","partition","age","language","level","major","net_acquaintances2017"))
  estimates <- c(-3.96,-0.04,-0.01,-0.15,0.20,-0.47,6.16)
  
  sims_M1 <- simulate_ERPM(theta = estimates,
                        nodes = nodes,
                        effects = effects, 
                        objects = objects,
                        partition = partition,
                        burnin = 500,
                        thining = 500,
                        num.steps = 2000,
                        mini.steps = "normalized",
                        neighborhood = 2,
                        sizes.allowed = 2:5,
                        sizes.simulated = 1:5) 
  
  # Model 1 - all violin plots
  theme_set(theme_minimal())
  numsims <- nrow(sims_M1$stats)
  
  # group sizes
  numsize2_1 <- c()
  numsize3_1 <- c()
  numsize4_1 <- c()
  numsize5_1 <- c()
  for(i in 1:numsims){
    p <- sims_M1$partitions[i,]
    numsize2_1 <- c(numsize2_1,length(which(table(sims_M1$partitions[i,])==2)))
    numsize3_1 <- c(numsize3_1,length(which(table(sims_M1$partitions[i,])==3)))
    numsize4_1 <- c(numsize4_1,length(which(table(sims_M1$partitions[i,])==4)))
    numsize5_1 <- c(numsize5_1,length(which(table(sims_M1$partitions[i,])==5)))
  }
  df <- data.frame(nsim = 1:numsims,
                   sim_stats=c(numsize2,numsize3,numsize4,numsize5),
                   size=c(rep(2,numsims),rep(3,numsims),rep(4,numsims),rep(5,numsims)),
                   obs_stat=c(rep(1,numsims),rep(1,numsims),rep(5,numsims),rep(7,numsims)),
                   mean_stat=c(rep(mean(numsize2),numsims),rep(mean(numsize3),numsims),rep(mean(numsize4),numsims),rep(mean(numsize5),numsims)))
  ggplot(df, aes(factor(size), sim_stats)) + 
    geom_violin() +
    geom_jitter(shape=16,position=position_jitter(0.3), color="grey", alpha=0.2) +
    stat_summary(fun.data="mean_sdl", geom="pointrange", size=0.5) +
    geom_point(aes(x=factor(size),y=mean_stat, colour = "sim", shape= "sim")) +
    geom_point(aes(x=factor(size),y=obs_stat, colour = "obs", shape= "obs"), stroke= 1.5) +
    labs(x = "Group size",
         y = "Group count",
         color="Legend",
         shape="Legend") +
    scale_color_manual(values = c(sim="black",obs="red"), labels=c("observed value","simulated average"))+
    scale_shape_manual(values = c(sim=19,obs=5), labels=c("observed value","simulated average")) +
    theme(legend.position = "right")
    
  
  # icc
  icc_1 <- c()
  for(i in 1:numsims){
    p <- sims_M1$partitions[i,]
    icc_1 <- c(icc_1,computeicc(nodes$age,p))
  }
  
  # densities
  denslan_1 <- c()
  denslev_1 <- c()
  densmaj_1 <- c()
  denstie_1 <- c()
  for(i in 1:numsims){
    p <- sims_M1$partitions[i,]
    denslan_1 <- c(denslan_1,computedensity(distlan,p)$average)
    denslev_1 <- c(denslev_1,computedensity(distlev,p)$average)
    densmaj_1 <- c(densmaj_1,computedensity(distmaj,p)$average)
    denstie_1 <- c(denstie_1,computedensity(net_acquaintances2017,p)$average)
  }
  
  # average number of majors
  avnummaj_1 <- c()
  for(i in 1:numsims){
    p <- sims_M1$partitions[i,]
    a <- unlist(lapply(1:max(p),
                       function(x){return(length(unique(nodes$major[which(p==x)])))}))
    avnummaj_1 <- c(avnummaj_1,mean(a))
  }
  
  # Model 2
  effects <- list( names = c("num_groups","sizes_squared","diff","same","same","count","tie"),
                   objects = c("partition","partition","age","language","level","major","net_acquaintances2017"))
  estimates <- c(-4.36,-0.11,-0.02,-0.11,0.25,0.32,5.71)
  
  sims_M2 <- simulate_ERPM(theta = estimates,
                        nodes = nodes,
                        effects = effects, 
                        objects = objects,
                        partition = partition,
                        burnin = 500,
                        thining = 500,
                        num.steps = 2000,
                        mini.steps = "normalized",
                        neighborhood = 2,
                        sizes.allowed = 2:5,
                        sizes.simulated = 1:5) 
  
  # Model 2 - all violin plots
  theme_set(theme_minimal())
  numsims <- nrow(sims_M2$stats)
  
  # group sizes
  numsize2_2 <- c()
  numsize3_2 <- c()
  numsize4_2 <- c()
  numsize5_2 <- c()
  for(i in 1:numsims){
    p <- sims_M2$partitions[i,]
    numsize2_2 <- c(numsize2_2,length(which(table(sims_M2$partitions[i,])==2)))
    numsize3_2 <- c(numsize3_2,length(which(table(sims_M2$partitions[i,])==3)))
    numsize4_2 <- c(numsize4_2,length(which(table(sims_M2$partitions[i,])==4)))
    numsize5_2 <- c(numsize5_2,length(which(table(sims_M2$partitions[i,])==5)))
  }
  
  # icc
  icc_2 <- c()
  for(i in 1:numsims){
    p <- sims_M2$partitions[i,]
    icc_2 <- c(icc_2,computeicc(nodes$age,p))
  }
  
  # densities
  denslan_2 <- c()
  denslev_2 <- c()
  densmaj_2 <- c()
  denstie_2 <- c()
  for(i in 1:numsims){
    p <- sims_M2$partitions[i,]
    denslan_2 <- c(denslan_2,computedensity(distlan,p)$average)
    denslev_2 <- c(denslev_2,computedensity(distlev,p)$average)
    densmaj_2 <- c(densmaj_2,computedensity(distmaj,p)$average)
    denstie_2 <- c(denstie_2,computedensity(net_acquaintances2017,p)$average)
  }
  
  # average number of majors
  avnummaj_2 <- c()
  for(i in 1:numsims){
    p <- sims_M2$partitions[i,]
    a <- unlist(lapply(1:max(p),
                       function(x){return(length(unique(nodes$major[which(p==x)])))}))
    avnummaj_2 <- c(avnummaj_2,mean(a))
  }
 
  df <- data.frame(nsim = 1:numsims,
                   sim_stats_1=c(numsize2_1,numsize3_1,numsize4_1,numsize5_1),
                   sim_stats_2=c(numsize2_2,numsize3_2,numsize4_2,numsize5_2),
                   size=c(rep(2,numsims),rep(3,numsims),rep(4,numsims),rep(5,numsims)),
                   obs_stat=c(rep(1,numsims),rep(1,numsims),rep(5,numsims),rep(7,numsims)),
                   mean_stat_1=c(rep(mean(numsize2_1),numsims),rep(mean(numsize3_1),numsims),rep(mean(numsize4_1),numsims),rep(mean(numsize5_1),numsims)),
                   mean_stat_2=c(rep(mean(numsize2_2),numsims),rep(mean(numsize3_2),numsims),rep(mean(numsize4_2),numsims),rep(mean(numsize5_2),numsims)))
  ggplot(df) + 
    geom_violin(aes(factor(size), sim_stats_1, linetype = "m1"), fill = NA) +
    stat_summary(aes(factor(size), sim_stats_1),fun.data="mean_sdl", geom="pointrange", size=0.5, shape=16) +
    geom_violin(aes(factor(size), sim_stats_2,  linetype = "m2"), fill = NA) +
    stat_summary(aes(factor(size), sim_stats_2),fun.data="mean_sdl", geom="pointrange", size=0.5, shape=15) +
    geom_point(aes(x=factor(size),y=mean_stat_1, colour = "sim1", shape= "sim1"), stroke= 1.5) +
    geom_point(aes(x=factor(size),y=mean_stat_2, colour = "sim2", shape= "sim2")) +
    geom_point(aes(x=factor(size),y=obs_stat, colour = "obs", shape= "obs"), stroke= 1.5) +
    labs(x = "Group size",
         y = "Group count",
         color="",
         shape="",
         linetype="") +
    scale_color_manual(values = c(obs="red", sim1="black", sim2="black"), labels=c("Observed value","Simulated average from M1","Simulated average from M2"))+
    scale_shape_manual(values = c(obs=5, sim1=16, sim2=15), labels=c("Observed value","Simulated average from M1","Simulated average from M2")) +
    scale_linetype_manual(values = c(m1 = 1, m2=2), labels=c("Simulated distribution for M1","Simulated distribution for M2")) +
    theme(legend.position = "right",
          legend.key.height = unit(0.4,"in"),
          text = element_text(size=20))
  
  
  df <- data.frame(nsim = 1:numsims,
                   sim_stats_1=icc_1,
                   sim_stats_2=icc_2,
                   stat=rep("age",numsims),
                   obs_stat=rep(0.1529943,numsims),
                   mean_stat_1=rep(mean(icc_1),numsims),
                   mean_stat_2=rep(mean(icc_2),numsims))
  ggplot(df) + 
    geom_violin(aes(factor(stat), sim_stats_1, linetype = "m1"), fill = NA) +
    stat_summary(aes(factor(stat), sim_stats_1),fun.data="mean_sdl", geom="pointrange", size=0.5, shape=16) +
    geom_violin(aes(factor(stat), sim_stats_2,  linetype = "m2"), fill = NA) +
    stat_summary(aes(factor(stat), sim_stats_2),fun.data="mean_sdl", geom="pointrange", size=0.5, shape=15) +
    geom_point(aes(x=factor(stat),y=mean_stat_1, colour = "sim1", shape= "sim1"), stroke= 1.5) +
    geom_point(aes(x=factor(stat),y=mean_stat_2, colour = "sim2", shape= "sim2")) +
    geom_point(aes(x=factor(stat),y=obs_stat, colour = "obs", shape= "obs"), stroke= 1.5) +
    labs(x = "",
         y = "Intraclass correlation coefficient (intragroup)",
         color="",
         shape="",
         linetype="") +
    scale_color_manual(values = c(obs="red", sim1="black", sim2="black"), labels=c("Observed value","Simulated average from M1","Simulated average from M2"))+
    scale_shape_manual(values = c(obs=5, sim1=16, sim2=15), labels=c("Observed value","Simulated average from M1","Simulated average from M2")) +
    scale_linetype_manual(values = c(m1 = 1, m2=2), labels=c("Simulated distribution for M1","Simulated distribution for M2")) +
    theme(legend.position = "right",
          legend.key.height = unit(0.4,"in"),
          text = element_text(size=20))
  
  
  df <- data.frame(nsim = 1:numsims,
                   sim_stats_1=c(denslan_1,denslev_1,densmaj_1,denstie_1),
                   sim_stats_2=c(denslan_2,denslev_2,densmaj_2,denstie_2),
                   stat=c(rep("language",numsims),rep("level",numsims),rep("major",numsims),rep("tie",numsims)),
                   obs_stat=c(rep(0.3119048,numsims),rep(0.4238095,numsims),rep(0.1833333,numsims),rep(0.202381,numsims)),
                   mean_stat_1=c(rep(mean(denslan_1),numsims),rep(mean(denslev_1),numsims),rep(mean(densmaj_1),numsims),rep(mean(denstie_1),numsims)),
                   mean_stat_2=c(rep(mean(denslan_2),numsims),rep(mean(denslev_2),numsims),rep(mean(densmaj_2),numsims),rep(mean(denstie_2),numsims)))
  ggplot(df) + 
    geom_violin(aes(factor(stat), sim_stats_1, linetype = "m1"), fill = NA) +
    stat_summary(aes(factor(stat), sim_stats_1),fun.data="mean_sdl", geom="pointrange", size=0.5, shape=16) +
    geom_violin(aes(factor(stat), sim_stats_2,  linetype = "m2"), fill = NA) +
    stat_summary(aes(factor(stat), sim_stats_2),fun.data="mean_sdl", geom="pointrange", size=0.5, shape=15) +
    geom_point(aes(x=factor(stat),y=mean_stat_1, colour = "sim1", shape= "sim1"), stroke= 1.5) +
    geom_point(aes(x=factor(stat),y=mean_stat_2, colour = "sim2", shape= "sim2")) +
    geom_point(aes(x=factor(stat),y=obs_stat, colour = "obs", shape= "obs"), stroke= 1.5) +
    labs(x = "Attribute",
         y = "Average density per group",
         color="",
         shape="",
         linetype="") +
    scale_color_manual(values = c(obs="red", sim1="black", sim2="black"), labels=c("Observed value","Simulated average from M1","Simulated average from M2"))+
    scale_shape_manual(values = c(obs=5, sim1=16, sim2=15), labels=c("Observed value","Simulated average from M1","Simulated average from M2")) +
    scale_linetype_manual(values = c(m1 = 1, m2=2), labels=c("Simulated distribution for M1","Simulated distribution for M2")) +
    theme(legend.position = "right",
          legend.key.height = unit(0.4,"in"),
          text = element_text(size=20))
  
  
  
  df <- data.frame(nsim = 1:numsims,
                   sim_stats_1=avnummaj_1,
                   sim_stats_2=avnummaj_2,
                   stat=rep("major",numsims),
                   obs_stat=rep(3,numsims),
                   mean_stat_1=rep(mean(avnummaj_1),numsims),
                   mean_stat_2=rep(mean(avnummaj_2),numsims))
  ggplot(df) + 
    geom_violin(aes(factor(stat), sim_stats_1, linetype = "m1"), fill = NA) +
    stat_summary(aes(factor(stat), sim_stats_1),fun.data="mean_sdl", geom="pointrange", size=0.5, shape=16) +
    geom_violin(aes(factor(stat), sim_stats_2,  linetype = "m2"), fill = NA) +
    stat_summary(aes(factor(stat), sim_stats_2),fun.data="mean_sdl", geom="pointrange", size=0.5, shape=15) +
    geom_point(aes(x=factor(stat),y=mean_stat_1, colour = "sim1", shape= "sim1"), stroke= 1.5) +
    geom_point(aes(x=factor(stat),y=mean_stat_2, colour = "sim2", shape= "sim2")) +
    geom_point(aes(x=factor(stat),y=obs_stat, colour = "obs", shape= "obs"), stroke= 1.5) +
    labs(x = "Attribute",
         y = "Average number of attributes per group",
         color="",
         shape="",
         linetype="") +
    scale_color_manual(values = c(obs="red", sim1="black", sim2="black"), labels=c("Observed value","Simulated average from M1","Simulated average from M2"))+
    scale_shape_manual(values = c(obs=5, sim1=16, sim2=15), labels=c("Observed value","Simulated average from M1","Simulated average from M2")) +
    scale_linetype_manual(values = c(m1 = 1, m2=2), labels=c("Simulated distribution for M1","Simulated distribution for M2")) +
    theme(legend.position = "right",
          legend.key.height = unit(0.4,"in"),
          text = element_text(size=20))
  
}


AIC2017 <- function() {
  
  participants2017 <- read.csv(file="PolyHack/participants-data2017.csv")
  networks2017 <- readRDS(file="PolyHack/networks2017.rds")
  nodes <- data.frame(label = participants2017$key,
                      age = participants2017$age_imputed,
                      language = participants2017$language_imputed,
                      level = participants2017$current.degree,
                      major = participants2017$category)
  net_acquaintances2017 <- networks2017$known_before_imputed
  net_samelanguage2017 <- networks2017$same_language_imputed
  net_samemajor2017 <- networks2017$same_major
  
  num.nodes <- 60
  partition <- participants2017$Team
  objects <- list()
  objects[[1]] <- list(name = "net_acquaintances2017", object = net_acquaintances2017)
  
  effects <- list( names = c("num_groups","sizes_squared","diff","same","same","same","tie"),
                   objects = c("partition","partition","age","language","level","major","net_acquaintances2017"))
  estimates_0 <- c(-4.06,0,0,0,0,0,0)
  estimates <- c(-3.96,-0.04,-0.01,-0.15,0.20,-0.47,6.16)
  
  aic_1 <- estimate_logL(partition,
                            nodes,
                            effects, 
                            objects,
                            theta = estimates,
                            theta_0 = estimates_0,
                            M = 40,
                            num.steps = 200,
                            burnin = 300,
                            thining = 10,
                            mini.steps = "normalized",
                            neighborhood = 2,
                            sizes.allowed = 2:5,
                            sizes.simulated = 1:5)
 
  
  effects <- list( names = c("num_groups","sizes_squared","diff","same","same","number_attributes","tie"),
                   objects = c("partition","partition","age","language","level","major","net_acquaintances2017"))
  estimates_0 <- c(-4.06,0,0,0,0,0,0)
  estimates <- c(-4.36,-0.11,-0.02,-0.11,0.25,0.32,5.71)
  
  aic_2 <- estimate_logL(partition,
                         nodes,
                         effects, 
                         objects,
                         theta = estimates,
                         theta_0 = estimates_0,
                         M = 40,
                         num.steps = 300,
                         burnin = 200,
                         thining = 10,
                         mini.steps = "normalized",
                         neighborhood = 2,
                         sizes.allowed = 2:5,
                         sizes.simulated = 1:5) 
}


AIC2018 <- function() {
  
  participants2018 <- read.csv(file="../PolyHack/participants-data2018.csv")
  networks2018 <- readRDS(file="../PolyHack/networks2018.rds")
  nodes <- data.frame(label = participants2018$key,
                      age = participants2018$age_imputed,
                      language = participants2018$language_imputed,
                      level = participants2018$current.degree,
                      major = participants2018$category)
  net_acquaintances2018 <- networks2018$known_before_imputed
  net_samelanguage2018 <- networks2018$same_language_imputed
  net_samemajor2018 <- networks2018$same_major
  
  num.nodes <- 58
  partition <- participants2018$Team
  objects <- list()
  objects[[1]] <- list(name = "net_acquaintances2018", object = net_acquaintances2018)
  
  effects <- list( names = c("num_groups","diff","same","same","same","tie"),
                   objects = c("partition","age","language","level","major","net_acquaintances2018"))
  estimates_0 <- c(-2.14,0,0,0,0,0)
  estimates <- c(-1.07,-0.03,0.53,0.03,-0.08,2.5)
  
  aic_1 <- estimate_logL(partition,
                         nodes,
                         effects, 
                         objects,
                         theta = estimates,
                         theta_0 = estimates_0,
                         M = 40,
                         num.steps = 300,
                         burnin = 10,
                         thining = 10,
                         mini.steps = "normalized",
                         neighborhood = 2,
                         sizes.allowed = 3:5,
                         sizes.simulated = 1:5)
  
}

library(ggplot2)
plot_final_graphs <- function() {
  
  # NUMBER OF GROUPS: ALL SIZES EVOLUTION
  allsizes <- readRDS("allsizes_numgroups_10nodes.rds")
  allsizes <- allsizes[5:37,]
  
  p1 <- seq(-8,8,0.5)
  #from MATLAB
  numgroups <-  c(1.148,1.224,1.327,1.453,1.595,1.741,1.886,2.031,2.187,2.366,2.580,2.835,3.132,3.478,3.876,4.333,4.850,
                  5.429,6.056,6.713,7.370,7.989,8.534,8.979,9.316,9.556,9.718,9.824,9.892,9.934,9.959,9.975,9.985)
  
  allsizes_smoothed <- allsizes
  for(i in 2:9){
    sizei <- allsizes[,i]
    lo <- predict(loess(sizei~p1),p1)
    plot(p1,sizei)
    lines(p1,lo, col="red", lwd=2)
    lo[lo<0] <- 0
    lo[lo>10] <- 10
    allsizes_smoothed[,i] <- lo
  }
  
  allsizes_smoothed <- cbind(allsizes_smoothed,numgroups)
  
  theme_set(theme_minimal())
  df <- data.frame(parameter=p1, group_count=as.vector(allsizes_smoothed),
                   group=rep(c(paste0("of size ", 1:10),"of any size"), each=length(p1)))
  df$group <- factor(df$group, levels = df$group)
  
  # legend right
  ggplot(data = df, aes(x=parameter, y=group_count)) + 
    geom_line(aes(colour=group,linetype=group), size = 1.2) + 
    scale_linetype_manual(name="Groups\n", values=c(rep(1,10),3)) +
    scale_color_manual(name="Groups\n", values=c("#530018","#a00021","#ca4a3d","#ef936f","#fcd2bb","#f5f5f5","#c7dfec","#81b8d6","#3680b6","#1b519c","#000000")) +
    xlab(expression("Parameter" ~ alpha)) +
    ylab("Expected number of groups") +
    theme(legend.key.size =  unit(c(0.4,0.5), "in"),
          axis.title.x = element_text(margin=margin(25,0,0,0)),
          axis.title.y = element_text(margin=margin(0,25,0,0)),
          text = element_text(size = 24)) 
  # legend bottom
  ggplot(data = df, aes(x=parameter, y=group_count)) + 
    geom_line(aes(colour=group,linetype=group), size = 1.2) + 
    scale_linetype_manual(name="Groups", values=c(rep(1,10),3)) +
    scale_color_manual(name="Groups", values=c("#530018","#a00021","#ca4a3d","#ef936f","#fcd2bb","#f5f5f5","#c7dfec","#81b8d6","#3680b6","#1b519c","#000000")) +
    xlab(expression("Parameter" ~ alpha)) +
    ylab("Expected number of groups") +
    theme(legend.position = "bottom",
          text = element_text(size = 14)) +guides(linetype=guide_legend(nrow=2,byrow=TRUE))
  
  
  # NUMBER OF GROUPS: 3 BARPLOTS FOR -2, 0, 2
  sizes1 <- allsizes[9,]
  sizes2 <- allsizes[13,]
  sizes3 <- allsizes[17,]
  
  df2 <- data.frame(size=rep(seq(1,10,1),3), size_count=c(sizes1,sizes2,sizes3),
                    group=c(rep("-4",10),rep("0",10),rep("2",10)))
  df2$group <- factor(df2$group, levels = df2$group)
  ggplot(data=df2, aes(x=size, y=size_count, fill=group)) +
    geom_bar(stat="identity", position=position_dodge())+
    scale_x_discrete(name="Group size", limits = seq(1,10,1)) +
    ylab("Expected number of groups") +
    scale_fill_brewer(name="Parameter\n", palette = 1, labels = c(expression(alpha~" = -4"),expression(alpha~" = -2"),expression(alpha~" = 0")))+
    theme(legend.key.size =  unit(c(0.4,0.4), "in"),
      axis.title.x = element_text(margin=margin(25,0,0,0)),
          axis.title.y = element_text(margin=margin(0,25,0,0)),
          text = element_text(size = 24),
          legend.text.align = 0)
  
  
  # NUMBER OF GROUPS: AVERAGE GROUP SIZE OF ONE NODE
  
  # with simulations
  expectedsizes <- rep(0,length(p1))
  for(i in 1:length(p1)){
    sizei <- allsizes_smoothed[i,]
    num <- 0
    den <- 0
    for(j in 1:10) {
      num <- num + allsizes_smoothed[i,j]*j^2
      den <- den + allsizes_smoothed[i,j]*j
    }
    expectedsizes[i] <- num/den
  }
  df3 <- data.frame(parameter=p1, size=expectedsizes)
  
  # calculated value
  library(gmp)
  probability_of_partition <- function(alpha,size,n){
    p <- exp(alpha*size) 
    all <- 0
    for(i in 1:n){
      all <- all + as.integer(Stirling2(n,i)) * exp(alpha*i) 
    }
    return(p/all)
  }
  
  expectedsizes <- rep(0,length(p1))
  for(i in 1:length(p1)){
    expsize <- 0
    for(j in 1:9) {
      p <- 0 
      for(k in 1:(10-j)) {
        p <- p + as.integer(Stirling2(n-j,k)) * probability_of_partition(p1[i],k+1,10)
      }
      p <- choose(9,j-1) * p
      expsize <- expsize + j*p
    }
    expsize <- expsize + 10*probability_of_partition(p1[i],1,10)
    expectedsizes[i] <- expsize
  }
  
  df3 <- data.frame(parameter=p1, size=expectedsizes)
  
  ggplot(data = df3, aes(x=parameter, y=size)) + 
    geom_line(size = 1.2, colour = "#3680b6") + 
    xlab(expression("Parameter" ~ alpha)) +
    ylab("Expected group size of a random node") +
    theme(axis.title.x = element_text(margin=margin(25,0,0,0)),
          axis.title.y = element_text(margin=margin(0,25,0,0)),
          text = element_text(size = 28))
  
  
  # NUMBER OF GROUPS: distribution of probabilities for one node to be in a group
  expectedpresences1 <- rep(0,10)
  expectedpresences2 <- rep(0,10)
  expectedpresences3 <- rep(0,10)
  den1 <- 0
  den2 <- 0
  den3 <- 0
  for(j in 1:10) {
    expectedpresences1[j] <- allsizes[9,j]*j
    expectedpresences2[j] <- allsizes[13,j]*j
    expectedpresences3[j] <- allsizes[17,j]*j
    den1 <- den1 + allsizes[9,j]*j
    den2 <- den2 + allsizes[13,j]*j
    den3 <- den3 + allsizes[17,j]*j
  }
  expectedpresences1 <- expectedpresences1/den1
  expectedpresences2 <- expectedpresences2/den2
  expectedpresences3 <- expectedpresences3/den3
  
  df4 <- data.frame(size=rep(seq(1,10,1),3), presence=c(expectedpresences1,expectedpresences2,expectedpresences3),
                    group=c(rep("-4",10),rep("-2",10),rep("0",10)))
  df4$group <- factor(df4$group, levels = df4$group)
  ggplot(data=df4, aes(x=size, y=presence, fill=group)) +
    geom_bar(stat="identity", position=position_dodge())+
    scale_x_discrete(name="Group size", limits = seq(1,10,1)) +
    ylab("Probability of the group size of a random node") +
    scale_fill_brewer(name="Parameter", palette = 1, labels = c(expression(alpha~" = -4"),expression(alpha~" = -2"),expression(alpha~" = 0")))+
    theme(legend.key.size =  unit(c(0.4,0.4), "in"),
          axis.title.x = element_text(margin=margin(25,0,0,0)),
          axis.title.y = element_text(margin=margin(0,25,0,0)),
          text = element_text(size = 24),
          legend.text.align = 0)
  
  
  
  # DISPERSION NUMBER OF GROUPS: 3 BARPLOTS FOR -2, 0, 2
  sizes1 <- c()
  sizes2 <- c()
  sizes3 <- c()
  
  probability_of_partition2a <- function(alpha1,alpha2,stat,n){
    p <- exp(alpha1*4 + alpha2*stat) 
    all <- normconstant2a(alpha1,alpha2,n)
    return(p/all)
  }
  normconstant2a <- function(alpha1,alpha2,n){
    if(n==0){
      k <- 1
    } else if(n==1){
      k <- exp(alpha1)
    } else {
      k <- 0
      for(i in 0:(n-1)) {
        k <- k + choose(n-1,i) * exp(alpha1) * factorial(n-i)^alpha2 * normconstant2a(alpha1,alpha2,i)                              
      }
    }
    return(k)
  }
  
  df2 <- data.frame(size=rep(seq(1,10,1),3), size_count=c(sizes1,sizes2,sizes3),
                    group=c(rep("-4",10),rep("0",10),rep("2",10)))
  df2$group <- factor(df2$group, levels = df2$group)
  ggplot(data=df2, aes(x=size, y=size_count, fill=group)) +
    geom_bar(stat="identity", position=position_dodge())+
    scale_x_discrete(name="Group size", limits = seq(1,10,1)) +
    ylab("Expected number of groups") +
    scale_fill_brewer(name="Parameter\n", palette = 1, labels = c(expression(alpha~" = -4"),expression(alpha~" = -2"),expression(alpha~" = 0")))+
    theme(legend.key.size =  unit(c(0.4,0.4), "in"),
          axis.title.x = element_text(margin=margin(25,0,0,0)),
          axis.title.y = element_text(margin=margin(0,25,0,0)),
          text = element_text(size = 24),
          legend.text.align = 0)
  
  
}