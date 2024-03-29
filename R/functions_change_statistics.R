######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Main Function to calculate sufficient statistics of partitions   ##
## Author: Marion Hoffman                                           ##
######################################################################


# Compute complete statistics for single partitions

#' Compute Statistics
#'
#' Function that computes the statistic vector for a given partition and a given model
#'
#' @param partition vector, A partition
#' @param nodes data frame, Node set
#' @param effects list with a vector "names", and a vector "objects", Effects/sufficient statistics
#' @param objects list with a vector "name", and a vector "object", Objects used for statistics calculation
#' @return the statistics
#' @importFrom stats dist rnorm qt pnorm
#' @export
computeStatistics <- function (partition, nodes, effects, objects){
  
  num.groups <- max(partition)
  num.nodes <- nrow(nodes)
  num.effects <- length(effects[[1]])
  statistics <- rep(0,num.effects)
  
  # create affiliation matrix for calculations
  affiliation<-matrix(0,nrow=num.nodes,ncol=num.groups)
  idx<-cbind(1:num.nodes,partition)
  affiliation[idx]<-1
  
  # calculate sizes
  sizes<-colSums(affiliation)
  
  # find isolates and groups
  isolates <- which(sizes==1)
  if(length(isolates) == 0) isolates <- NULL
  groups <- which(sizes>1)
  if(length(groups) == 0) groups <- NULL
  # isolates <- c()
  # groups <- c()
  # for(g in 1:num.groups){
  #   if(length(which(partition == g)) == 1) {
  #     isolates <- c(isolates,g)
  #   } else {
  #     groups <- c(groups,g)
  #   }
  # }
  
  
  adjacency <- affiliation %*% t(affiliation)
  
  num_dyads <- (sizes*(sizes - 1)) / 2
  num_dyads[num_dyads == 0] <- 1
  adjacency_norm <- affiliation %*% t(affiliation) / num_dyads[partition]
  diag(adjacency_norm) <- 0
  
  #weights <- unlist(lapply(sizes,function(x){return(rep(1/x,x))}))
  #weighted_adjacency <- adjacency * weights
  
  for(e in 1:num.effects) {
    
    effect.name <- effects$names[e]
    object.name <- effects$objects[e]
    object2.name <- effects$objects2[e]
    
    # --------- ISOLATES -----------
    if(effect.name == "isolates") {
      statistics[e] <- length(isolates)
    }
    
    # --------- NUM GROUPS -----------
    if(effect.name == "num_groups") {
      statistics[e] <- num.groups 
    }
    
    # --------- NUM TIES -----------
    # TODO DIVISER PAR DEUX
    if(effect.name == "num_ties") {
      sum <- 0
      for(g in 1:num.groups) { 
        sum <- sum + sizes[g]*(sizes[g]-1)
      }
      statistics[e] <- sum
    }
    
    # --------- NUM GROUPS 3 -----------
    if(effect.name == "num_groups_3") {
      sum <- 0
      for(g in 1:num.groups){ 
        size <- sizes[g]
        sum <- sum + (size == 3)
      }
      statistics[e] <- sum
    }
    
    # --------- NUM GROUPS 4 -----------
    if(effect.name == "num_groups_4") {
      sum <- 0
      for(g in 1:num.groups){
        size <- sizes[g]
        sum <- sum + (size == 4)
      }
      statistics[e] <- sum
    }
    
    # --------- NUM GROUPS 5 -----------
    if(effect.name == "num_groups_5") {
      sum <- 0
      for(g in 1:num.groups){
        size <- sizes[g]
        sum <- sum + (size == 5)
      }
      statistics[e] <- sum
    }
    
    # --------- NUM GROUPS 6 -----------
    if(effect.name == "num_groups_6") {
      sum <- 0
      for(g in 1:num.groups){ 
        size <- sizes[g] 
        sum <- sum + (size == 6)
      }
      statistics[e] <- sum
    }
    
    # --------- NUM GROUPS x PRESENT NODES -----------
    if(effect.name == "num_groups_x_num_nodes") {
      statistics[e] <- num.groups * num.nodes 
    }
    
    # --------- NUM GROUPS x LOG OF PRESENT NODES -----------
    if(effect.name == "num_groups_x_log_num_nodes") {
      statistics[e] <- num.groups * log(num.nodes) 
    }
    
    # --------- NUM TRIANGLES -----------
    if(effect.name == "num_triangles") {
      sum <- 0
      for(g in 1:num.groups){ 
        size <- sizes[g] 
        if(size == 3) sum <- sum + 1
        if(size > 3) sum <- sum + dim(combn(size,3))[2]
      }
      statistics[e] <- sum
    }
    
    # --------- NUM FOURS -----------
    if(effect.name == "num_fours") {
      sum <- 0
      for(g in 1:num.groups){ 
        size <- sizes[g] 
        if(size == 4) sum <- sum + 1
        if(size > 4) sum <- sum + dim(combn(size,4))[2]
      }
      statistics[e] <- sum
    }
    
    # --------- NUM FIVES -----------
    if(effect.name == "num_fives") {
      sum <- 0
      for(g in 1:num.groups){ 
        size <- sizes[g]
        if(size == 5) sum <- sum + 1
        if(size > 5) sum <- sum + dim(combn(size,5))[2]
      }
      statistics[e] <- sum
    }
    
    # --------- ALTERNATED CLIQUES -----------
    if(effect.name == "alt_cliques") {
      sum <- 0
      lambda <- 2
      allcounts <- rep(0,length(partition)-2)
      for(g in 1:num.groups){
        size <- sizes[g] 
        for(s in 3:length(partition)){
          if(size == s) allcounts[s-2] <- allcounts[s-2] + 1
          if(size > s) allcounts[s-2] <- allcounts[s-2] + dim(combn(size,s))[2]
        }
      }
      for(s in 3:length(partition)){
        sum <- sum + (-1)^(s+1) * allcounts[s-2] / (lambda^(s-3))
      }
      statistics[e] <- sum
    }
    
    # --------- SIZES_SQUARED -----------
    if(effect.name == "sizes_squared") {
      statistics[e] <- sum(sizes^2)
    }
    
    # --------- SIZES_SQUARED_NORM -----------
    if(effect.name == "sizes_squared_norm") {
      statistics[e] <- sum(sizes^2) / num.groups
    }
    
    # --------- PRODUCT SIZES -----------
    if(effect.name == "product_sizes") {
      product <- 1
      for(g in 1:num.groups){ 
        product <- product * sizes[g] 
      }
      statistics[e] <- product
    }
    
    # --------- SUM LOG FACTORIALS -----------
    if(effect.name == "sum_log_factorials") {
      sum <- 1
      for(g in 1:max(partition)){
        s <- sizes[g] 
        if(s>2) {
          sum <- sum + log( factorial(s-1) )
        }
      }
      statistics[e] <- sum
    }
    
    # --------- TIE -----------
    if(effect.name == "tie") {
      if(length(groups) > 0) {
        for(o in 1:length(objects)){
          if(objects[[o]][[1]] == object.name){
            net <- objects[[o]][[2]]
          }
        }
        statistics[e] <- sum(1/2 * adjacency * net)
      }else{
        statistics[e] <- 0
      }
      
    }
    
    
    # --------- TIE_X_DIFF -----------
    if(effect.name == "tie_X_diff") {
      if(length(groups) > 0) {
        for(o in 1:length(objects)){
          if(objects[[o]][[1]] == object.name){
            net <- objects[[o]][[2]]
          }
        }      
        #d <- as.matrix(dist(nodes[,object2.name]))
        d<-dist(nodes[,object2.name])
        d<-as.matrix(d)
        statistics[e] <- sum(1/2 * adjacency * net * d)
      }else{
        statistics[e] <- 0
      }
      
    }
    
    # --------- BIPARTITE TIE -----------
    if(effect.name == "bipartite_tie") {
      if(length(groups) > 0) {
        for(o in 1:length(objects)){
          if(objects[[o]][[1]] == object.name){
            binet <- objects[[o]][[2]]
          }
        }
        stat <- 0
        for(g in 1:dim(binet)[2]){
          members <- which(binet[,g] == 1)
          stat <- stat + sum(1/2 * adjacency[members,members])
        }
        statistics[e] <- stat
      }else{
        statistics[e] <- 0
      }
      
    }
    
    # --------- BIPARTITE GROUP -----------
    if(effect.name == "bipartite_group") {
      if(length(groups) > 0) {
        for(o in 1:length(objects)){
          if(objects[[o]][[1]] == object.name){
            binet <- objects[[o]][[2]]
          }
        }
        stat <- 0
        for(g in 1:dim(binet)[2]){
          members <- which(binet[,g] == 1)
          if(length(unique(partition[members])) == 1) {
            stat <- stat + 1
          }
        }
        statistics[e] <- stat
      }else{
        statistics[e] <- 0
      }
      
    }
    
    # --------- ATTRIBUTE ISOLATION -----------
    if(effect.name == "attisolation") {
      if(length(isolates) > 0) {
        att <- which(colnames(nodes) == object.name)
        sum.att <- 0
        for(g in isolates){
          member <- which(partition == g)
          sum.att <- sum.att + (nodes[member,att])
        }
        statistics[e] <- sum.att
      } else {
        statistics[e] <- 0
      }
    }
    
    # --------- ATTRIBUTE GROUPS -----------
    if(effect.name == "attgroups") {
      if(length(groups) > 0) {
        att <- which(colnames(nodes) == object.name)
        sum.att <- 0
        for(g in groups){
          members <- which(partition == g)
          sum.att <- sum.att + sum(nodes[members,att])
        }
        statistics[e] <- sum.att
      }else{
        statistics[e] <- 0
      }
    }
    
    # --------- ALTER -----------
    if(effect.name == "alter") {
      att <- which(colnames(nodes) == object.name)
      sum <- 0
      
      for(a in 1:num.nodes){
        sum <- sum + nodes[a,att]* sizes[partition[a]] 
      }
      
      statistics[e] <- sum
    }
    
    # --------- HOMOPHILY:SAME -----------
    if(effect.name == "same") {
      if(length(groups) > 0) {
        # num.same_total <- 0
        # for(g in groups){
        #   members <- which(partition == g)
        #   num.same <- 0
        #   for(i in 1:(length(members)-1)){
        #     for(j in (i+1):length(members)){
        #       num.same <- num.same + (nodes[members[i],att] == nodes[members[j],att]) 
        #     }
        #   }
        #   num.same_total <- num.same_total + num.same/length(members)
        # }
        # statistics[e] <- num.same_total
        #d <- as.matrix(dist(nodes[,object.name]))
        d<-dist(nodes[,object.name])
        d<-as.matrix(d)
        d <- d==0
        diag(d) <- 0
        statistics[e] <- sum(1/2 * adjacency * d)
      } else {
        statistics[e] <- 0
      }
      
    }
    
    # --------- HOMOPHILY:SAME NORMALIZED -----------
    if(effect.name == "same_norm") {
      if(length(groups) > 0) {
        #d <- as.matrix(dist(nodes[,object.name]))
        d<-dist(nodes[,object.name])
        d<-as.matrix(d)
        d <- d==0
        diag(d) <- 0
        statistics[e] <- sum(1/2 * adjacency_norm * d)
      } else {
        statistics[e] <- 0
      }
      
    }
    
    # --------- HOMOPHILY:DIFF -----------
    if(effect.name == "diff") {
      if(length(groups) > 0){
        # diff_total <- 0
        # for(g in groups){
        #   members <- which(partition == g)
        #   diff <- 0
        #   for(i in 1:(length(members)-1)){
        #     for(j in (i+1):length(members)){
        #       diff <- diff + abs(nodes[members[i],att] - nodes[members[j],att]) 
        #     }
        #   }
        #   diff_total <- diff_total + diff/length(members)
        # }
        # statistics[e] <- diff_total
        #d <- as.matrix(dist(nodes[,object.name]))
        d<-dist(nodes[,object.name])
        d<-as.matrix(d)
        statistics[e] <- sum(1/2 * adjacency * d)
      }else{
        statistics[e] <- 0
      }
      
    }
    
    # --------- HOMOPHILY:DIFF NORMALIZED -----------
    if(effect.name == "diff_norm") {
      if(length(groups) > 0){
        #d <- as.matrix(dist(nodes[,object.name]))
        d<-dist(nodes[,object.name])
        d<-as.matrix(d)
        stat <- stat + sum(1/2 * adjacency_norm * d)
      }else{
        statistics[e] <- 0
      }
    }
    
    # --------- HOMOPHILY:DIFF PER INDIVIDUAL -----------
    if(effect.name == "diff_ind") {
      if(length(groups) > 0){
        stat <- 0
        for(a in 1:num.nodes){
          g <- partition[a]
          others <- which(partition == g & 1:num.nodes != a)
          if(length(others) > 0) {
            diffs <- abs(nodes[a,object.name] - nodes[others,object.name])
            stat <- stat + min(diffs)
          }
        }
        statistics[e] <- stat
      }else{
        statistics[e] <- 0
      }
    }
    
    # --------- HOMOPHILY:DIFF PER INDIVIDUAL NORMALIZED-----------
    if(effect.name == "diff_ind_norm") {
      if(length(groups) > 0){
        for(a in 1:num.nodes){
          g <- partition[a]
          others <- which(partition == g & 1:num.nodes != a)
          if(length(others) > 0) {
            diffs <- abs(nodes[a,object.name] - nodes[others,object.name])
            stat <- stat + min(diffs) / (length(others)+1)
          }
        }
      }else{
        statistics[e] <- 0
      }
    }
    
    # --------- HOMOPHILY:RANGE -----------
    if(effect.name == "range") {
      if(length(groups) > 0){
        sum <- 0
        for(g in groups){
          members <- which(partition == g)
          sum <- sum + abs(max(nodes[members,object.name]) - min(nodes[members,object.name]))
        }
        statistics[e] <- sum
      }else{
        statistics[e] <- 0
      }
      
    }
    
    # --------- HOMOPHILY:PROPORTION -----------
    if(effect.name == "proportion") {
      if(length(groups) > 0){
        attribute <- as.numeric(factor(nodes[,object.name]))
        minatt <- min(attribute)
        maxatt <- max(attribute)
        sum <- 0
        for(g in groups){
          members <- which(partition == g)
          nmax <- sum(attribute[members] == maxatt)
          nmin <- sum(attribute[members] == minatt)
          sum <- sum + (max(nmax,nmin) - min(nmax,nmin)) / length(members)
        }
        statistics[e] <- sum
      }else{
        statistics[e] <- 0
      }
      
    }
    
    # --------- HOMOPHILY:NUM GROUPS SAME -----------
    if(effect.name == "num_groups_same") {
      if(length(groups) > 0){
        attribute <- as.numeric(factor(nodes[,object.name]))
        sum <- 0
        for(g in groups){
          members <- which(partition == g)
          sum <- sum + (length(unique(attribute[members])) == 1)
        }
        statistics[e] <- sum
      }else{
        statistics[e] <- 0
      }
      
    }
    
    # ---------GROUP: NUMBER_ATTRIBUTES -----------
    if(effect.name == "number_attributes") {
      if(length(groups) > 0) {
        # total_sum <- 0
        # lev <- levels(nodes[,att])
        # for(g in groups){
        #   members <- which(partition == g)
        #   allattributes <- c()
        #   sum <- 0
        #   for(i in 1:length(members)){
        #     l <- lev[nodes[members[i],att]]
        #     if(!l %in% allattributes){
        #       allattributes <- c(allattributes,l)
        #       sum <- sum + 1
        #     }
        #   }
        #   total_sum <- total_sum + sum/length(members)
        # }
        # statistics[e] <- total_sum
        d <- unlist(lapply(1:num.groups,
                           function(x){return(length(unique(nodes[which(partition==x),object.name])))}))
        statistics[e] <- sum(d)
      } else {
        statistics[e] <- 0
      }
      
    }
    
  }
  
  return(statistics)
  
}





# Compute complete statistics for multiple partitions

#' Compute Statistics multiple
#'
#' Function that computes the statistic vector for given (multiple) partitions and a given model
#'
#' @param partitions Observed partitions
#' @param presence.tables to indicate which nodes were present when
#' @param nodes Node set (data frame)
#' @param effects Effects/sufficient statistics (list with a vector "names", and a vector "objects")
#' @param objects Objects used for statistics calculation (list with a vector "name", and a vector "object")
#' @param single.obs equal NULL by default
#' @return A list
#' @importFrom stats dist
#' @export
computeStatistics_multiple <- function(partitions, presence.tables, nodes, effects, objects, single.obs = NULL){

  num.nodes <- nrow(nodes)
  num.obs <- ncol(presence.tables)
  num.effects <- length(effects[[1]])
  statistics <- matrix(0,num.effects,num.obs)

  # find isolates, groups, and sizes
  nums.groups <- rep(0,num.obs)
  isolates <- list()
  groups <- list()
  sizes <- list()
  for(o in 1:num.obs){
    p <- partitions[,o]
    p <- p[as.logical(presence.tables[,o])]

    nums.groups[o] <- max(p)
    isolates[[o]] <- as.vector(which(table(p)==1))
    groups[[o]] <- as.vector(which(table(p)>1))
    sizes[[o]] <- as.vector(table(p))
  }

  # create adjacency matrices
  adjacencies <- list()
  adjacencies_norm <- list()
  for(o in 1:num.obs){
    p <- partitions[,o]
    p <- p[as.logical(presence.tables[,o])]
    affiliation <- as.matrix(table(data.frame(actor = 1:length(p), group= p)))

    adjacencies[[o]] <- affiliation %*% t(affiliation)
    diag(adjacencies[[o]]) <- 0

    num_dyads <- (sizes[[o]]*(sizes[[o]] - 1)) / 2
    num_dyads[num_dyads == 0] <- 1
    adjacencies_norm[[o]] <- affiliation %*% t(affiliation) / num_dyads[p]
    diag(adjacencies_norm[[o]]) <- 0
  }



  for(e in 1:num.effects) {

    effect.name <- effects$names[e]
    object.name <- effects$objects[e]
    object2.name <- effects$objects2[e]

    # --------- ISOLATES -----------
    if(effect.name == "isolates") {
      for(o in 1:num.obs){
        statistics[e,o] <-length(isolates[[o]])
      }
    }

    # --------- NUM GROUPS -----------
    if(effect.name == "num_groups") {
      for(o in 1:num.obs){
        statistics[e,o] <- nums.groups[o]
      }
    }

    # --------- NUM GROUPS 3 -----------
    if(effect.name == "num_groups_3") {
      for(o in 1:num.obs){
        for(g in 1:nums.groups[o]){
          size <- length(which(partitions[,o]==g))
          statistics[e,o] <- (size == 3)
        }
      }
    }

    # --------- NUM GROUPS 4 -----------
    if(effect.name == "num_groups_4") {
      for(o in 1:num.obs){
        for(g in 1:nums.groups[o]){
          size <- length(which(partitions[,o]==g))
          statistics[e,o] <- (size == 4)
        }
      }
    }

    # --------- NUM GROUPS 5 -----------
    if(effect.name == "num_groups_5") {
      for(o in 1:num.obs){
        for(g in 1:nums.groups[o]){
          size <- length(which(partitions[,o]==g))
          statistics[e,o] <- (size == 5)
        }
      }
    }

    # --------- NUM GROUPS 6 -----------
    if(effect.name == "num_groups_6") {
      for(o in 1:num.obs){
        for(g in 1:nums.groups[o]){
          size <- length(which(partitions[,o]==g))
          statistics[e,o] <- (size == 6)
        }
      }
    }

    # --------- NUM GROUPS x PRESENT NODES -----------
    if(effect.name == "num_groups_x_num_nodes") {
      for(o in 1:num.obs){
        statistics[e,o] <- nums.groups[o] * sum(presence.tables[,o])
      }
    }

    # --------- NUM GROUPS x LOG OF PRESENT NODES -----------
    if(effect.name == "num_groups_x_log_num_nodes") {
      for(o in 1:num.obs){
        statistics[e,o] <- nums.groups[o] * log(sum(presence.tables[,o]))
      }
    }
    
    # --------- SIZES_SQUARED -----------
    if(effect.name == "sizes_squared") {
      for(o in 1:num.obs){
        statistics[e,o] <- sum(sizes[[o]]^2)
      }
    }

    # --------- SIZES_SQUARED_NORM -----------
    if(effect.name == "sizes_squared_norm") {
      for(o in 1:num.obs){
        statistics[e,o] <- sum(sizes[[o]]^2) / nums.groups[o]
      }
    }
    
    # --------- TIE -----------
    if(effect.name == "tie") {
      for(ob in 1:length(objects)){
        if(objects[[ob]][[1]] == object.name){
          net <- objects[[ob]][[2]]
        }
      }
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0) {
          net2 <- net[as.logical(presence.tables[,o]),as.logical(presence.tables[,o])]
          statistics[e,o] <- sum(1/2 * adjacencies[[o]] * net2)
        }
      }
    }

    # --------- TIE X DIFF -----------
    if(effect.name == "tie_X_diff") {
      for(ob in 1:length(objects)){
        if(objects[[ob]][[1]] == object.name){
          net <- objects[[ob]][[2]]
        }
      }
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0) {
          att <- which(colnames(nodes) == object2.name)
          d <- as.matrix(dist(nodes[,att]))
          d <- d[as.logical(presence.tables[,o]),as.logical(presence.tables[,o])]
          net2 <- net[as.logical(presence.tables[,o]),as.logical(presence.tables[,o])]
          statistics[e,o] <- sum(1/2 * adjacencies[[o]] * net2 * d)
        }
      }
    }

    # --------- VARIABLE TIE -----------
    if(effect.name == "tie_var") {
      for(ob in 1:length(objects)){
        if(objects[[ob]][[1]] == object.name){
          nets <- objects[[ob]][[2]]
        }
      }
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0) {
          net <- nets[[o]]
          net2 <- net[as.logical(presence.tables[,o]),as.logical(presence.tables[,o])]
          statistics[e,o] <- sum(1/2 * adjacencies[[o]] * net2)
        }
      }
    }

    # --------- VARIABLE TIE X DIFF -----------
    if(effect.name == "tie_var_X_diff") {
      for(ob in 1:length(objects)){
        if(objects[[ob]][[1]] == object.name){
          nets <- objects[[ob]][[2]]
        }
      }
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0) {
          att <- which(colnames(nodes) == object2.name)
          d <- as.matrix(dist(as.numeric(nodes[,att])))
          d <- d[as.logical(presence.tables[,o]),as.logical(presence.tables[,o])]
          net <- nets[[o]]
          net2 <- net[as.logical(presence.tables[,o]),as.logical(presence.tables[,o])]
          statistics[e,o] <- sum(1/2 * adjacencies[[o]] * net2 * d)
        }
      }
    }

    # --------- BIPARTITE TIE -----------
    if(effect.name == "bipartite_tie") {
      for(o in 1:length(objects)){
        if(objects[[o]][[1]] == object.name){
          binet <- objects[[o]][[2]]
        }
      }
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0) {
          binet2 <- binet[as.logical(presence.tables[,o]),]
          for(g in 1:dim(binet2)[2]){
            members <- which(binet2[,g] == 1)
            statistics[e,o] <- statistics[e,o] + sum(1/2 * adjacencies[[o]][members,members])
          }
        }
      }
    }

    # --------- BIPARTITE GROUP -----------
    if(effect.name == "bipartite_group") {
      for(o in 1:length(objects)){
        if(objects[[o]][[1]] == object.name){
          binet <- objects[[o]][[2]]
        }
      }
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0) {
          for(g in 1:dim(binet)[2]){
            members <- which(binet[,g] == 1)
            pmembers <- partitions[members,o]
            if(length(unique(pmembers[!is.na(pmembers)])) == 1) {
              statistics[e,o] <- statistics[e,o] + 1
            }
          }
        }
      }
    }

    # --------- ATTRIBUTE ISOLATION -----------
    if(effect.name == "attisolation") {
      for(o in 1:num.obs){
        if(length(isolates[[o]]) > 0) {
          att <- which(colnames(nodes) == object.name)
          for(g in isolates[[o]]){
            member <- which(partitions[,o] == g)
            statistics[e,o] <- statistics[e,o] + (nodes[member,att])
          }
        }
      }
    }
    
    # --------- ALTER -----------
    if(effect.name == "alter") {
      att <- which(colnames(nodes) == object.name)
      for(o in 1:num.obs){
        for(a in which(presence.tables[,o]==1)){
          statistics[e,o] <- statistics[e,o] + nodes[a,att]*length(which(partitions[,o] == partitions[a,o]))
        }
      }
    }

    # --------- HOMOPHILY:SAME -----------
    if(effect.name == "same") {
      att <- which(colnames(nodes) == object.name)
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0) {
          att.nodes <- factor(nodes[as.logical(presence.tables[,o]),att])
          d <- as.matrix(dist(as.numeric(att.nodes)))
          d <- d==0
          diag(d) <- 0
          statistics[e,o] <- sum(1/2 * adjacencies[[o]] * d)
        }
      }
    }

    # --------- HOMOPHILY:SAME_VAR -----------
    if(effect.name == "same_var") {
      for(ob in 1:length(objects)){
        if(objects[[ob]][[1]] == object.name){
          atts <- objects[[ob]][[2]]
        }
      }
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0) {
          att.nodes <- atts[as.logical(presence.tables[,o]),o]
          d <- as.matrix(dist(as.numeric(att.nodes)))
          d <- d==0
          diag(d) <- 0
          statistics[e,o] <- sum(1/2 * adjacencies[[o]] * d)
        }
      }
    }

    # --------- HOMOPHILY:SAME NORMALIZED -----------
    if(effect.name == "same_norm") {
      att <- which(colnames(nodes) == object.name)
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0) {
          att.nodes <- factor(nodes[as.logical(presence.tables[,o]),att])
          d <- as.matrix(dist(as.numeric(att.nodes)))
          d <- d==0
          diag(d) <- 0
          statistics[e,o] <- sum(1/2 * adjacencies_norm[[o]] * d)
        }
      }
    }

    # --------- HOMOPHILY:DIFF -----------
    if(effect.name == "diff") {
      att <- which(colnames(nodes) == object.name)
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0){
          d <- as.matrix(dist(as.numeric(nodes[as.logical(presence.tables[,o]),att])))
          statistics[e,o] <- sum(1/2 * adjacencies[[o]] * d)
        }
      }
    }

    # --------- HOMOPHILY:DIFF NORMALIZED -----------
    if(effect.name == "diff_norm") {
      att <- which(colnames(nodes) == object.name)
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0){
          d <- as.matrix(dist(as.numeric(nodes[as.logical(presence.tables[,o]),att])))
          statistics[e,o] <- sum(1/2 * adjacencies_norm[[o]] * d)
        }
      }
    }

    # --------- HOMOPHILY:DIFF PER INDIVIDUAL -----------
    if(effect.name == "diff_ind") {
      att <- which(colnames(nodes) == object.name)
      for(o in 1:num.obs){
        for(a in 1:num.nodes){
          g <- partitions[a,o]
          others <- which(partitions[-a,o] == g)
          if(length(others) > 0) {
            diffs <- abs(as.numeric(nodes[a,att]) - as.numeric(nodes[others,att]))
            statistics[e,o] <- statistics[e,o] + min(diffs)
          }
        }
      }
    }

    # --------- HOMOPHILY:DIFF PER INDIVIDUAL NORMALIZED -----------
    if(effect.name == "diff_ind_norm") {
      att <- which(colnames(nodes) == object.name)
      for(o in 1:num.obs){
        for(a in 1:num.nodes){
          g <- partitions[a,o]
          others <- which(partitions[-a,o] == g)
          if(length(others) > 0) {
            diffs <- abs(as.numeric(nodes[a,att]) - as.numeric(nodes[others,att]))
            statistics[e,o] <- statistics[e,o] + min(diffs) / (length(others)+1)
          }
        }
      }
    }


    # --------- HOMOPHILY:RANGE -----------
    if(effect.name == "range") {
      att <- which(colnames(nodes) == object.name)
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0) {
          for(g in groups[[o]]){
            members <- which(partitions[,o] == g)
            statistics[e,o] <- statistics[e,o] + abs(max(nodes[members,att]) - min(nodes[members,att]))
          }
        }
      }
    }
    
    # --------- HOMOPHILY:PROPORTION -----------
    if(effect.name == "proportion") {
      att <- which(colnames(nodes) == object.name)
      attribute <- as.numeric(factor(nodes[,att]))
      minatt <- min(attribute)
      maxatt <- max(attribute)
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0) {
          for(g in groups[[o]]){
            members <- which(partitions[,o] == g)
            nmax <- sum(attribute[members] == maxatt)
            nmin <- sum(attribute[members] == minatt)
            statistics[e,o] <- statistics[e,o] + (max(nmax,nmin) - min(nmax,nmin)) / length(members)
          }
        }
      }
    }

    # --------- HOMOPHILY:NUM GROUPS SAME -----------
    if(effect.name == "num_groups_same") {
      att <- which(colnames(nodes) == object.name)
      attribute <- as.numeric(factor(nodes[,att]))
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0){
          for(g in groups[[o]]){
            members <- which(partitions[,o] == g)
            statistics[e,o] <- statistics[e,o] + (length(unique(attribute[members])) == 1)
          }
        }
      }
    }


    # ---------GROUP: NUMBER_ATTRIBUTES -----------
    if(effect.name == "number_attributes") {
      att <- which(colnames(nodes) == object.name)
      for(o in 1:num.obs){
        if(length(groups[[o]]) > 0) {
          d <- unlist(lapply(1:nums.groups[o],
                             function(x){return(length(unique(nodes[which(partitions[,o]==x),att])))}))
          statistics[e,o] <- sum(d)
        }
      }
    }


    # --------- INERTIA MINUS 1 -----------
    if(effect.name == "inertia_1") {
      for(o in 2:num.obs){
        if(length(groups[[o]]) > 0) {
          adj1 <- matrix(0,num.nodes,num.nodes)
          adj2 <- matrix(0,num.nodes,num.nodes)
          adj1[as.logical(presence.tables[,o-1]),as.logical(presence.tables[,o-1])] <- adjacencies[[o-1]]
          adj2[as.logical(presence.tables[,o]),  as.logical(presence.tables[,o])]   <- adjacencies[[o]]
          adj12 <- adj1 * adj2
          statistics[e,o] <- sum(1/2 * adj12)
        }
      }
    }

    # --------- INERTIA TOTAL -----------
    if(effect.name == "inertia_total") {
      adj_total <- matrix(0,num.nodes,num.nodes)
      adj_total[as.logical(presence.tables[,1]),as.logical(presence.tables[,1])] <- adjacencies[[1]]
      for(o in 2:num.obs){
        if(length(groups[[o]]) > 0) {
          adj2 <- matrix(0,num.nodes,num.nodes)
          adj2[as.logical(presence.tables[,o]),  as.logical(presence.tables[,o])]   <- adjacencies[[o]]
          adj_total <- adj_total + adj2
          statistics[e,o] <- sum(1/2 * adj_total)
        }
      }
    }

    # --------- INTERACTION: INERTIA TOTAL X DIFF -----------
    if(effect.name == "inertia_total_X_diff") {
      att <- which(colnames(nodes) == object2.name)
      adj_total <- matrix(0,num.nodes,num.nodes)
      adj_total[as.logical(presence.tables[,1]),as.logical(presence.tables[,1])] <- adjacencies[[1]]
      for(o in 2:num.obs){
        if(length(groups[[o]]) > 0) {
          d <- as.matrix(dist(as.numeric(nodes[,att])))
          adj2 <- matrix(0,num.nodes,num.nodes)
          adj2[as.logical(presence.tables[,o]),  as.logical(presence.tables[,o])]   <- adjacencies[[o]]
          adj_total <- adj_total + adj2
          statistics[e,o] <- sum(1/2 * adj_total * d)
        }
      }
    }

  }

  return(statistics)

}

