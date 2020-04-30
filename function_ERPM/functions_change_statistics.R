# Computes complete statistics for all partitions
computeStatistics <- function (partition, nodes, effects, objects){
  
  num.groups <- max(partition)
  num.nodes = nrow(nodes)
  num.effects <- length(effects[[1]])
  statistics <- rep(0,num.effects)
  
  # find isolates and groups
  isolates <- as.vector(which(table(partition)==1))
  if(length(isolates) == 0) isolates <- NULL
  groups <- as.vector(which(table(partition)>1))
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
  
  # calculate sizes
  sizes <- as.vector(table(partition))
  
  # create affiliation matrix for calculations
  affiliation <- as.matrix(table(data.frame(actor = 1:num.nodes, group= partition)))
  adjacency <- affiliation %*% t(affiliation)
  weights <- unlist(lapply(sizes,function(x){return(rep(1/x,x))}))
  weighted_adjacency <- adjacency * weights
  
  for(e in 1:num.effects) {
    
    effect.name <- effects$names[e]
    object.name <- effects$objects[e]
    
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
      for(g in 1:max(partition)){
        sum <- sum + length(which(partition==g))*(length(which(partition==g))-1)
      }
      statistics[e] <- sum
    }
    
    # --------- NUM GROUPS 3 -----------
    if(effect.name == "num_groups_3") {
      sum <- 0
      for(g in 1:max(partition)){
        size <- length(which(partition==g))
        sum <- sum + (size == 3)
      }
      statistics[e] <- sum
    }
    
    # --------- NUM GROUPS 4 -----------
    if(effect.name == "num_groups_4") {
      sum <- 0
      for(g in 1:max(partition)){
        size <- length(which(partition==g))
        sum <- sum + (size == 4)
      }
      statistics[e] <- sum
    }
    
    # --------- NUM GROUPS 5 -----------
    if(effect.name == "num_groups_5") {
      sum <- 0
      for(g in 1:max(partition)){
        size <- length(which(partition==g))
        sum <- sum + (size == 5)
      }
      statistics[e] <- sum
    }
    
    # --------- NUM TRIANGLES -----------
    if(effect.name == "num_triangles") {
      sum <- 0
      for(g in 1:max(partition)){
        size <- length(which(partition==g))
        if(size == 3) sum <- sum + 1
        if(size > 3) sum <- sum + dim(combn(size,3))[2]
      }
      statistics[e] <- sum
    }
    
    # --------- NUM FOURS -----------
    if(effect.name == "num_fours") {
      sum <- 0
      for(g in 1:max(partition)){
        size <- length(which(partition==g))
        if(size == 4) sum <- sum + 1
        if(size > 4) sum <- sum + dim(combn(size,4))[2]
      }
      statistics[e] <- sum
    }
    
    # --------- NUM FIVES -----------
    if(effect.name == "num_fives") {
      sum <- 0
      for(g in 1:max(partition)){
        size <- length(which(partition==g))
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
      for(g in 1:max(partition)){
        size <- length(which(partition==g))
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
    
    # --------- DEGREE2 -----------
    if(effect.name == "degree2") {
      sum <- 0
      for(a in 1:length(partition)){
        sum <- sum + length(which(partition == partition[a]))^2
      }
      statistics[e] <- sum
    }
    
    # --------- AV_DEGREE -----------
    if(effect.name == "av_degree") {
      sum <- 0
      for(a in 1:length(partition)){
        sum <- sum + length(which(partition == partition[a]))
      }
      statistics[e] <- sum/num.nodes
    }
    
    # --------- AV_DEGREE2 -----------
    if(effect.name == "av_degree2") {
      sum <- 0
      for(a in 1:length(partition)){
        sum <- sum + length(which(partition == partition[a]))^2
      }
      statistics[e] <- sum/num.nodes
    }
    
    # --------- SIZE -----------
    if(effect.name == "size") {
      statistics[e] <- num.nodes/(length(groups) + length(isolates))
    }
    
    # --------- SIZE2 -----------
    if(effect.name == "size2") {
      sum <- 0
      for(g in 1:max(partition)){
        sum <- sum + length(which(partition==g))^2
      }
      statistics[e] <- sum/(length(groups) + length(isolates))
    }
    
    # --------- SIZE GROUPS -----------
    if(effect.name == "size_groups") {
      if(length(groups) > 0){
        statistics[e] <- (num.nodes-length(isolates))/length(groups)
      } else {
        statistics[e] <- 0
      }
    }
    
    # --------- PRODUCT SIZES -----------
    if(effect.name == "product_sizes") {
      product <- 1
      for(g in 1:max(partition)){
        product <- product * length(which(partition==g))
      }
      statistics[e] <- product
    }
    
    # --------- SUM LOG FACTORIALS -----------
    if(effect.name == "sum_log_factorials") {
      sum <- 1
      for(g in 1:max(partition)){
        s <- length(which(partition==g))
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
        # num.ties_total <- 0
        # for(g in groups){
        #   members <- which(partition == g)
        #   num.ties <- 0
        #   for(i in 1:(length(members)-1)){
        #     for(j in (i+1):length(members)){
        #       num.ties <- num.ties + net[members[i],members[j]] + net[members[j],members[i]]
        #     }
        #   }
        #   num.ties_total <- num.ties_total + num.ties/length(members)
        # }
        # statistics[e] <- num.ties_total
        statistics[e] <- sum(1/2 * adjacency * net)
      }else{
        statistics[e] <- 0
      }
      
    }
    
    # --------- ATTRIBUTE ISOLATION -----------
    if(effect.name == "attisolation") {
      if(length(isolates) > 0) {
        att <- which(colnames(nodes) == object.name)
        mean.att <- 0
        for(g in isolates){
          member <- which(partition == g)
          mean.att <- mean.att + (nodes[member,att])
        }
        statistics[e] <- mean.att/length(isolates)
      } else {
        statistics[e] <- 0
      }
    }
    
    # --------- ATTRIBUTE GROUPS -----------
    if(effect.name == "attgroups") {
      if(length(groups) > 0) {
        att <- which(colnames(nodes) == object.name)
        mean.att <- 0
        for(g in groups){
          members <- which(partition == g)
          mean.att <- mean.att + sum(nodes[members,att])
        }
        statistics[e] <- mean.att/(num.nodes-length(isolates))
      }else{
        statistics[e] <- 0
      }
    }
    
    # --------- ALTER -----------
    if(effect.name == "alter") {
      att <- which(colnames(nodes) == object.name)
      sum <- 0
      
      for(a in 1:num.nodes){
        sum <- sum + nodes[a,att]*length(which(partition == partition[a]))
      }
      
      statistics[e] <- sum/num.nodes
    }
    
    # --------- HOMOPHILY:SAME -----------
    if(effect.name == "same") {
      if(length(groups) > 0) {
        att <- which(colnames(nodes) == object.name)
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
        d <- as.matrix(dist(as.numeric(nodes[,att])))
        d <- d==0
        diag(d) <- 0
        statistics[e] <- sum(1/2 * adjacency * d)
      } else {
        statistics[e] <- 0
      }
     
    }
    
    # --------- HOMOPHILY:DIFF -----------
    if(effect.name == "diff") {
      if(length(groups) > 0){
        att <- which(colnames(nodes) == object.name)
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
        d <- as.matrix(dist(as.numeric(nodes[,att])))
        statistics[e] <- sum(1/2 * adjacency * d)
      }else{
        statistics[e] <- 0
      }
      
    }
    
    
    # ---------GROUP: NUMBER_ATTRIBUTES -----------
    if(effect.name == "number_attributes") {
      if(length(groups) > 0) {
        att <- which(colnames(nodes) == object.name)
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
                           function(x){return(length(unique(nodes[which(partition==x),att])))}))
        statistics[e] <- sum(d)
      } else {
        statistics[e] <- 0
      }
      
    }
    
  }
  
  return(statistics)
  
}



# Computes change statistics when merge or division
computeChangeStatistics <- function (old.z, old.partition, new.partition, old_g1, old_g2, new_g1, new_g2, nodes, effects, objects){
  
  num.groups <- max(new.partition)
  num.nodes = nrow(nodes)
  num.effects <- length(effects[[1]])
  statistics <- rep(0,num.effects)

  # define type of change
  if(old_g2 == 0) change <- "division"
  if(new_g2 == 0) change <- "merge"
  
  # find isolates and groups
  isolates <- c()
  groups <- c()
  for(g in 1:num.groups){
    if(length(which(new.partition == g)) == 1) {
      isolates <- c(isolates,g)
    } else {
      groups <- c(groups,g)
    }
  }
  
  for(e in 1:num.effects) {
    
    effect.name <- effects$names[e]
    object.name <- effects$objects[e]
    
    # --------- ISOLATES -----------
    if(effect.name == "isolates") {
      statistics[e] <- length(isolates)
    }
    
    # --------- NUM GROUPS -----------
    if(effect.name == "num_groups") {
      statistics[e] <- length(groups) 
    }
    
    # --------- SIZE -----------
    if(effect.name == "size") {
      statistics[e] <- num.nodes/(length(groups) + length(isolates))
    }
    
    # --------- SIZE GROUPS -----------
    if(effect.name == "size_groups") {
      if(length(groups) > 0){
        statistics[e] <- (num.nodes-length(isolates))/length(groups)
      } else {
        statistics[e] <- 0
      }
    }
    
    # --------- TIE -----------
    if(effect.name == "tie") {
      if(length(groups) > 0) {
        for(o in 1:length(objects)){
         if(objects[[o]][1] == object.name){
            net <- objects[[o]][2]
          }
        }
        if(change == "merge"){
          num.ties <- old.z[e]*((num.nodes-length(isolates))+1)
          members1 <- which(old.partition == old_g1)
          members2 <- which(old.partition == old_g2)
          for(m1 in members1){
            for(m2 in members2){
              num.ties <- num.ties + net[m1,m2] + net[m2,m1]
            }
          }
        }else if(change == "division"){
          num.ties <- old.z[e]*((num.nodes-length(isolates))-1)
          members1 <- which(new.partition == new_g1)
          members2 <- which(new.partition == new_g2)
          for(m1 in members1){
            for(m2 in members2){
              num.ties <- num.ties - net[m1,m2] - net[m2,m1]
            }
          }
        }
        statistics[e] <- num.ties/(num.nodes-length(isolates))
      } else {
        statistics[e] <- 0
      }
      
    }
    
    # --------- ATTRIBUTE ISOLATION -----------
    if(effect.name == "attisolation") {
      if(length(isolates) > 0) {
        att <- which(colnames(nodes) == object.name)
        if(change == "merge"){
          m_old1 <- which(old.partition == old_g1)
          m_old2 <- which(old.partition == old_g2)
          mean.att <- old.z[e]*(length(isolates) + (length(m_old1)==1) + (length(m_old2)==1))
          if(length(m_old1) == 1) mean.att <- mean.att - nodes[m_old1,att] 
          if(length(m_old2) == 1) mean.att <- mean.att - nodes[m_old2,att] 
        }else if(change == "division"){
          m_new1 <- which(new.partition == new_g1)
          m_new2 <- which(new.partition == new_g2)
          mean.att <- old.z[e]*(length(isolates) - (length(m_new1)==1) - (length(m_new2)==1))
          if(length(m_new1) == 1) mean.att <- mean.att + nodes[m_new1,att] 
          if(length(m_new2) == 1) mean.att <- mean.att + nodes[m_new2,att] 
        }
        statistics[e] <- mean.att/length(isolates)
      } else {
        statistics[e] <- 0
      }
      
    }
    
    # --------- ATTRIBUTE GROUPS -----------
    if(effect.name == "attgroups") {
      if(length(groups) > 0) {
        att <- which(colnames(nodes) == object.name)
        if(change == "merge"){
          m_old1 <- which(old.partition == old_g1)
          m_old2 <- which(old.partition == old_g2)
          mean.att <- old.z[e]*((num.nodes-length(isolates)) - (length(m_old1)==1) - (length(m_old2)==1))
          if(length(m_old1) == 1) mean.att <- mean.att + nodes[m_old1,att] 
          if(length(m_old2) == 1) mean.att <- mean.att + nodes[m_old2,att] 
        }else if(change == "division"){
          m_new1 <- which(new.partition == new_g1)
          m_new2 <- which(new.partition == new_g2)
          mean.att <- old.z[e]*((num.nodes-length(isolates)) + (length(m_new1)==1) + (length(m_new2)==1))
          if(length(m_new1) == 1) mean.att <- mean.att - nodes[m_new1,att] 
          if(length(m_new2) == 1) mean.att <- mean.att - nodes[m_new2,att] 
        }
        statistics[e] <- mean.att/(num.nodes-length(isolates))
      }else{
        statistics[e] <- 0
      }
      
    }
    
    # --------- HOMOPHILY:SAME -----------
    if(effect.name == "same") {
      if(length(groups) > 0) {
        att <- which(colnames(nodes) == object.name)
        if(change == "merge"){
          num.same <- old.z[e]*((num.nodes-length(isolates))+1)
          members1 <- which(old.partition == old_g1)
          members2 <- which(old.partition == old_g2)
          for(m1 in members1){
            for(m2 in members2){
              num.same <- num.same + (nodes[m1,att] == nodes[m2,att])
            }
          }
        }else if(change == "division"){
          num.same <- old.z[e]*((num.nodes-length(isolates))-1)
          members1 <- which(new.partition == new_g1)
          members2 <- which(new.partition == new_g2)
          for(m1 in members1){
            for(m2 in members2){
              num.same <- num.same - (nodes[m1,att] == nodes[m2,att])
            }
          }
        }
        statistics[e] <- num.same/(num.nodes-length(isolates))
      }else{
        statistics[e] <- 0
      }
      
    }
    
    # --------- HOMOPHILY:DIFF -----------
    if(effect.name == "diff") {
      if(length(groups) > 0) {
        att <- which(colnames(nodes) == object.name)
        if(change == "merge"){
          diff <- old.z[e]*((num.nodes-length(isolates))+1)
          members1 <- which(old.partition == old_g1)
          members2 <- which(old.partition == old_g2)
          for(m1 in members1){
            for(m2 in members2){
              diff <- diff + abs(nodes[m1,att] - nodes[m2,att])
            }
          }
        }else if(change == "division"){
          diff <- old.z[e]*((num.nodes-length(isolates))-1)
          members1 <- which(new.partition == new_g1)
          members2 <- which(new.partition == new_g2)
          for(m1 in members1){
            for(m2 in members2){
              diff <- diff - abs(nodes[m1,att] - nodes[m2,att])
            }
          }
        }
        statistics[e] <- diff/(num.nodes-length(isolates))
      } else {
        statistics[e] <- 0
      }
      
    }
    
  }
  
  return(statistics)
  
}