s <- list(p,p)
s[[1]][1]
p <- c(1,2,2,3,3,4,4,4,5)

p1 <- c(1,2,2,3,3,4,4,4,5)
p2 <- c(1,2,4,3,1,4,3,4,1)
p3 <- c(1,1,1,3,3,4,2,4,2)




at3 <- c(1,0,0,0,1,1,0,0,1)
at <- c(3,5,23,2,1,0,3,9,2)
at2 <- c(3,5,20,2,1,0,0,9,0)

tapply(at,p,var)
table(p)
help(shuffle)
# Average size of groups in a partition

compute_average_size <- function(partition){
  av_size <- mean(table(partition))
  return(av_size)
}
compute_average_size(p2)

# standard deviation of groups' sizes

compute_sd_size <- function(partition){
  sd_size <- sd(table(partition))
  return(sd_size)
}
compute_sd_size(p)

# Proportion of isolates

compute_prop_isolate <- function(partition){
  prop <- sum(table(partition)==1) / length(partition)
  return(prop)
}
sum(table(p)==1)/9
compute_prop_isolate(p)

# Correlation degree/ group size XXXXX

compute_correlation_degree_size <- function(partition, attribute){
#Laisser tomber finalement
}

# Intra class correlation

compute_correlation_intra_class <- function(partition, attribute){

  sum_between <- 0
  sum_within <- 0

  for (g in 1:max(partition,na.rm = F)){
    members <- which(partition == g)

    var_b <- mean(at[members]) - mean(at)
    sum_between <- sum_between + var_b^2

    var_w <- sum((at[members] - mean(at[members]))^2)
    sum_within <- sum_within + var_w
  }

  between <- sum_between / g-1
  within <- sum_within /(length(partition)-g)

  coef <- between/(between+within)

  return(coef)
}

compute_correlation_intra_class(p,at)

computeicc(at,p)
# difference average and isolates

compute_diff_average_isolate <- function(partition,attribute){

  isolates <- match(which(table(partition)==1),partition)
  diff <- mean(attribute) - mean(attribute[isolates])

  return(diff)
}

match(which(table(p)==1),p)
compute_diff_average_isolate(p,at)

# Within groups correlation 2 var

compute_correlation_two_attributes_intra <- function(partition, attribute1, attribute2, group){

  members <- which(partition == group)

  sum_att1 <- attribute1[members] - mean(attribute1[members])
  sum_att2 <- attribute2[members] - mean(attribute2[members])

  sum2_att1 <- sum((attribute1[members] - mean(attribute1[members]))^2)
  sum2_att2 <- sum((attribute2[members] - mean(attribute2[members]))^2)

  sum_att <- sum(sum_att1*sum_att2)
  coef <- sum_att / (sqrt(sum2_att1) * sqrt(sum2_att2))

  return(coef)
}
compute_correlation_two_attributes(p,at,at,4)
compute_correlation_two_attributes(p,at,at2,4)


#correlation between the group averags of the two covariates XXXXX

compute_correlation_two_attributes_inter <- function(partition, attribute1,attribute2) {

  for (g in 1:max(partition,na.rm = F)) {
    members <- which(partition == g)

  }
}


# General correlation with size

computecorrelationwithsize_gen <- function(attribute,teams, categorical){

  n <- length(attribute)
  m <- max(teams)

  sizes <- 0

  for(i in 1:n){
    t <- teams[i]
    sizes[i] <- length(which(teams==t))
  }

  if (categorical == F){
    c <- cor(attribute,sizes)
    return(c)
    } else{

    model <- summary(lm(sizes ~ attribute))
    return(cor.ratio = model$r.squared)
  }
 
}
# p.value = model$coefficients[,4][2])
computecorrelationwithsize_gen(at3,p,T)
computecorrelationwithsize_gen(at,p,F)


# General compute proportion
# prop par group a ajouter 

computeproportion_gen <- function(attribute, teams, class){

  sum_g <- 0
  x <- 0

  if (class == "all"){
    attribute <- as.factor(attribute)

    for (level in levels(attribute)){
      x[level]<-0
    }

    for(g in 1:max(teams, na.rm = T)){
      members <- which(teams == g)
      table_g <- table(attribute[members])
      prop_g <- prop.table(table_g)
      for (level in dimnames(table_g)){
        x[level] <- x[level] + prop_g[level]
      }

    }
    return(x/max(teams, na.rm = T))

  }

  else{
    bin <- ifelse(attribute==class, 1,0)
    for(g in 1:max(teams, na.rm = T)){
      members <- which(teams == g)
      sum_g <- sum_g + (sum(bin[members] == 1)/length(members))
    }
    return(sum_g/max(teams, na.rm = T))

  }

}


computeproportion_gen(at3,p,'all')
computeproportion_gen(at3,p,1)



# Mettre les types en arg ( de cot )


# Pareil que compute _ proportion 
compute_number_ind_binary_withtypes <- function(partition_obs, types_obs, attribute, type) {
  
  partition_obs2 <- partition_obs
  out_groups <- which(types_obs != type)
  partition_obs2[partition_obs2 %in% out_groups] <- NA
  partition_obs2[!is.na(partition_obs2)] <- order_groupids(partition_obs2[!is.na(partition_obs2)])
  
  numinds_obs <- sum(attribute * !is.na(partition_obs2))
  avgroupnuminds_obs <- sum(attribute * !is.na(partition_obs2))  / (max(partition_obs2, na.rm = T))

  res <- list(number_individuals= numinds_obs,
              averagepergroup_number_individuals = avgroupnuminds_obs)
  
  return(res)
}

#### SAME PAIRS ####

compute_number_same_pairs <- function(partition, attribute) {
  
  dist<- as.matrix(dist(partition, diag = T, upper = T))
  net_partition <- dist
  net_partition[dist == 0] <- 1
  net_partition[dist != 0] <- 0
  net_partition[is.na(dist)] <- 0
  diag(net_partition) <- 0
  
  dist <- as.matrix(dist(as.numeric(factor(attribute)), diag = T, upper = T))
  distances <- dist
  distances[dist == 0] <- 1
  distances[dist != 0] <- 0
  diag(distances) <- 0
  
  samepairs <- sum(distances * net_partition) / 2
  
  res <- samepairs
  
  
  return(res)
}

compute_ind_same_pairs <- function(partition, attribute) {
  
  dist<- as.matrix(dist(partition, diag = T, upper = T))
  net_partition <- dist
  net_partition[dist == 0] <- 1
  net_partition[dist != 0] <- 0
  net_partition[is.na(dist)] <- 0
  diag(net_partition) <- 0
  
  dist <- as.matrix(dist(as.numeric(factor(attribute)), diag = T, upper = T))
  distances <- dist
  distances[dist == 0] <- 1
  distances[dist != 0] <- 0
  diag(distances) <- 0
  
  indsamepairs <- sum( rowSums(distances * partition) > 0 )
  
  res <- indsamepairs
  
  
  return(res)
}

compute_avg_same_pairs <- function(partition, attribute) {
  
  dist<- as.matrix(dist(partition, diag = T, upper = T))
  net_partition <- dist
  net_partition[dist == 0] <- 1
  net_partition[dist != 0] <- 0
  net_partition[is.na(dist)] <- 0
  diag(net_partition) <- 0
  
  dist <- as.matrix(dist(as.numeric(factor(attribute)), diag = T, upper = T))
  distances <- dist
  distances[dist == 0] <- 1
  distances[dist != 0] <- 0
  diag(distances) <- 0

  avgroupsameties <- sum(distances * net_partition) / (2*max(partition, na.rm = T))
  
  res <- avgroupsameties
  
  return(res)
}



#### TIES #### 

# For dyadic attributes

# Mettre le type plus tard

compute_number_ties <- function(partition, net_attribute) {
  
  dist<- as.matrix(dist(partition, diag = T, upper = T))
  net_partition <- dist
  net_partition[dist == 0] <- 1
  net_partition[dist != 0] <- 0
  net_partition[is.na(dist)] <- 0
  diag(net_partition) <- 0
  
  numties <- sum(net_attribute * net_partition) / 2
  
  res <- numties
  
  return(res)
}

compute_avg_number_ties <- function(partition, net_attribute) {
  
  dist<- as.matrix(dist(partition, diag = T, upper = T))
  net_partition <- dist
  net_partition[dist == 0] <- 1
  net_partition[dist != 0] <- 0
  net_partition[is.na(dist)] <- 0
  diag(net_partition) <- 0
  
  avgroupnumties <- sum(net_attribute * net_partition) / (2*max(partition, na.rm = T))

  res <- avgroupnumties
  
  return(res)
}


compute_ind_ties_dyadic <- function(partition, net_attribute) {
  
  dist<- as.matrix(dist(partition, diag = T, upper = T))
  net_partition <- dist
  net_partition[dist == 0] <- 1
  net_partition[dist != 0] <- 0
  net_partition[is.na(dist)] <- 0
  diag(net_partition) <- 0

  indnumties <- sum( rowSums(net_attribute * net_partition) > 0 )
  
  res <- indnumties
  
  return(res)
}


# mettre de cote
compute_number_ties_dyadicwithtypes <- function(partition_obs, types_obs, net_attribute) {
  
  partition_obs2 <- partition_obs
  out_groups <- which(types_obs != type)
  partition_obs2[partition_obs2 %in% out_groups] <- NA
  partition_obs2[!is.na(partition_obs2)] <- order_groupids(partition_obs2[!is.na(partition_obs2)])
  
  dist<- as.matrix(dist(partition_obs2, diag = T, upper = T))
  net_partition_obs2 <- dist
  net_partition_obs2[dist == 0] <- 1
  net_partition_obs2[dist != 0] <- 0
  net_partition_obs2[is.na(dist)] <- 0
  diag(net_partition_obs2) <- 0

  numties_obs <- sum(net_attribute * net_partition_obs2) / 2
  avgroupnumties_obs <- sum(net_attribute * net_partition_obs2) / (2*max(partition_obs2, na.rm = T))
  indnumties_obs <- sum( rowSums(net_attribute * net_partition_obs2) > 0 )
  
  res <- list(number_ties = numties_obs,
              averagepergroup_number_ties = avgroupnumties_obs,
              individual_number_ties = indnumties_obs)
  
  return(res)
}

#### SIMILAR PAIRS ####

# Continuous attribute 

compute_number_similar_pairs <- function(partition, attribute, threshold) {
  
  dist<- as.matrix(dist(partition, diag = T, upper = T))
  net_partition <- dist
  net_partition[dist == 0] <- 1
  net_partition[dist != 0] <- 0
  net_partition[is.na(dist)] <- 0
  diag(net_partition) <- 0
  
  dist <- as.matrix(dist(attribute, diag = T, upper = T))
  distances <- dist
  distances[dist <= threshold] <- 1
  distances[dist > threshold] <- 0
  diag(distances) <- 0
  
  sameties <- sum(distances * net_partition) / 2
  
  res <- sameties
  
  return(res)
}

compute_avg_similar_pairs <- function(partition, attribute, threshold) {
  
  dist<- as.matrix(dist(partition, diag = T, upper = T))
  net_partition <- dist
  net_partition[dist == 0] <- 1
  net_partition[dist != 0] <- 0
  net_partition[is.na(dist)] <- 0
  diag(net_partition) <- 0
  
  dist <- as.matrix(dist(attribute, diag = T, upper = T))
  distances <- dist
  distances[dist <= threshold] <- 1
  distances[dist > threshold] <- 0
  diag(distances) <- 0

  avgroupsameties <- sum(distances * net_partition) / (2*max(partition, na.rm = T))
  
  res <- avgroupsameties
  
  return(res)
}


compute_ind_similar_pairs <- function(partition, attribute, threshold) {
  
  dist<- as.matrix(dist(partition, diag = T, upper = T))
  net_partition <- dist
  net_partition[dist == 0] <- 1
  net_partition[dist != 0] <- 0
  net_partition[is.na(dist)] <- 0
  diag(net_partition) <- 0
  
  dist <- as.matrix(dist(attribute, diag = T, upper = T))
  distances <- dist
  distances[dist <= threshold] <- 1
  distances[dist > threshold] <- 0
  diag(distances) <- 0

  indsameties <- sum( rowSums(distances * net_partition) > 0 )
  
  res <- indsameties
  
  return(res)
}


# En construction

#### Fontion Globale CUG ####
 
# Pour l'instant avec juste un sample contenant des partitions simplement 
# l'element partition_sample est un vecteur de vecteur

#Pour tester  faire : CUG(obs, fun = function(p){compute...(p,..,..,..)}, part_sample)


CUG <- function(observation, fun, partition_sample) {
  
  if (is.null(partition_sample)) {
      partition_sample = list()
      for (i in 1:10)
        partition_sample[[i]] <- sample(observation, replace = T)
  }
  
  stat_obs <- fun(observation)
  nb_partions <- length(partition_sample)
  
  stat_sample <- 0:nb_partions
  for(i in 1:nb_partions) {
      stat_sample[i] <- fun(partition_sample[[i]])
  }

  minnum <- quantile(stat_sample, probs = 0.05)
  maxnum <- quantile(stat_sample, probs = 0.95)
 
  p1 <- sum(stat_sample <= stat_obs) / nb_partions
  p2 <- sum(stat_sample >= stat_obs) / nb_partions
  if(p1 == 0) {
    p_num <- p2
  } else if(p2 == 0) {
    p_num <- p1
  } else {
    p_num <- 2*min(p1,p2)
  }

  res <- list(observed = stat_obs,
              mean_sample = mean(stat_sample),
              sd_sample = sd(stat_sample),
              min_CI_95 = minnum,
              max_CI_95 = maxnum,
              p.value = p_num)
  
  return(res)
}

#Works
sample <- list(p1,p2,p3)
CUG(p,fun=function(part){compute_number_same_pairs(part,at3)}, sample)

x = list()
x[[1]] <- sample(p, replace = T)
