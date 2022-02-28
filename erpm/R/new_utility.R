######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Utility functions and CUG                                        ##
######################################################################


##### CORRELATION FUNCTIONS #####


#' Intra class correlation
#'
#'This function computes the intra class correlation :correlation
#'of attributes for 2 randomly drawn individuals in the same group.
#'
#' @param partition A partition (vector)
#' @param attribute A vector containing the values of the attribute
#' @return A number corresponding to the ICC
#' @examples
#' p <- c(1,2,2,3,3,4,4,4,5)
#' at <- c(3,5,23,2,1,0,3,9,2)
#' icc(p, at)


icc <- function(partition, attribute){

  sum_between <- 0
  sum_within <- 0

  for (g in 1:max(partition,na.rm = F)){
    members <- which(partition == g)

    var_b <- mean(attribute[members]) - mean(attribute)
    sum_between <- sum_between + var_b^2

    var_w <- sum((attribute[members] - mean(attribute[members]))^2)
    sum_within <- sum_within + var_w
  }

  between <- sum_between / g-1
  within <- sum_within /(length(partition)-g)

  coef <- between/(between+within)

  return(coef)
}



#' Within groups correlation
#'
#'This function computes the correlation between the two attributes for individuals in the same group.
#'
#' @param partition A partition (vector)
#' @param attribute1 A vector containing the values of the first attribute
#' @param attribute2 A vector containing the values of the second attribute
#' @param group A number indicating the selected group
#' @return A number corresponding to the correlation coefficient
#' @examples
#' p <- c(1,2,2,3,3,4,4,4,5)
#' at <- c(3,5,23,2,1,0,3,9,2)
#' at2 <- c(3,5,20,2,1,0,0,9,0)
#' correlation_within(p,at,at2,4)


correlation_within <- function(partition, attribute1, attribute2, group){

  members <- which(partition == group)

  sum_att1 <- attribute1[members] - mean(attribute1[members])
  sum_att2 <- attribute2[members] - mean(attribute2[members])

  sum2_att1 <- sum((attribute1[members] - mean(attribute1[members]))^2)
  sum2_att2 <- sum((attribute2[members] - mean(attribute2[members]))^2)

  sum_att <- sum(sum_att1*sum_att2)
  coef <- sum_att / (sqrt(sum2_att1) * sqrt(sum2_att2))

  return(coef)
}



#' Between groups correlation
#'
#'This function computes the correlation between the group averages of the two attributes.
#'
#' @param partition A partition (vector)
#' @param attribute1 A vector containing the values of the first attribute
#' @param attribute2 A vector containing the values of the second attribute
#' @return A number corresponding to the correlation coefficient
#' @examples
#' p <- c(1,2,2,3,3,4,4,4,5)
#' at <- c(3,5,23,2,1,0,3,9,2)
#' at2 <- c(3,5,20,2,1,0,0,9,0)
#' correlation_between(p,at,at2)


correlation_between <- function(partition, attribute1, attribute2) {

  mean1 <- mean(attribute1)
  mean2 <- mean(attribute2)
  sum_mean <- 0
  sum1 <- 0
  sum2 <- 0
  for (g in 1:max(partition,na.rm = F)) {
    members <- which(partition == g)
    sum_mean <- sum_mean + (mean(attribute1[members])-mean1)*(mean(attribute2[members])-mean2)
    sum1 <- sum1 + (mean(attribute1[members])-mean1)^2
    sum2 <- sum2 + (mean(attribute2[members])-mean2)^2
  }
  res <- sum_mean/(sqrt(sum1)*sqrt(sum2))
  return(res)
}



#' Correlation with size
#'
#'This function computes the correlation between an attribute and the size of the groups.
#'
#' @param partition A partition (vector)
#' @param attribute A vector containing the values of the attribute
#' @param categorical Boolean (True or False) indicating if the attribute is categorical
#' @return A number corresponding to the correlation coefficient if the attribute is numerical or
#' the correlation ratio if the attribute is categorical.
#' @examples
#' p <- c(1,2,2,3,3,4,4,4,5)
#' at <- c(3,5,23,2,1,0,3,9,2)
#' at2 <- c(3,5,20,2,1,0,0,9,0)
#' correlation_with_size(p,at,F)


correlation_with_size <- function(partition, attribute, categorical){

  n <- length(attribute)
  m <- max(partition)

  sizes <- 0

  for(i in 1:n){
    t <- partition[i]
    sizes[i] <- length(which(partition==t))
  }

  if (categorical == F){
    c <- cor(attribute,sizes)
    return(c)
    } else{

    model <- summary(lm(sizes ~ attribute))
    return(cor.ratio = model$r.squared)
  }

}


## UTILITY FUNCTIONS ------


#' Statistics on the size of groups in a partition
#'
#'This function computes the average or the standard
#'deviation of the size of groups in a partition.
#'
#' @param partition A partition (vector)
#' @param stat The statistic to compute : 'avg' for average and 'sd' for standard deviation
#' @return A number corresponding to the correlation coefficient if the attribute is numerical or
#' the correlation ratio if the attribute is categorical.
#' @examples
#' p <- c(1,2,2,3,3,4,4,4,5)
#' correlation_with_size(p,'avg')
#' correlation_with_size(p,'sd')

stat_size <- function(partition, stat){

  if (stat == 'avg'){
  av_size <- mean(table(partition))
  return(av_size)
  }
  if (stat == 'sd'){
  sd_size <- sd(table(partition))
  return(sd_size)
  }
}




#' Proportion of isolates
#'
#'This function computes the proportion of individuals not joining others.
#'
#' @param partition A partition (vector)
#' @return A number corresponding to proportion of individuals alone.
#' @examples
#' p <- c(1,2,2,3,3,4,4,4,5)
#' prop_isolate(p)

prop_isolate <- function(partition){
  prop <- sum(table(partition)==1) / length(partition)
  return(prop)
}




#' Number of individuals having an attribute
#'
#'This function computes the total number of individuals being in a category of an attribute in a partition.
#'It also computes the sum of the proportion in each group of individuals being in
#'a category.
#'
#' @param partition A partition (vector)
#' @param attribute A vector containing the values of the attribute
#' @param stat The statistic to compute : 'avg' for the sum of proportion and 'sum' for the total number
#' @param category The category to consider or category = 'all' if all categories have to be considered
#' @return The statisic chosen in stat depending on the value of category. If category = 'all', returns a vector.
#' @examples
#' p <- c(1,2,2,3,3,4,4,4,5)
#' at <- c(1,0,0,0,1,1,0,0,1)
#' number_categories(p,at,sum,'all')


number_categories <- function(partition, attribute, stat, category){

  sum_g <- 0
  x <- 0

  if (category == "all"){
    attribute <- as.factor(attribute)
    for (level in levels(attribute)){
      x[level]<-0
    }


    if (stat == 'sum'){
      for(g in 1:max(partition, na.rm = T)){
        members <- which(partition == g)
        table_g <- table(attribute[members])
        for (level in dimnames(table_g)){
          x[level] <- x[level] + table_g[level]
        }
      }
      return(x)
    }


    if (stat == 'avg'){
      for(g in 1:max(partition, na.rm = T)){
        members <- which(partition == g)
        table_g <- table(attribute[members])
        prop_g <- prop.table(table_g)
        for (level in dimnames(table_g)){
          x[level] <- x[level] + prop_g[level]
        }
      }
      return(x[-1])
    }

  }

  else {
    bin <- ifelse(attribute==category, 1,0)

    if (stat == 'sum'){
      for(g in 1:max(partition, na.rm = T)){
        members <- which(partition == g)
        sum_g <- sum_g + sum(bin[members] == 1)
      }
      return(sum_g)
    }

    if (stat == 'avg'){
      for(g in 1:max(partition, na.rm = T)){
        members <- which(partition == g)
        sum_g <- sum_g + (sum(bin[members] == 1)/length(members))
      }
      return(sum_g)
    }

  }

}



#' Density of intra ties
#'
#'This function computes the average "density" of dyadic attributes within groups and the
#'dispersion of the "density" of dyadic attributes within groups.
#'
#' @param partition A partition (vector)
#' @param matrix A matrix containing the values of the dyadic attribute
#' @param stat The statistic to compute : 'avg' for the average and 'sd' for the standard deviation
#' @return The statisic chosen in stat.
#' @examples
#' p <- c(1,2,2,3,3,4,4,4,5)
#' mat <- matrix(1:81, ncol = 9)
#' density(p,mat,'avg')


density <- function(partition, matrix, stat){

  m <- max(partition,na.rm=T)
  avd <- 0
  std <- 0

  for(h in 1:m){
    members <- which(partition==h)

    if(length(members)>1){
      beforeties <- 0
      allties <- length(members)*(length(members)-1)/2

      for(i in 1:(length(members)-1)){
        for(j in (i+1):length(members)){
          beforeties <- beforeties + matrix[members[i],members[j]]
        }
      }

      avd <- avd + beforeties/allties
      std <- std + (beforeties/allties)^2
    }
  }

  avd <- avd/m
  std <- sqrt(1/m*std - avd^2)
  if (stat == 'avg'){
   return(average=avd)
  }
  if (stat == 'sd'){
    return(standard.deviation=std)
  }
}


#' Range of attribute in groups
#'
#'This function computes the sum or the average range of an attribute for groups in a partition.
#'
#' @param partition A partition (vector)
#' @param attribute A vector containing the values of the attribute
#' @param stat The statistic to compute : 'avg' for the average and 'sum' for the sum
#' @return The statisic chosen in stat
#' @examples
#' p <- c(1,2,2,3,3,4,4,4,5)
#' at <- c(3,5,23,2,1,0,3,9,2)
#' density(p,at,'avg')

range <- function(partition, attribute, stat ){

  sum <- 0
  for(g in 1:max(partition, na.rm = T)){
    members <- which(partition == g)
    sum <- sum + max(attribute[members], na.rm = T) - min(attribute[members], na.rm = T)
  }

  if (stat == 'sum'){
    return(sum = sum )
  }
  if (stat == 'avg'){
    return(average = sum/max(partition, na.rm = T))
  }
}



#' Same pairs of individuals in a partition
#'
#'This function computes the total number, the average number having the same value of a categorical variable
#'and the number of individuals a partition.
#'
#' @param partition A partition (vector)
#' @param attribute A vector containing the values of the attribute
#' @param stat The statistic to compute : 'avg' for the average, 'sum' for the sum,  'ind' for the number of ties
#' each individual has in its group.
#' @return The statisic chosen in stat
#' @examples
#' p <- c(1,2,2,3,3,4,4,4,5)
#' at <- c(0,1,1,1,1,0,0,0,0)
#' same_pairs(p,at,'avg')


same_pairs <- function(partition, attribute, stat ) {

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
  indsamepairs <- sum( rowSums(distances * net_partition) > 0 )
  avgroupsamepairs <- sum(distances * net_partition) / (2*max(partition, na.rm = T))

  if (stat == 'sum'){
    return(samepairs)
  }
  if (stat == 'avg'){
    return(avgroupsamepairs)
  }

  if (stat == 'ind'){
    return(indsamepairs)
  }

}



#' Same pairs of individuals in a partition
#'
#'This function computes the number of ties.
#'
#' @param partition A partition (vector)
#' @param attribute A vector containing the values of the attribute
#' @param stat The statistic to compute : 'avg' for the average, 'sum' for the sum, 'ind' for the number of ties
#' each individual has in its group.
#' @return The statisic chosen in stat
#' @examples
#' p <- c(1,2,2,3,3,4,4,4,5)
#' at <- c(0,1,1,1,1,0,0,0,0)
#' number_ties(p,at,'avg')

number_ties <- function(partition, dyadic_attribute) {

  net_partition <- dyadic_attribute
  net_partition[dyadic_attribute == 0] <- 1
  net_partition[dyadic_attribute != 0] <- 0
  net_partition[is.na(dyadic_attribute)] <- 0
  diag(net_partition) <- 0

  numties <- sum(attribute * net_partition) / 2
  avgroupnumties <- sum(attribute * net_partition) / (2*max(net_partition, na.rm = T))
  indnumties <- sum( rowSums(attribute * net_partition) > 0 )


  if (stat == 'sum'){
    return(numties)
  }
  if (stat == 'avg'){
    return(avgroupnumties)
  }

  if (stat == 'ind'){
    return(indnumties)
  }

}


#' Similar pairs of individuals in a partition
#'
#'This function computes the total number, the average number having the close values of a numerical variable
#'and the number of individuals a partition.
#'
#' @param partition A partition (vector)
#' @param attribute A vector containing the values of the attribute
#' @param stat The statistic to compute : 'avg' for the average, 'sum' for the sum, 'ind' for individuals
#' @param threshold Threshold to determine if 2 individuals attributes values are close
#' @return The statisic chosen in stat
#' @examples
#' p <- c(1,2,2,3,3,4,4,4,5)
#' at <- c(3,5,23,2,1,0,3,9,2)
#' similar_pairs(p,at,1,'avg')

similar_pairs <- function(partition, attribute, stat,  threshold) {

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
  avgroupsameties <- sum(distances * net_partition) / (2*max(partition, na.rm = T))
  indsameties <- sum( rowSums(distances * net_partition) > 0 )

  if (stat == 'sum'){
    return(sameties)
  }
  if (stat == 'avg'){
    return(avgroupsameties)
  }

  if (stat == 'ind'){
    return(indsameties)
  }
}



#### Fontion Globale CUG ####


#' CUG
#'
#'This function computes ??
#'
#' @param observation The observed partition
#' @param fun An utility function or statistic to compute
#' @param permutations A predicted partition or partition to compare to observation. NULL by default
#' @param num.permutations Number of permutations to do to generate a partition if permutation = NULL. 1000 by default
#' @return The value of mean, standard deviation, p-value, and confidence interval of the estimator of the statistic in fun
#' @examples
#' p <- c(1,2,2,3,3,4,4,4,5)
#'CUG(p,fun=function(part){compute_number_same_pairs(part,at3,'avg')})


CUG <- function(observation, fun, permutations=NULL, num.permutations=1000) {

  if (is.null(permutations)) {
    permutations = list()
      for (i in 1:num.permutations)
        permutations[[i]] <- order_groupids(sample(observation, replace = F))
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
