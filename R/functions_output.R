######################################################################
## Simulation and estimation of Exponential Random Partition Models ## 
## Functions to plot the estimation algorithm results               ##
## Author: Marion Hoffman                                           ##
######################################################################




print_results <- function(result){
  
  num.effects <- length(result$effect)
  
  effect <- result$effect
  object <- result$object
  est <- result$est
  std.err <- result$std.err
  conv <- result$conv
  t <- est / std.err
  sig <- rep("", num.effects)
  sig[abs(t) > qnorm(1 - 0.05/2)] <- "*"
  sig[abs(t) > qnorm(1 - 0.01/2)] <- "**"
  sig[abs(t) > qnorm(1 - 0.001/2)] <- "***"
  
  print( data.frame(effect, object, est, std.err, sig, t, conv) )
}

print_results_bayesian <- function(result){
  
  num.effects <- length(result$effect)
  
  effect <- result$effect
  object <- result$object
  post.mean <- result$post.mean
  post.sd <- result$post.sd
  cred.min <- post.mean - 1.95996*post.sd
  cred.max <- post.mean + 1.95996*post.sd
  
  print( data.frame(effect, object, post.mean, post.sd, cred.min, cred.max) )

  
}
