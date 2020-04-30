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
  
  # TODO: numerical estimation of the log likelihood?
  #aic <- -2*result$log.likelihood + 2 * num.effects
  # TODO: how many observations are there??
  #aicc <- aic + 2*num.effects*(num.effects + 1)/(result$n.events - num.effects - 1)
  #bic <- -2 * result$log.likelihood + num.effects * log(result$n.events)
  
  print( data.frame(effect, object, est, std.err, sig, t, conv) )
  #cat(" ", paste("Log likelihood", round(result$log.likelihood, 4), "\n"))
  #cat(" ", paste(ifelse(max(abs(result$conv))>0.1, "Converged", "Not converged"), "with max abs. score of", 
  #               round(max(abs(result$scores)), 5)), "\n")
  #cat(" ", paste("AIC ", round(AIC(result), 5)  )) 
                 #"\n  AICc", round(aicc, 5), 
                 #"\n  BIC ", round(BIC(result), 5)), "\n")
  
}
