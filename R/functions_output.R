######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Functions to plot the estimation algorithm results               ##
## Author: Marion Hoffman                                           ##
######################################################################



#' Print estimation results
#'
#'
#' @param x output of the estimate function
#' @param ... For internal use only.
#'
#' @return a data frame
#' @importFrom stats qnorm
#' @export
print.results.list.erpm <- function(x, ...) {

  result <- x$results

  num.effects <- length(result$effect)

  effect <- result$effect
  object <- result$object
  est <- result$est
  std.err <- result$std.err
  conv <- result$conv
  t <- est / std.err
  sig <- rep("", num.effects)
  sig[abs(t) > qnorm(1 - 0.05 * 0.5)] <- "*"
  sig[abs(t) > qnorm(1 - 0.01 * 0.5)] <- "**"
  sig[abs(t) > qnorm(1 - 0.001 * 0.5)] <- "***"

  print(data.frame(effect, object, est, std.err, sig, t, conv))

}


#' Print results of estimation of phase 3
#'
#'
#' @param x output of the estimate function
#' @param ... For internal use only.
#'
#' @return a data frame
#' @importFrom stats qnorm
#' @export
print.results.p3.erpm <- function(x, ...) {

  result <- x

  num.effects <- length(result$effect)

  effect <- result$effect
  object <- result$object
  est <- result$est
  std.err <- result$std.err
  conv <- result$conv
  t <- est / std.err
  sig <- rep("", num.effects)
  sig[abs(t) > qnorm(1 - 0.05 * 0.5)] <- "*"
  sig[abs(t) > qnorm(1 - 0.01 * 0.5)] <- "**"
  sig[abs(t) > qnorm(1 - 0.001 * 0.5)] <- "***"

  print(data.frame(effect, object, est, std.err, sig, t, conv))
}


#' Print results of bayesian estimation (beta version)
#'
#' @param x output of the bayesian estimate function
#' @param ... For internal use only.
#'
#' @return a data frame
#' @export
print.results.bayesian.erpm <- function(x, ...) {

  result <- x$results

  num.effects <- length(result$effect)

  effect <- result$effect
  object <- result$object
  post.mean <- result$post.mean
  post.sd <- result$post.sd
  cred.min <- post.mean - 1.95996 * post.sd
  cred.max <- post.mean + 1.95996 * post.sd

  print(data.frame(effect, object, post.mean, post.sd, cred.min, cred.max))


}
