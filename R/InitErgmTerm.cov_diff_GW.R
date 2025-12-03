#' ERGM term: cov_diff_GW (geometrically weighted range over k-actor subsets)
#' @name InitErgmTerm.cov_diff_GW
#' @aliases cov_diff_GW
#' @note InitErgmTerm.cov_diff_GW.R
#'
#' @description
#' \code{cov_diff_GW} is an ERGM term for bipartite networks that builds a
#' geometrically weighted combination of \code{cov_diff}-type statistics over
#' all subset sizes \eqn{k \ge 2} inside each group. The bipartite network is
#' interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side, representing groups).
#' }
#'
#' Each actor in the actor mode carries a numeric covariate value \eqn{x_i}. For
#' each group and each subset size \eqn{k \ge 2}, \code{cov_diff_GW} considers
#' all \eqn{k}-actor subsets inside that group and computes on each subset the
#' max–min range of the covariate. These \code{cov_diff}-type contributions are
#' then combined over all \eqn{k \ge 2} using a geometric weight controlled by
#' \eqn{\lambda > 1}.
#'
#' For a given value of \eqn{\lambda > 1}, let
#' \eqn{c_k} denote the \code{cov_diff}-style statistic of order \eqn{k}, i.e.
#' the sum of ranges over all \eqn{k}-actor subsets in all groups. The term
#' \code{cov_diff_GW} is defined as
#' \deqn{
#'   T_{\mathrm{GW}}(\lambda)
#'   =
#'   \sum_{k \ge 2} \left(-\frac{1}{\lambda}\right)^{k-1} c_k.
#' }
#' The initializer supports a vector of \eqn{\lambda} values and returns one
#' scalar statistic per value.
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic, declared in the
#' compiled code under a symbol compatible with \code{name = "cov_diff_GW"}.
#' The R initializer below:
#' \itemize{
#'   \item enforces that the network is bipartite via \code{nw \%n\% "bipartite"};
#'   \item extracts an actor-level covariate from a vertex attribute or from a
#'         literal vector, coercing it to numeric and failing fast on \code{NA};
#'   \item normalizes the \code{lambda} argument into a numeric vector with all
#'         entries strictly greater than 1;
#'   \item packs the actor covariate and the vector of \eqn{\lambda} values into
#'         a compact \code{INPUT_PARAM} layout for the C layer.
#' }
#'
#' The ERGM infrastructure will call the C change-statistic whenever a toggle
#' affects an edge between an actor and a group. The C code recomputes the
#' local contribution of the affected group(s) to each geometrically weighted
#' statistic \eqn{T_{\mathrm{GW}}(\lambda_\ell)} for \eqn{\ell = 1,\dots,L}.
#'
#' Internally, \code{cov_diff_GW} passes its parameters and covariate to the C
#' layer through \code{INPUT_PARAM} with layout
#' \deqn{
#'   \text{INPUT\_PARAM}
#'   =
#'   \bigl(
#'     n_1,\,
#'     L,\,
#'     \lambda_1,\dots,\lambda_L,\,
#'     x_1,\dots,x_{n_1}
#'   \bigr),
#' }
#' where:
#' \itemize{
#'   \item \eqn{n_1} is the actor-mode size (number of actors);
#'   \item \eqn{L} is the number of \eqn{\lambda} values;
#'   \item \eqn{\lambda_1,\dots,\lambda_L} are the strictly positive geometric
#'         parameters (each greater than 1);
#'   \item \eqn{x_1,\dots,x_{n_1}} is the numeric covariate on the actor mode.
#' }
#'
#' @section Mathematical definition:
#' Let:
#' \itemize{
#'   \item \eqn{A} be the set of actor-mode nodes;
#'   \item \eqn{G} be the set of group-mode nodes;
#'   \item \eqn{B} the bipartite adjacency between actors and groups;
#'   \item \eqn{x_i \in \mathbb{R}} the numeric covariate value of actor
#'         \eqn{i \in A};
#'   \item \eqn{A_g = \{ i \in A : B_{ig}=1\}} the set of actors in group
#'         \eqn{g \in G};
#'   \item \eqn{n_g = |A_g|} the size of group \eqn{g};
#'   \item \eqn{\mathcal{C}_k(g)} the family of all \eqn{k}-element subsets of
#'         \eqn{A_g}.
#' }
#'
#' For each \eqn{k \ge 2} and each subset \eqn{S \in \mathcal{C}_k(g)}, define
#' the range
#' \deqn{
#'   R(S) = \max_{i \in S} x_i - \min_{i \in S} x_i.
#' }
#' For each \eqn{k \ge 2}, define
#' \deqn{
#'   c_k(B;x)
#'   =
#'   \sum_{g \in G}
#'   \sum_{S \in \mathcal{C}_k(g)} R(S),
#' }
#' with the convention that groups with \eqn{n_g < k} contribute zero. Then for
#' each \eqn{\lambda > 1}, the \code{cov_diff_GW} statistic is
#' \deqn{
#'   T_{\mathrm{GW}}(B;x;\lambda)
#'   =
#'   \sum_{k \ge 2}
#'     \left(-\frac{1}{\lambda}\right)^{k-1} c_k(B;x).
#' }
#' When a vector \eqn{\lambda_1,\dots,\lambda_L} is provided, the ERGM term
#' returns the vector
#' \deqn{
#'   \bigl(
#'     T_{\mathrm{GW}}(B;x;\lambda_1),\dots,
#'     T_{\mathrm{GW}}(B;x;\lambda_L)
#'   \bigr).
#' }
#'
#' @section Usage:
#' Typical usage with {ergm} on a bipartite network \code{nw}:
#' \preformatted{
#'   # Single lambda (default = 2)
#'   summary(nw ~ cov_diff_GW(cov = "x_attr", lambda = 2))
#'
#'   # Multiple lambda values
#'   summary(nw ~ cov_diff_GW(cov = "x_attr", lambda = c(1.5, 2, 3)))
#' }
#'
#' When using the ERPM wrapper on a partition-based workflow, the term can be
#' invoked indirectly as:
#' \preformatted{
#'   erpm(partition ~ cov_diff_GW(cov = "x_attr", lambda = 2))
#' }
#' provided that the wrapper translates the partition into a bipartite network
#' with a consistent actor mode and group mode.
#'
#' @note
#' The network must be strictly bipartite:
#' \itemize{
#'   \item the actor mode size is given by \code{nw \%n\% "bipartite"} and must be
#'         a positive integer;
#'   \item the group mode is the complementary set of nodes;
#'   \item the covariate is read on the actor mode only.
#' }
#'
#' The covariate is coerced to numeric if possible (logical and integer types
#' are accepted), and any \code{NA} on the actor mode is rejected in a
#' fail-fast manner. All \code{lambda} values must be finite and strictly
#' greater than 1. The initializer is vectorized in \code{lambda}: a single
#' term in the model can return multiple statistics, one per \code{lambda}
#' value.
#'
#' Debugging output for the initializer can be enabled via:
#' \preformatted{
#'   options(ERPM.cov_diff_GW.debug = TRUE)
#' }
#' which prints internal sizes and parameter summaries.
#'
#' @examples
#' \dontrun{
#'   library(network)
#'   library(ergm)
#'
#'   # Build a small bipartite network: 4 actors, 2 groups
#'   n_actors <- 4
#'   n_groups <- 2
#'   n_total  <- n_actors + n_groups
#'
#'   # Adjacency: actors 1..4, groups 5..6
#'   adj <- matrix(0, n_total, n_total)
#'
#'   # Group 5 has actors 1,2,3
#'   adj[1, 5] <- adj[5, 1] <- 1
#'   adj[2, 5] <- adj[5, 2] <- 1
#'   adj[3, 5] <- adj[5, 3] <- 1
#'
#'   # Group 6 has actors 2,4
#'   adj[2, 6] <- adj[6, 2] <- 1
#'   adj[4, 6] <- adj[6, 4] <- 1
#'
#'   nw <- network(adj, directed = FALSE, matrix.type = "adjacency")
#'   nw \%n\% "bipartite" <- n_actors  # actor mode size
#'
#'   # Numeric covariate on the actor mode
#'   x_attr <- c(1.0, 3.0, 5.0, 10.0)
#'   set.vertex.attribute(nw, "x_attr", c(x_attr, rep(NA_real_, n_groups)))
#'
#'   # Single lambda
#'   summary(nw ~ cov_diff_GW(cov = "x_attr", lambda = 2))
#'
#'   # Multiple lambda values, producing a vector of statistics
#'   summary(nw ~ cov_diff_GW(cov = "x_attr", lambda = c(1.5, 2, 3)))
#'
#'   # Simple ERGM including the term
#'   fit <- ergm(nw ~ cov_diff_GW(cov = "x_attr", lambda = 2))
#'   summary(fit)
#' }
#'
#' @section Tests:
#' Self-tests for \code{cov_diff_GW} construct small bipartite networks with
#' known groups and numeric covariates and compare, for each \code{lambda},
#' \itemize{
#'   \item the ERGM summary
#'         \code{summary(nw ~ cov_diff_GW(cov = ..., lambda = lambda_vec))};
#'   \item a direct evaluation of
#'         \eqn{T_{\mathrm{GW}}(B;x;\lambda)} via
#'         \eqn{\sum_{k \ge 2} (-1/\lambda)^{k-1} c_k(B;x)},
#'         where each \eqn{c_k} is computed explicitly by enumerating all
#'         \eqn{k}-actor subsets inside each group.
#' }
#' Additional checks verify that toggling an actor–group tie changes the
#' statistic by the local increment predicted by the C change-statistic, and
#' that all \code{lambda} values are handled consistently in the vectorized
#' interface.
#'
#' @keywords ERGM term bipartite groups covariate range geometric weighting
#' @md
#'
#' @export
InitErgmTerm.cov_diff_GW <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_diff_GW"

  # Debug helpers:
  # - dbg: logical flag, controlled by an R option;
  # - dbgcat(): emits prefixed debug messages when dbg is TRUE.
  dbg    <- isTRUE(getOption("ERPM.cov_diff_GW.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[cov_diff_GW][DEBUG]", ..., "\n", sep = "")

  # Run standard ERGM term checks and parse user arguments:
  # - enforce bipartite network;
  # - accept 'cov' and 'lambda' with flexible types;
  # - let {ergm} handle generic validations (missing args, etc.).
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov",                              "lambda"),
    vartypes      = c("character,numeric,logical,vector", "numeric,vector"),
    defaultvalues = list(NULL,                             2),
    required      = c(TRUE,                                FALSE)
  )

  # ----- 1) Actor-mode size n1 -----------------------------------------------
  # n1 is the number of actors, retrieved from the bipartite network attribute.
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L) stop(termname, ": réseau non biparti strict.")
  dbgcat("n1 = ", n1)

  # ----- 2) Extract actor covariate (length >= n1) ---------------------------
  # The covariate may be:
  # - the name of a vertex attribute (preferred);
  # - a literal vector. In both cases we keep only the first n1 entries for actors.
  cov_raw <- a$cov
  if (is.character(cov_raw) && length(cov_raw) == 1L) {
    cov_vec   <- network::get.vertex.attribute(nw, cov_raw)
    cov_label <- cov_raw
    if (is.null(cov_vec))
      stop(termname, ": attribut inexistant: ", sQuote(cov_raw), ".")
    dbgcat("cov source = vertex attribute ", sQuote(cov_label))
  } else {
    cov_vec   <- cov_raw
    cov_label <- "cov"
    dbgcat("cov source = vector literal")
  }

  if (length(cov_vec) < n1)
    stop(termname, ": longueur de la covariée < n1.")

  cov_vec <- cov_vec[seq_len(n1)]
  dbgcat("cov length = ", length(cov_vec),
         " | head = ", paste(utils::head(as.character(cov_vec), 6L), collapse = ","))

  # ----- 3) Numeric coercion + fail-fast on NA -------------------------------
  # The covariate must be coercible to numeric. Logical and integer inputs
  # are accepted; any NA on the actor mode leads to an immediate error.
  if (is.logical(cov_vec) || is.integer(cov_vec)) {
    cov_vec <- as.numeric(cov_vec)
  }
  if (!is.numeric(cov_vec))
    stop(termname, ": la covariée doit être convertissable en numérique.")

  if (anyNA(cov_vec))
    stop(termname, ": NA non autorisé dans la covariée du mode acteurs.")

  # ----- 4) Lambda: numeric, > 1, possibly vectorized ------------------------
  # 'lambda' is interpreted as one or more geometric parameters. Each value
  # must be finite and strictly greater than 1. The initializer is vectorized
  # in lambda: L = length(lambda_vec) statistics are returned.
  lambda_raw <- a$lambda
  if (!is.numeric(lambda_raw) || length(lambda_raw) < 1L)
    stop(termname, ": 'lambda' doit être numérique (scalaire ou vecteur).")

  lambda_vec <- as.double(lambda_raw)
  if (any(!is.finite(lambda_vec)))
    stop(termname, ": 'lambda' doit être fini.")
  if (any(lambda_vec <= 1))
    stop(termname, ": toutes les valeurs de 'lambda' doivent être > 1.")

  L <- length(lambda_vec)
  dbgcat("lambda_vec = {", paste(signif(lambda_vec, 5L), collapse = ", "), "} (L = ", L, ")")

  # ----- 5) Coefficient names -------------------------------------------------
  # One coefficient name is created per lambda value, encoding the covariate
  # label and the lambda value for interpretability in model summaries.
  coef.names <- sprintf(
    "cov_diff_GW[%s]_lambda%.5g",
    cov_label,
    lambda_vec
  )
  dbgcat("coef.names = ", paste(coef.names, collapse = " | "))

  # ----- 6) Build INPUT_PARAM for the C layer --------------------------------
  # Layout of INPUT_PARAM:
  #   [1]           = n1            (actor-mode size)
  #   [2]           = L             (number of lambda values)
  #   [3..(2+L)]    = lambda_vec    (L entries)
  #   [3+L .. end]  = cov_vec[1:n1] (numeric covariate on the actor mode)
  inputs <- c(
    as.double(n1),
    as.double(L),
    as.double(lambda_vec),
    as.double(cov_vec)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " L=", L,
         " | lambda[1:6]=", paste(utils::head(signif(lambda_vec, 5L), 6L), collapse = ","),
         " | cov[1:6]=",    paste(utils::head(signif(cov_vec,    5L), 6L), collapse = ","))

  # ----- 7) Return ERGM term specification -----------------------------------
  # The 'name' field must match the C symbol implementing the change-statistic.
  # The emptynwstats vector has length L, one entry per lambda value.
  list(
    name         = "cov_diff_GW",
    coef.names   = coef.names,
    inputs       = inputs,              # n1, L, lambda[L], x[n1]
    dependence   = TRUE,
    minval       = -Inf,                # alternating combination can be negative
    maxval       = Inf,
    emptynwstats = numeric(L)
  )
}
