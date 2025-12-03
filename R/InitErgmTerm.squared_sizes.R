# ==============================================================================
# File    : InitErgmTerm.squared_sizes.R
# Purpose : Declare the ERGM term 'squared_sizes' for bipartite networks
#           (aggregation over the group mode).
# Project : ERPM / ERGM extensions
# ============================================================================

#' ERGM term: squared_sizes (group-mode degrees raised to a power)
#'
#' @name InitErgmTerm.squared_sizes
#' @aliases squared_sizes
#' @note InitErgmTerm.squared_sizes.R
#'
#' @description
#' \code{squared_sizes} is an ERGM term for bipartite networks that aggregates
#' powers of group sizes. The network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side of the bipartite graph).
#' }
#'
#' For each node in the group mode with degree \eqn{\deg(g)} (number of adjacent
#' actors), the term can accumulate contributions of the form
#' \deqn{
#'   \deg(g)^{\text{pow}},
#' }
#' restricted to a set of admissible group sizes. The argument \code{sizes}
#' specifies a (possibly multi-valued) set of group sizes over which the
#' contributions are aggregated into a single statistic.
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic, declared in the
#' compiled code under the name \code{c_squared_sizes}. The R initializer:
#' \itemize{
#'   \item enforces that the network is bipartite (actor mode vs group mode);
#'   \item parses the user arguments \code{sizes} and \code{pow};
#'   \item checks that \code{sizes} are integer group sizes \eqn{\ge 1};
#'   \item checks that \code{pow} is a single integer and at least 1;
#'   \item builds a compact \code{inputs} vector for the C code encoding:
#'         \code{pow}, the number of admissible sizes, then the list of sizes.
#' }
#'
#' The ERGM infrastructure will call the C change-statistic whenever a toggle
#' affects an actor–group edge. On the C side, the code is responsible for
#' updating the scalar statistic based on how the degrees of group-mode
#' nodes change.
#'
#' @section Mathematical definition:
#' Let:
#' \itemize{
#'   \item \eqn{A} denote the set of actor-mode nodes;
#'   \item \eqn{G} denote the set of group-mode nodes;
#'   \item \eqn{y} the adjacency matrix between actors and groups;
#'   \item \eqn{\deg(g)} the degree of group node \eqn{g \in G}, i.e. the number
#'         of actors in that group;
#'   \item \eqn{S} the set of admissible group sizes specified by \code{sizes}.
#' }
#' For a given pair \code{(sizes, pow)}, the statistic is
#' \deqn{
#'   T(y)
#'   =
#'   \sum_{g \in G}
#'     \mathbf{1}\big[\deg(g) \in S\big]\,
#'     \deg(g)^{\text{pow}}.
#' }
#' When \code{sizes} is left unspecified, it defaults to all admissible
#' group sizes \eqn{1, 2, \dots, n_1}, where \eqn{n_1} is the actor-mode size,
#' so that \eqn{T(y) = \sum_{g \in G} \deg(g)^{\text{pow}}}.
#'
#' @section Usage:
#' Typical usage with {ergm} (constraints explicit):
#' \preformatted{
#'   # Sum of squared group sizes over all non-empty groups
#'   summary(nw ~ squared_sizes, constraints = ~ b1part)
#'
#'   # Sum of cubic sizes over groups of size 2 only
#'   summary(
#'     nw ~ squared_sizes(sizes = 2, pow = 3),
#'     constraints = ~ b1part
#'   )
#'
#'   # Single statistic aggregating groups of size 1 or 5 (squared sizes)
#'   summary(
#'     nw ~ squared_sizes(sizes = c(1, 5), pow = 2),
#'     constraints = ~ b1part
#'   )
#' }
#'
#' When using the ERPM wrapper:
#' \preformatted{
#'   # On an explicit bipartite network
#'   erpm(nw ~ squared_sizes)
#'
#'   # From a partition vector (actor -> group id)
#'   erpm(partition ~ squared_sizes)
#' }
#'
#' @param sizes numeric, integer(s) \eqn{\ge 1}. Admissible group sizes.
#'   If \code{sizes} is \code{NULL} or missing, it defaults to
#'   \code{1:network.bipartite(nw)}, i.e. all non-empty group sizes up to the
#'   number of actors. Can be a scalar or a vector; all values are aggregated
#'   into a \emph{single} statistic.
#' @param pow numeric, integer \eqn{\ge 1}. Power applied to the group
#'   degrees. Defaults to \code{2}. Must be of length 1; the same exponent is
#'   applied to all sizes in \code{sizes}.
#'
#' @return
#' A standard {ergm} term initialization list with components:
#' \itemize{
#'   \item \code{name}         = \code{"squared_sizes"};
#'   \item \code{coef.names}   = a single character string encoding the set
#'         of sizes and the power;
#'   \item \code{inputs}       = numeric vector
#'         \code{c(pow, length(sizes), sizes)};
#'   \item \code{dependence}   = \code{TRUE};
#'   \item \code{emptynwstats} = scalar \code{0}.
#' }
#'
#' @note
#' \itemize{
#'   \item The network must be bipartite and interpreted as actors versus groups.
#'   \item The actor mode is identified by \code{nw \%n\% "bipartite"}.
#'   \item The group mode is the complement of the actor mode.
#'   \item The initializer enforces bipartiteness via \code{check.ErgmTerm(..., bipartite = TRUE)}.
#'   \item The term is declared as dependent (\code{dependence = TRUE}), since it
#'         depends on the current pattern of actor–group memberships.
#' }
#'
#' @examples
#' \dontrun{
#'   library(network)
#'   library(ergm)
#'
#'   # Build a small bipartite network: 4 actors, 3 groups
#'   n_actors <- 4
#'   n_groups <- 3
#'   n_total  <- n_actors + n_groups
#'
#'   adj <- matrix(0, n_total, n_total)
#'   # Actors = 1..4, Groups = 5..7
#'   adj[1, 5] <- adj[5, 1] <- 1
#'   adj[2, 5] <- adj[5, 2] <- 1
#'   adj[3, 6] <- adj[6, 3] <- 1
#'   adj[4, 6] <- adj[6, 4] <- 1
#'   # Group 7 remains empty
#'
#'   nw <- network(adj, directed = FALSE, matrix.type = "adjacency")
#'   nw %n% "bipartite" <- n_actors  # actor-mode size
#'
#'   # Sum of squared group sizes over all non-empty groups
#'   summary(nw ~ squared_sizes, constraints = ~ b1part)
#'
#'   # Single statistic: squared sizes for groups of size 1 or 3
#'   summary(
#'     nw ~ squared_sizes(sizes = c(1, 3),
#'                        pow   = 2),
#'     constraints = ~ b1part
#'   )
#'
#'   # ERGM fit using squared_sizes
#'   fit <- ergm(nw ~ squared_sizes, constraints = ~ b1part)
#'   summary(fit)
#' }
#'
#' @note
#' Self-tests for \code{squared_sizes} construct small bipartite networks with
#' known group degrees and compare:
#' \itemize{
#'   \item the ERGM summary \code{summary(nw ~ squared_sizes(...), constraints = ~ b1part)};
#'   \item a direct computation of
#'         \code{sum( (deg_g \%in\% sizes) * deg_g^pow )};
#' }
#' They also verify that toggling a single actor–group edge changes the
#' statistic \eqn{T(y)} by the expected local increment, in agreement with
#' the C change-statistic \code{c_squared_sizes}.
#'
#' @keywords ERGM term bipartite groups degree power
#' @md
#'
#' @export
InitErgmTerm.squared_sizes <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "squared_sizes"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  # Global option:
  #   options(ERPM.squared_sizes.debug = TRUE/FALSE)
  # When TRUE, the initializer prints diagnostic messages to the console.
  dbg    <- isTRUE(getOption("ERPM.squared_sizes.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[squared_sizes][DEBUG]", ..., "\n", sep = "")

  dbgcat("InitErgmTerm.squared_sizes called with args: ",
         paste(names(arglist), collapse = ", "))

  # Guard:  typo 'size' instead of 'sizes'
  if ("size" %in% names(arglist) && !"sizes" %in% names(arglist)) {
    ergm_Init_stop(
      sQuote(termname),
      ": argument 'size' is not supported; did you mean 'sizes'?"
    )
  }

  # ---------------------------------------------------------------------------
  # Base validation and structural requirements
  # ---------------------------------------------------------------------------
  # - Enforce a bipartite network (actor mode / group mode).
  # - Declare the allowed arguments: sizes, pow.
  # - Provide default values: sizes = NULL (all sizes), pow = 2.
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,                     # require a bipartite actor–group representation
    varnames      =    c("sizes",  "pow"),
    vartypes      =    c("numeric","numeric"),
    defaultvalues = list(NULL,     2),        # default: all non-empty groups, squared sizes
    required      =    c(FALSE,    FALSE)
  )

  sizes <- a$sizes
  pow   <- a$pow

  # ---------------------------------------------------------------------------
  # Default sizes: all non-empty group sizes up to the actor-mode size
  # ---------------------------------------------------------------------------
  # If sizes is NULL or empty, we take all possible sizes 1..n1.
  if (is.null(sizes) || length(sizes) == 0L) {
    # On lit explicitement l'attribut "bipartite" du network, pour éviter tout
    # masquage éventuel d'une fonction network.bipartite() dans l'environnement.
    n1 <- network::get.network.attribute(nw, "bipartite")
    if (is.null(n1) || is.na(n1)) {
      ergm_Init_stop(
        sQuote(termname),
        ": network must have a valid 'bipartite' attribute (actor-mode size)."
      )
    }
    n1 <- as.integer(n1)
    if (!is.finite(n1) || n1 < 1L) {
      ergm_Init_stop(
        sQuote(termname),
        ": 'bipartite' attribute must be a positive finite integer."
      )
    }
    dbgcat("bipartite attribute (n1) =", n1,
           " | using default sizes = 1:", n1)
    sizes <- seq_len(n1)
  } else {
    dbgcat("user-specified sizes (raw) = ",
           paste(sizes, collapse = ","))
  }

  # ---------------------------------------------------------------------------
  # Constraints and normalization on sizes
  # ---------------------------------------------------------------------------
  sizes <- as.numeric(sizes)
  if (any(is.na(sizes))) {
    ergm_Init_stop(sQuote(termname), ": 'sizes' must not contain NA.")
  }
  if (any(sizes < 1 | sizes != as.integer(sizes))) {
    ergm_Init_stop(sQuote(termname), ": 'sizes' must be integer >= 1.")
  }
  sizes <- as.integer(sizes)
  dbgcat("sizes (validated) = ", paste(sizes, collapse = ","))

  # ---------------------------------------------------------------------------
  # Constraints on pow (scalar)
  # ---------------------------------------------------------------------------
  if (length(pow) == 0L) {
    pow <- 2L
  }
  if (length(pow) > 1L) {
    ergm_Init_stop(sQuote(termname), ": 'pow' must be of length 1.")
  }
  if (any(pow < 1 | pow != as.integer(pow))) {
    ergm_Init_stop(sQuote(termname), ": 'pow' must be an integer >= 1.")
  }
  pow <- as.integer(pow[1L])
  dbgcat("pow (validated)   = ", pow)

  # ---------------------------------------------------------------------------
  # Trivial case: no sizes => no statistic
  # ---------------------------------------------------------------------------
  if (length(sizes) == 0L) {
    dbgcat("no sizes after validation -> returning NULL term")
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Coefficient name
  # ---------------------------------------------------------------------------
  # Build a readable coefficient name encoding:
  # - the set of group sizes 'sizes';
  # - the power 'pow', when pow != 1.
  if (length(sizes) == 1L) {
    size_tag <- paste0("size", sizes)
  } else {
    size_tag <- paste0("sizes", paste(sizes, collapse = "_"))
  }
  coef.name <- paste0(
    "squared_", size_tag,
    ifelse(pow != 1L, paste0("_pow", pow), "")
  )
  dbgcat("coef.names = ", coef.name)

  # ---------------------------------------------------------------------------
  # Compact INPUT_PARAM vector for the C change-statistic
  # ---------------------------------------------------------------------------
  # Layout:
  #   inputs = c(
  #     pow,
  #     length(sizes),
  #     sizes[1], ..., sizes[length(sizes)]
  #   )
  inputs <- c(
    as.double(pow),
    as.double(length(sizes)),
    as.double(sizes)
  )
  dbgcat("inputs length = ", length(inputs),
         " (1 aggregated stat over ", length(sizes), " sizes)")

  # ---------------------------------------------------------------------------
  # Standard ERGM term initialization return value
  # ---------------------------------------------------------------------------
  # - name         : must match the C symbol (without the 'c_' prefix parsed by {ergm}).
  # - coef.names   : single label for the aggregated statistic.
  # - inputs       : numeric vector passed to the C change-statistic as INPUT_PARAM.
  # - dependence   : TRUE since the term depends on the network configuration.
  # - emptynwstats : empty network -> scalar 0.
  list(
    name         = "squared_sizes",          # must match CHANGESTAT_FN(c_squared_sizes) (without 'c_')
    coef.names   = coef.name,
    inputs       = inputs,
    dependence   = TRUE,
    emptynwstats = 0                         # empty network => 0
  )
}