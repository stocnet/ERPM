#' ERGM term: cov_fulldiff (within-group range of a numeric covariate)
#' @name InitErgmTerm.cov_fulldiff
#' @aliases cov_fulldiff
#' @note InitErgmTerm.cov_fulldiff.R
#'
#' @description
#' \code{cov_fulldiff} is an ERGM term for bipartite networks that measures,
#' for each group, the dispersion of a numeric actor covariate through the
#' max–min range, optionally restricted to a subset of group sizes.
#'
#' The bipartite network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side, representing groups).
#' }
#'
#' Each actor in the actor mode carries a numeric covariate value \eqn{x_i}. For
#' each group node, we look at the actors adjacent to that group and compute
#' the within-group range \eqn{x_g^{\max} - x_g^{\min}}. An optional \code{size}
#' filter restricts which group sizes contribute to the statistic.
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic, declared in the
#' compiled code under a symbol compatible with \code{name = "cov_fulldiff"}.
#' The R initializer below:
#' \itemize{
#'   \item enforces that the network is bipartite via \code{nw \%n\% "bipartite"};
#'   \item extracts an actor-level covariate from a vertex attribute or from a
#'         literal vector, coercing it to numeric and failing fast on \code{NA};
#'   \item normalizes the optional \code{size} argument into a sorted set of
#'         strictly positive integers;
#'   \item packs the actor covariate and the size filter into a compact
#'         \code{INPUT_PARAM} vector for the C layer.
#' }
#'
#' The ERGM infrastructure will call the C change-statistic whenever a toggle
#' affects an edge between an actor and a group. The C code then recomputes the
#' local contribution of the affected group(s) to the total max–min range.
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
#'   \item \eqn{S} the set of allowed group sizes, derived from the
#'         \code{size} argument (\eqn{S = \mathbb{N}^*} if \code{size = NULL}).
#' }
#'
#' Define, for each group \eqn{g \in G} with \eqn{n_g \ge 1},
#' \deqn{
#'   x_g^{\max} = \max_{i \in A_g} x_i,
#'   \qquad
#'   x_g^{\min} = \min_{i \in A_g} x_i.
#' }
#' The statistic is
#' \deqn{
#'   T(B;x)
#'   =
#'   \sum_{g \in G}
#'     \mathbf{1}[\,n_g \in S\,]
#'     \bigl(x_g^{\max} - x_g^{\min}\bigr),
#' }
#' i.e. the sum of within-group ranges of the covariate, restricted to group
#' sizes in \eqn{S}. When \code{size = NULL}, \eqn{S} is the set of all
#' positive integers and every non-empty group contributes.
#'
#' @section Usage:
#' Typical usage with {ergm} on a bipartite network \code{nw}:
#' \preformatted{
#'   # Range of a numeric actor covariate within each group (all group sizes)
#'   summary(nw ~ cov_fulldiff(cov = "x_attr"))
#'
#'   # Restrict to groups of size 3 or 4
#'   summary(nw ~ cov_fulldiff(cov = "x_attr", size = c(3, 4)))
#' }
#'
#' When using the ERPM wrapper on a partition-based workflow, the term can be
#' invoked indirectly as:
#' \preformatted{
#'   erpm(partition ~ cov_fulldiff(cov = "x_attr", size = 2:5))
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
#' fail-fast manner. The \code{size} argument, if supplied, must contain
#' strictly positive integers; a zero-length \code{size} vector is considered
#' invalid and should be replaced by \code{NULL}.
#'
#' Internally, \code{cov_fulldiff} passes the size filter and the actor
#' covariate to the C layer through \code{INPUT_PARAM} with layout
#' \code{c(n1, L, sizes[1:L], x[1:n1])}.
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
#'   # Adjacency: actors 1..4, groups 5..7
#'   adj <- matrix(0, n_total, n_total)
#'   # Group 5 has actors 1 and 2
#'   adj[1, 5] <- adj[5, 1] <- 1
#'   adj[2, 5] <- adj[5, 2] <- 1
#'   # Group 6 has actors 3 and 4
#'   adj[3, 6] <- adj[6, 3] <- 1
#'   adj[4, 6] <- adj[6, 4] <- 1
#'   # Group 7 is empty
#'
#'   nw <- network(adj, directed = FALSE, matrix.type = "adjacency")
#'   nw \%n\% "bipartite" <- n_actors  # actor mode size
#'
#'   # Numeric covariate on the actor mode
#'   x_attr <- c(1.0, 2.0, 10.0, 11.0)
#'   set.vertex.attribute(nw, "x_attr", c(x_attr, rep(NA_real_, n_groups)))
#'
#'   # Group 5 range: 2.0 - 1.0 = 1.0
#'   # Group 6 range: 11.0 - 10.0 = 1.0
#'   # Group 7 empty: contributes 0
#'   summary(nw ~ cov_fulldiff(cov = "x_attr"))
#'
#'   # Restrict to groups of size 2
#'   summary(nw ~ cov_fulldiff(cov = "x_attr", size = 2))
#'
#'   # Fit a simple ERGM including the term
#'   fit <- ergm(nw ~ cov_fulldiff(cov = "x_attr"))
#'   summary(fit)
#' }
#'
#' @section Tests:
#' Self-tests for \code{cov_fulldiff} construct small bipartite networks with
#' known group memberships and numeric covariates, and compare:
#' \itemize{
#'   \item the ERGM summary \code{summary(nw ~ cov_fulldiff(...))};
#'   \item direct evaluations of
#'         \eqn{\sum_g \mathbf{1}[n_g \in S] (x_g^{\max} - x_g^{\min})}
#'         for various choices of the size filter \eqn{S}.
#' }
#' Additional checks verify that toggling an actor–group tie updates the
#' statistic by the expected local change in the max–min range of the affected
#' group, including cases where the group becomes empty or changes size.
#'
#' @keywords ERGM term bipartite groups covariate range
#' @md
#'
#' @export
InitErgmTerm.cov_fulldiff <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_fulldiff"

  # Debug helpers:
  # - dbg: logical flag, controlled by an R option;
  # - dbgcat(): emits prefixed debug messages when dbg is TRUE.
  dbg    <- isTRUE(getOption("ERPM.cov_fulldiff.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[cov_fulldiff][DEBUG]", ..., "\n", sep = "")

  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov",                              "size"),
    vartypes      = c("character,numeric,logical,vector", "numeric,integer"),
    defaultvalues = list(NULL,                             NULL),
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
    stop(termname, ": la covariée doit être coercible en numérique.")

  if (anyNA(cov_vec))
    stop(termname, ": NA non autorisé dans la covariée du mode acteurs.")

  # ----- 4) Size filter S (argument 'size') ----------------------------------
  # The 'size' argument defines the set S of allowed group sizes:
  # - NULL => no restriction (all group sizes allowed, L=0);
  # - numeric => converted to a sorted, unique vector of positive integers.
  sizes <- a$size
  if (is.null(sizes)) {
    L         <- 0L
    sizes_vec <- numeric(0)
    dbgcat("size filter = <ALL> (L=0)")
    size_label <- "_all"
  } else {
    if (!is.numeric(sizes))
      stop(termname, ": 'size' doit être numérique.")
    if (length(sizes) == 0L)
      stop(termname, ": 'size' vide (integer(0)) interdit. Utilisez NULL pour toutes tailles.")

    sizes <- as.integer(round(sizes))
    if (any(!is.finite(sizes)) || any(sizes <= 0L))
      stop(termname, ": 'size' doit contenir des entiers positifs.")

    sizes     <- sort(unique(sizes))
    L         <- length(sizes)
    sizes_vec <- as.double(sizes)
    dbgcat("size filter = {", paste(sizes, collapse = ","), "} (L=", L, ")")
    size_label <- sprintf("_S{%s}", paste(sizes, collapse = ","))
  }

  # ----- 5) Coefficient name --------------------------------------------------
  # Coefficient name encodes the covariate label and the size filter for
  # interpretability in model summaries.
  coef.name <- sprintf("cov_fulldiff[%s]%s", cov_label, size_label)
  dbgcat("coef.name = ", coef.name)

  # ----- 6) Build INPUT_PARAM for the C layer --------------------------------
  # Layout of INPUT_PARAM:
  #   [1]   = n1            (actor-mode size)
  #   [2]   = L             (number of size values)
  #   [3..] = sizes_vec     (L entries; may be empty if L=0)
  #   [...] = cov_vec[1:n1] (numeric covariate on the actor mode)
  inputs <- c(
    as.double(n1),
    as.double(L),
    sizes_vec,
    as.double(cov_vec)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " L=", L,
         " | cov[1:6]=", paste(utils::head(signif(cov_vec, 5L), 6L), collapse = ","))

  # ----- 7) Return ERGM term specification -----------------------------------
  list(
    name         = "cov_fulldiff",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, L, sizes[L], x[n1]
    dependence   = TRUE,
    minval       = 0,
    maxval       = Inf,
    emptynwstats = 0
  )
}
