# ==============================================================================
# File    : R/InitErgmTerm.dyadcov.R
# Term    : dyadcov
# Project : ERPM / ERGM extensions
# ==============================================================================

#' ERGM term: dyadcov (clique-based dyadic covariate)
#'
#' @name InitErgmTerm.dyadcov
#' @aliases dyadcov
#' @note InitErgmTerm.dyadcov.R
#'
#' @description
#' \code{dyadcov} is an ERGM term for bipartite networks that aggregates a
#' numeric dyadic covariate \eqn{Z = (z_{ij})} over cliques of actors within
#' each group. The network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side of the bipartite graph).
#' }
#'
#' For each group node in the group mode, we consider the set of adjacent actors
#' and all \eqn{k}-cliques of actors within that group. For a given clique
#' \eqn{C} of size \eqn{k}, the contribution is the product of dyadic covariates
#' over all actor pairs in the clique. For each unordered pair \eqn{\{i,j\}} with
#' \eqn{i<j}, the term uses the sum of both orientations:
#' \deqn{
#'   \prod_{i<j,\, i,j \in C} (z_{ij} + z_{ji}).
#' }
#' If \eqn{Z} is symmetric, each factor reduces to \eqn{2 z_{ij}}.
#'
#' Three variants are supported, controlled by the \code{normalize} argument
#' (with deprecated aliases \code{normalized} and \code{norm}):
#' \itemize{
#'   \item \emph{no normalisation} (raw sum) ;
#'   \item \emph{global normalisation} (\code{normalize = "global"}), factor
#'         \eqn{1 / n_g} per group;
#'   \item \emph{by-group normalisation} (\code{normalize = "by_group"}), factor
#'         \eqn{1 / \binom{n_g}{k}} per group.
#' }
#' Backward compatibility is preserved with the historical interface based on
#' TRUE/FALSE and the labels \code{"none"}, \code{"size"} and \code{"cliques"}.
#'
#' @details
#' The dyadic covariate \eqn{Z} is a real actor-by-actor matrix of size
#' \eqn{n_A \times n_A}, where \eqn{n_A} is the size of the actor mode
#' (\code{n1 = nw \%n\% "bipartite"}). The matrix is allowed to be non-symmetric;
#' both \eqn{z_{ij}} and \eqn{z_{ji}} are read whenever available, through the
#' symmetric combination \eqn{z_{ij} + z_{ji}}. The initializer:
#' \itemize{
#'   \item enforces bipartiteness and retrieves the actor-mode size \eqn{n_A};
#'   \item accepts \code{dyadcov} either as a literal matrix or as the name of a
#'         network-level attribute (\code{nw \%n\% "..."});
#'   \item truncates \code{dyadcov} to its top-left \code{n1 x n1} block if larger;
#'   \item checks that \code{dyadcov} is numeric and free of \code{NA};
#'   \item optionally reports (in debug mode) how far the matrix is from being
#'         symmetric, but does not require symmetry;
#'   \item parses \code{clique_size} into an integer \eqn{k \ge 2};
#'   \item parses \code{normalize} (or its aliases) into a normalisation mode.
#' }
#'
#' @section Core statistic:
#' The core statistic (for a fixed \eqn{k}) can be written in terms of the
#' partition \eqn{p} and the dyadic covariate \eqn{Z}. Let
#' \eqn{S_g^{(k)}(Z)} be the raw clique-based sum in group \eqn{g}:
#' \deqn{
#'   S_g^{(k)}(Z)
#'     = \sum_{C \in C_k(g)} \prod_{i<j,\, i,j \in C} (z_{ij} + z_{ji}),
#' }
#' and \eqn{n_g} the number of actors adjacent to group \eqn{g}. Then:
#'
#'   - no normalisation (backward-compatible with \code{normalize = FALSE}) :
#' \deqn{
#'       T^{(k)}(p; Z)
#'         = \sum_g S_g^{(k)}(Z)
#' }
#'
#'   - global normalisation (backward-compatible with \code{normalize = "global"}
#'     and \code{normalize = TRUE}) :
#' \deqn{
#'       T^{(k)}_{\mathrm{global}}(p; Z)
#'         = \sum_g \mathbf{1}(n_g \ge k) \frac{1}{n_g}
#'             S_g^{(k)}(Z)
#' }
#'
#'   - by-group normalisation:
#' \deqn{
#'       T^{(k)}_{\mathrm{by\_group}}(p; Z)
#'         = \sum_g \mathbf{1}(n_g \ge k) \frac{1}{\binom{n_g}{k}}
#'             S_g^{(k)}(Z).
#' }
#'
#' For \eqn{k = 2} and no normalisation, this reduces to:
#' \deqn{
#'   T^{(2)}(p; Z)
#'     = \sum_g \sum_{i<j,\, i,j \in A(g)} (z_{ij} + z_{ji}),
#' }
#' which matches the definition of \code{dyadcov_full(dyadcov, size = NULL)}
#' aggregating the dyadic covariate over all ordered pairs of co-grouped actors.
#' For \eqn{k = 2} and global normalisation, each group contribution is
#' additionally scaled by \eqn{1 / n_g}.
#'
#' @section Implementation and change-statistic:
#' The term is implemented as a native ERGM C change-statistic, exposed under the
#' symbol \code{c_dyadcov}. The R initializer below:
#' \itemize{
#'   \item packages \code{n1}, \code{k}, the normalisation mode and the flattened
#'         \code{dyadcov} matrix into \code{INPUT_PARAM};
#'   \item declares the term as dependent (\code{dependence = TRUE}) with no
#'         finite \code{minval}/\code{maxval};
#'   \item sets the empty-network statistic to \code{0}.
#' }
#'
#' On each toggle of an actor–group edge, the C change-statistic recomputes the
#' local contribution for the unique group affected by the toggle, by updating
#' the clique-based sums according to the new set of actors in that group and,
#' when a normalisation is requested, dividing by either \eqn{n_g} or
#' \eqn{\binom{n_g}{k}} depending on the chosen mode.
#'
#' @section Arguments:
#' The initializer is not called directly by users; it is invoked automatically
#' by {ergm} when the term \code{dyadcov(...)} appears on the right-hand side of
#' a model formula.
#'
#' @param dyadcov matrix or character. Either:
#'   \itemize{
#'     \item a numeric matrix of size at least \code{n1 x n1}, where
#'           \code{n1 = nw \%n\% "bipartite"} is the actor-mode size; or
#'     \item the name of a network-level attribute containing such a matrix
#'           (retrieved as \code{nw \%n\% dyadcov}).
#'   }
#'   In both cases the matrix is truncated, if necessary, to its top-left
#'   \code{n1 x n1} block.
#' @param clique_size numeric or integer scalar. Target clique size \eqn{k}.
#'   It is rounded and coerced to an integer and must satisfy \eqn{k \ge 2}.
#'   Defaults to \code{2}.
#' @param normalize logical or character scalar. Controls the normalisation of
#'   the group-level sums \eqn{S_g^{(k)}(Z)}:
#'   \itemize{
#'     \item \code{FALSE} or \code{"none"}: no normalisation (raw sum) ;
#'     \item \code{"global"} or \code{TRUE} or legacy \code{"size"}:
#'           global normalisation by \eqn{1 / n_g} when \eqn{n_g \ge k};
#'     \item \code{"by_group"} or legacy \code{"cliques"}:
#'           normalisation by \eqn{1 / \binom{n_g}{k}} when \eqn{n_g \ge k}.
#'   }
#'   Deprecated aliases \code{normalized} and \code{norm} are also accepted and
#'   mapped internally to \code{normalize}.
#'
#' @return
#' A standard {ergm} term initialization list with components:
#' \itemize{
#'   \item \code{name}         = \code{"dyadcov"};
#'   \item \code{coef.names}   = a single coefficient name encoding
#'         \code{dyadcov} label, \code{k} and the normalization mode;
#'   \item \code{inputs}       = numeric vector
#'         \code{c(n1, k, norm_mode, as.double(Z))};
#'   \item \code{dependence}   = \code{TRUE};
#'   \item \code{minval}       = \code{-Inf};
#'   \item \code{maxval}       = \code{Inf};
#'   \item \code{emptynwstats} = \code{0}.
#' }
#'
#' @note
#' \itemize{
#'   \item The network must be bipartite and interpreted as actors versus groups.
#'   \item The actor mode is identified by \code{nw \%n\% "bipartite"} and must
#'         be a positive integer.
#'   \item The dyadic matrix must be numeric and free of \code{NA} values.
#'         Any \code{NA} triggers a fail-fast error. Symmetry is not required;
#'         when the matrix is not symmetric, both \eqn{z_{ij}} and \eqn{z_{ji}}
#'         are used through the sum \eqn{z_{ij} + z_{ji}}.
#'   \item Debug logging is controlled by the global option
#'         \code{options(ERPM.dyadcov.debug = TRUE/FALSE)}. When \code{TRUE},
#'         the initializer prints diagnostic messages prefixed by
#'         \code{"[dyadcov][DEBUG]"} and may report how far the dyadic matrix
#'         is from symmetry.
#' }
#'
#' @examples
#' \dontrun{
#'   library(network)
#'   library(ergm)
#'
#'   # -------------------------------------------------------------------------
#'   # Build a bipartite network: 4 actors, 2 groups
#'   # -------------------------------------------------------------------------
#'   n_actors <- 4
#'   n_groups <- 2
#'   n_total  <- n_actors + n_groups
#'
#'   adj <- matrix(0, n_total, n_total)
#'
#'   # Actors = 1..4, Groups = 5..6
#'   # Group 5 has actors 1, 2, 3
#'   adj[1, 5] <- adj[5, 1] <- 1
#'   adj[2, 5] <- adj[5, 2] <- 1
#'   adj[3, 5] <- adj[5, 3] <- 1
#'   # Group 6 has actors 3, 4
#'   adj[3, 6] <- adj[6, 3] <- 1
#'   adj[4, 6] <- adj[6, 4] <- 1
#'
#'   nw <- network(adj, directed = FALSE, matrix.type = "adjacency")
#'   nw %n% "bipartite" <- n_actors  # actor-mode size
#'
#'   # -------------------------------------------------------------------------
#'   # Build a dyadic covariate on actors (symmetric or not)
#'   # -------------------------------------------------------------------------
#'   Z <- matrix(
#'     c(
#'       0, 1,  2,  3,
#'       10, 0, 4,  5,
#'       20, 40, 0, 6,
#'       30, 50, 60, 0
#'     ),
#'     n_actors, n_actors, byrow = TRUE
#'   )
#'   nw %n% "Z_example" <- Z
#'
#'   # -------------------------------------------------------------------------
#'   # Inspect the dyadcov statistic
#'   # -------------------------------------------------------------------------
#'   # k = 2, raw sum over all cliques of size 2 (pairs) within groups,
#'   # using (z_ij + z_ji) for each pair {i,j}
#'   summary(nw ~ dyadcov("Z_example", clique_size = 2, normalize = FALSE),
#'           constraints = ~ b1part)
#'
#'   # k = 3, global normalisation by 1 / n_g on 3-cliques
#'   summary(nw ~ dyadcov("Z_example", clique_size = 3, normalize = "global"),
#'           constraints = ~ b1part)
#'
#'   # k = 3, by-group normalisation by 1 / C(n_g, 3)
#'   summary(nw ~ dyadcov("Z_example", clique_size = 3,
#'                        normalize = "by_group"),
#'           constraints = ~ b1part)
#'
#'   # -------------------------------------------------------------------------
#'   # Example with the ERPM wrapper (network already bipartite)
#'   # -------------------------------------------------------------------------
#'   # erpm(nw ~ dyadcov("Z_example", clique_size = 2))
#' }
#'
#' @section Tests:
#' Self-tests for \code{dyadcov} (not shown here) typically:
#' \itemize{
#'   \item build small bipartite networks with a known partition of actors into
#'         groups;
#'   \item define dyadic matrices \code{Z} (symmetric or not) so that clique
#'         products can be computed analytically using \code{z_ij + z_ji};
#'   \item compare \code{summary(nw ~ dyadcov(...), constraints = ~ b1part)} with
#'         a direct implementation of
#'         \code{sum_g f(n_g) * sum_{C in C_k(g)} prod_{i<j in C}(Z[i,j]+Z[j,i])},
#'         where \code{f(n_g)} is \code{1}, \code{1 / n_g} or
#'         \code{1 / choose(n_g, k)} depending on the chosen normalisation;
#'   \item verify that toggling a single actor–group edge changes the statistic
#'         by the local difference between the "before" and "after" clique sums
#'         (with or without normalisation), in agreement with the C
#'         change-statistic \code{c_dyadcov}.
#' }
#'
#' @keywords ERGM term bipartite dyadic covariate cliques
#' @md

InitErgmTerm.dyadcov <- function(nw, arglist, ...) {
  termname <- "dyadcov"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  # Global option:
  #   options(ERPM.dyadcov.debug = TRUE/FALSE)
  # When TRUE, the initializer prints diagnostic messages to the console.
  dbg    <- isTRUE(getOption("ERPM.dyadcov.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[dyadcov][DEBUG]", ..., "\n", sep = "")

  # ---------------------------------------------------------------------------
  # Base ERGM term validation and argument parsing
  # ---------------------------------------------------------------------------
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("dyadcov",                 "clique_size",
                      "normalize",               "normalized", "norm"),
    vartypes      = c("matrix,character",        "numeric,integer",
                      "logical,character",       "logical,character",
                      "logical,character"),
    defaultvalues = list(NULL,                   2,
                         NULL,                   NULL,         NULL),
    required      = c(TRUE,                      FALSE,
                      FALSE,                     FALSE,        FALSE)
  )

  # ---------------------------------------------------------------------------
  # Actor-mode size (n1) from the bipartite attribute
  # ---------------------------------------------------------------------------
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L)
    stop(termname, ": réseau non biparti strict (attribut %n% 'bipartite' manquant ou invalide).")
  dbgcat("n1 = ", n1)

  # ---------------------------------------------------------------------------
  # Retrieve the dyadic matrix (n1 x n1)
  # ---------------------------------------------------------------------------
  dyad_raw   <- a$dyadcov
  dyad_label <- NULL

  if (is.character(dyad_raw) && length(dyad_raw) == 1L) {
    # Case: name of a network-level attribute storing the matrix
    dyad_mat   <- nw %n% dyad_raw
    dyad_label <- dyad_raw
    if (is.null(dyad_mat))
      stop(termname, ": attribut de niveau réseau inexistant: ", sQuote(dyad_raw), ".")
    dbgcat("dyadcov source = network attribute ", sQuote(dyad_label))
  } else {
    # Case: matrix passed literally as an argument
    dyad_mat   <- dyad_raw
    dyad_label <- "dyadcov"
    dbgcat("dyadcov source = matrix literal")
  }

  if (!is.matrix(dyad_mat))
    stop(termname, ": 'dyadcov' doit être une matrice ou le nom d'un attribut de niveau réseau.")

  nr <- nrow(dyad_mat)
  nc <- ncol(dyad_mat)

  if (nr < n1 || nc < n1)
    stop(termname, ": dimensions de la matrice dyadique (", nr, "x", nc,
         ") insuffisantes pour n1 = ", n1, ".")

  # If larger than needed, restrict to the top-left n1 x n1 block
  if (nr > n1 || nc > n1) {
    dyad_mat <- dyad_mat[seq_len(n1), seq_len(n1), drop = FALSE]
    dbgcat("dyadcov truncated to ", n1, "x", n1)
  }

  # ---------------------------------------------------------------------------
  # Numeric coercion and fail-fast on NA
  # ---------------------------------------------------------------------------
  if (!is.numeric(dyad_mat))
    stop(termname, ": la matrice dyadique doit être numérique.")

  if (anyNA(dyad_mat))
    stop(termname, ": NA non autorisé dans la matrice dyadique.")

  # ---------------------------------------------------------------------------
  # Optional symmetry diagnostics (non-blocking)
  # ---------------------------------------------------------------------------
  tol <- 1e-8
  max_asym <- max(abs(dyad_mat - t(dyad_mat)))
  if (dbg && max_asym > tol) {
    dbgcat("warning: dyadcov matrix is not symmetric, max |Z - t(Z)| = ",
           signif(max_asym, 5L))
  }

  dbgcat("dyadcov dim = ", paste(dim(dyad_mat), collapse = "x"),
         " | sample = ",
         paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L), collapse = ","))

  # ---------------------------------------------------------------------------
  # clique_size (k)
  # ---------------------------------------------------------------------------
  k_raw <- a$clique_size
  if (is.null(k_raw) || length(k_raw) == 0L) {
    k <- 2L
  } else {
    if (!is.numeric(k_raw))
      stop(termname, ": 'clique_size' doit être numérique ou entier.")
    k <- as.integer(round(k_raw[1L]))
  }

  if (!is.finite(k) || k < 2L)
    stop(termname, ": 'clique_size' doit être un entier >= 2.")

  dbgcat("clique_size (k) = ", k)

  # ---------------------------------------------------------------------------
  # normalize / normalized / norm (normalisation mode)
  # ---------------------------------------------------------------------------
  norm_raw <- NULL
  if (!is.null(a$normalize)) {
    norm_raw <- a$normalize
  } else if (!is.null(a$normalized)) {
    norm_raw <- a$normalized
  } else if (!is.null(a$norm)) {
    norm_raw <- a$norm
  } else {
    norm_raw <- FALSE
  }

  norm_mode  <- 0L
  norm_label <- "none"

  if (is.logical(norm_raw)) {
    # backward-compatible behaviour
    if (isTRUE(norm_raw[1L])) {
      norm_mode  <- 1L   # global 1/n_g
      norm_label <- "global"
    } else {
      norm_mode  <- 0L
      norm_label <- "none"
    }
  } else if (is.character(norm_raw)) {
    val <- norm_raw[1L]
    val <- match.arg(val, c("none",
                            "global", "by_group",
                            "size", "cliques"))
    if (val == "none") {
      norm_mode  <- 0L
      norm_label <- "none"
    } else if (val %in% c("global", "size")) {
      norm_mode  <- 1L   # global 1/n_g
      norm_label <- "global"
    } else { # "by_group" or "cliques"
      norm_mode  <- 2L   # 1 / C(n_g, k)
      norm_label <- "by_group"
    }
  } else {
    stop(termname, ": 'normalize' doit être logique ou caractère ",
         "(\"none\", \"global\", \"by_group\").")
  }

  dbgcat("normalized mode = ", norm_label, " (code=", norm_mode, ")")

  # ---------------------------------------------------------------------------
  # Coefficient name
  # ---------------------------------------------------------------------------
  # Encode the dyadic matrix label, the clique size, and the normalization mode.
  base_label <- sprintf("dyadcov[%s]_k%d", dyad_label, k)
  coef.name  <- switch(
    as.character(norm_mode),
    "0" = base_label,
    "1" = paste0(base_label, "_global"),   # 1 / n_g
    "2" = paste0(base_label, "_bygrp"),    # 1 / C(n_g, k)
    base_label
  )
  dbgcat("coef.name = ", coef.name)

  # ---------------------------------------------------------------------------
  # Build INPUT_PARAM for the C side
  # ---------------------------------------------------------------------------
  # Layout:
  #   inputs = c(
  #     as.double(n1),
  #     as.double(k),
  #     as.double(norm_mode),
  #     as.double(Z[1]), ..., as.double(Z[n1*n1])
  #   )
  # where as.double(matrix) uses column-major order, consistent with C indexing.
  inputs <- c(
    as.double(n1),
    as.double(k),
    as.double(norm_mode),
    as.double(dyad_mat)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " k=", k, " norm_mode=", norm_mode,
         " | Z[1:6]=", paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L),
                             collapse = ","))

  # ---------------------------------------------------------------------------
  # Standard ERGM term initialization return value
  # ---------------------------------------------------------------------------
  list(
    name         = "dyadcov",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, k, norm_mode, then Z[n1*n1]
    dependence   = TRUE,
    minval       = -Inf,
    maxval       = Inf,
    emptynwstats = 0
  )
}