# ==============================================================================
# File    : R/InitErgmTerm.cliques.R
# Purpose : Declare the ERGM term 'cliques' for bipartite actor-group networks
#           (counts k-actor cliques induced by group sizes).
# Project : ERPM / ERGM extensions
# ==============================================================================

#' ERGM term: cliques (k-actor cliques via group sizes)
#'
#' @name InitErgmTerm.cliques
#' @aliases cliques
#' @note InitErgmTerm.cliques.R
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' \code{cliques} is an ERGM term for bipartite actor-group networks that counts
#' k-actor cliques induced by group memberships. The bipartite network is
#' interpreted as:
#' \itemize{
#'   \item an \emph{actor mode}, whose size is given by \code{nw \%n\% "bipartite"};
#'   \item a \emph{group mode}, consisting of the remaining nodes that represent
#'         groups.
#' }
#'
#' Each group node in the group mode has a degree \eqn{n_g} (number of adjacent
#' actors). Interpreting each group as forming a complete clique among its
#' actors, the total number of k-actor cliques is
#' \deqn{
#'   T_k(y) = \sum_{g \in \text{group mode}} \binom{n_g}{k}.
#' }
#' For k = 1 this reduces to the number of groups of size 1.
#'
#' The term \code{cliques} computes this statistic directly from group sizes,
#' without explicitly materializing the actor-actor projection.
#'
#' The initializer supports:
#' \itemize{
#'   \item \eqn{k \ge 1}, where:
#'     \itemize{
#'       \item for \eqn{k \ge 2}, \eqn{T_k(y)} is the number of k-actor cliques;
#'       \item for \eqn{k = 1}, \eqn{T_1(y)} is the number of groups of size 1.
#'     }
#'   \item an optional group-size-based normalization that rescales \eqn{T_k(y)}
#'         by the group sizes \eqn{n_g}.
#' }
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic.
#'
#' IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
#' - This term MUST support multi-toggle moves (swap/split/merge represented as
#'   a list of membership-edge toggles).
#' - Therefore, the compiled change-statistic is implemented using the
#'   D_CHANGESTAT_FN API (multi-toggle).
#' - On the R side, we MUST advertise this to ergm by returning `d_func = TRUE`.
#'   Otherwise ergm will try to call the changestat as a one-toggle C_CHANGESTAT_FN,
#'   causing a signature mismatch and typically a segfault.
#'
#' Compiled symbol naming convention:
#' - Recommended: implement the C function as `d_cliques` via
#'   D_CHANGESTAT_FN(d_cliques).
#' - Avoid exposing a symbol named `c_cliques` with a D-signature, because ergm
#'   may resolve it as a one-toggle entrypoint and crash.
#'
#' The R initializer:
#' \itemize{
#'   \item enforces bipartite network;
#'   \item normalizes argument names (positional, clique_size -> k);
#'   \item accepts one or several values of \code{k} (vectorized interface);
#'   \item optionally selects group-size-normalized mode via a sign-flag on scale;
#'   \item packs \code{k} and \code{scale} into \code{INPUT_PARAM} for the C layer.
#' }
#'
#' @keywords ERGM term bipartite groups cliques
#' @md
#'
#' @export
InitErgmTerm.cliques <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cliques"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  # Global option:
  #   options(ERPM.cliques.debug = TRUE/FALSE)
  dbg    <- isTRUE(getOption("ERPM.cliques.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[cliques][DEBUG]", ..., "\n", sep = "")

  dbgcat("InitErgmTerm.cliques called with args: ", paste(names(arglist), collapse = ", "))

  # ---------------------------------------------------------------------------
  # Normalize user arguments so that the initializer consistently sees 'k'
  # ---------------------------------------------------------------------------
  # - cliques(1)             -> k = 1
  # - cliques(clique_size=1) -> k = 1
  # - cliques(k=1)           -> k = 1
  if (length(arglist) == 1L) {
    nm <- names(arglist)
    if (is.null(nm) || isTRUE(nm[1L] == "")) {
      # Single positional argument: cliques(1)
      arglist <- list(k = arglist[[1L]])
    }
  }

  # Backward compatibility: accept 'clique_size' and rename it to 'k'
  if (!is.null(names(arglist)) && "clique_size" %in% names(arglist)) {
    arglist[["k"]] <- arglist[["clique_size"]]
    arglist[["clique_size"]] <- NULL
  }

  # Guard: common typo "size" instead of "k"
  if (!is.null(names(arglist)) && "size" %in% names(arglist) && !"k" %in% names(arglist)) {
    ergm_Init_stop(sQuote(termname), ": argument 'size' is not supported; did you mean 'k'?")
  }

  # ---------------------------------------------------------------------------
  # Base validation and structural requirements
  # ---------------------------------------------------------------------------
  # - Enforce a bipartite network (actor mode / group mode).
  # - Declare allowed arguments: k, normalized.
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("k", "normalized"),
    vartypes      = c("numeric", "logical"),
    defaultvalues = list(2, FALSE),
    required      = c(FALSE, FALSE)
  )

  k  <- a$k
  nz <- a$normalized

  # ---------------------------------------------------------------------------
  # Validate k and normalized
  # ---------------------------------------------------------------------------
  if (length(k) < 1L)
    ergm_Init_stop(sQuote(termname), ": specify at least one value of k.")
  if (any(is.na(k)))
    ergm_Init_stop(sQuote(termname), ": 'k' must not contain NA.")
  if (any(k != as.integer(k) | k < 1))
    ergm_Init_stop(sQuote(termname), ": 'k' must contain integers >= 1.")
  if (length(nz) != 1L || is.na(nz))
    ergm_Init_stop(sQuote(termname), ": 'normalized' must be a scalar boolean.")

  k <- as.integer(k)
  nz <- isTRUE(nz)

  # ---------------------------------------------------------------------------
  # Retrieve the actor-mode size N_A from the bipartite attribute
  # ---------------------------------------------------------------------------
  n1 <- network::get.network.attribute(nw, "bipartite")
  if (is.null(n1) || is.na(n1))
    ergm_Init_stop(sQuote(termname), ": non-bipartite network or missing 'bipartite' attribute.")
  n1 <- as.integer(n1)
  if (!is.finite(n1) || n1 < 1L)
    ergm_Init_stop(sQuote(termname), ": invalid bipartite attribute (must be positive integer).")

  dbgcat("bipartite attribute (n1) =", n1)
  dbgcat("k (validated) =", paste(k, collapse = ","), " | normalized =", nz)

  # ---------------------------------------------------------------------------
  # Prepare scaling / mode flags for the C layer
  # ---------------------------------------------------------------------------
  # - if normalized = FALSE: scale_j > 0 (raw T_k statistic);
  # - if normalized = TRUE : scale_j < 0 (group-size-normalized statistic),
  #   with the absolute value used as an extra multiplicative factor.
  scale <- rep(1, length(k))
  if (nz) scale <- rep(-1, length(k))

  # ---------------------------------------------------------------------------
  # Coefficient names and INPUT_PARAM layout
  # ---------------------------------------------------------------------------
  # - one coefficient per k;
  # - INPUT_PARAM = (k_1, scale_1, k_2, scale_2, ...).
  coef.names <- if (nz) paste0("cliques_k", k, "_grp") else paste0("cliques_k", k)
  inputs <- c(rbind(as.integer(k), as.double(scale)))

  dbgcat("coef.names =", paste(coef.names, collapse = " | "))
  dbgcat("inputs length =", length(inputs), " (", length(k), " stats)")

  # ---------------------------------------------------------------------------
  # Standard ERGM term initialization return value
  # ---------------------------------------------------------------------------
  # IMPORTANT:
  # - `d_func = TRUE` tells ergm to call the multi-toggle (D_) changestat entrypoint.
  # - Without it, ergm assumes a one-toggle C_ changestat and will call the function
  #   with the wrong signature if you compiled only a D_ function (=> segfault).
  list(
    name         = "cliques",
    coef.names   = coef.names,
    inputs       = inputs,
    dependence   = TRUE,
    d_func       = TRUE,                 # <-- REQUIRED for D_CHANGESTAT_FN (multi-toggle)
    emptynwstats = numeric(length(k))
  )
}