# ==============================================================================
# File    : InitErgmTerm.cov_diff.R
# Purpose : Declare the ERGM term 'cov_diff' for bipartite networks
#           (range over k-actor subsets within groups).
# Project : ERPM / ERGM extensions
# ============================================================================

#' ERGM term: cov_diff (range over k-actor subsets within groups)
#' @name InitErgmTerm.cov_diff
#' @aliases cov_diff
#' @note InitErgmTerm.cov_diff.R
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' \code{cov_diff} is an ERGM term for bipartite networks that measures, for each
#' group, the dispersion of a numeric actor covariate over all \eqn{k}-actor
#' subsets within that group. The bipartite network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side, representing groups).
#' }
#'
#' Each actor in the actor mode carries a numeric covariate value \eqn{x_i}.
#' For a fixed integer \eqn{k \ge 2}, and for each group, the term considers all
#' \eqn{k}-actor subsets of that group and computes, on each subset, the
#' max-min range of the covariate. The statistic can be used either:
#' \itemize{
#'   \item in its raw form (sum over all subsets in all groups);
#'   \item in a by-group normalized form (average range per \eqn{k}-subset in each group);
#'   \item in a global, size-normalized form (average range per actor in each group).
#' }
#' The choice is controlled by the \code{normalize} argument (with aliases
#' \code{normalized} and \code{norm}).
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic.
#'
#' IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
#' - This term MUST support multi-toggle proposals (swap/split/merge represented
#'   as a list of membership edge toggles).
#' - Therefore, the compiled change-statistic is implemented as a D_ entrypoint
#'   (D_CHANGESTAT_FN), typically named `d_cov_diff`.
#' - On the R side, we MUST advertise this to ergm by returning `d_func = TRUE`.
#'   Otherwise ergm will try to call the changestat as a one-toggle C_ function,
#'   causing a signature mismatch and typically a segfault.
#'
#' Compiled symbol naming convention:
#' - Recommended: implement the C function as `d_cov_diff` via
#'   D_CHANGESTAT_FN(d_cov_diff).
#' - Avoid exposing a symbol named `c_cov_diff` with a D-signature, because
#'   ergm may resolve it as the one-toggle entrypoint and crash.
#'
#' The R initializer:
#' \itemize{
#'   \item enforces that the network is bipartite (actor mode / group mode);
#'   \item parses the user arguments \code{cov}, \code{clique_size}, and normalization;
#'   \item validates the actor covariate (numeric coercible, no NA on actor mode);
#'   \item checks \code{clique_size} is a single integer \eqn{\ge 2};
#'   \item normalizes the \code{normalize} argument into \code{norm_mode ∈ \{0,1,2\}};
#'   \item builds a compact \code{inputs} vector for the C code encoding:
#'         \code{n1}, \code{k}, \code{norm_mode}, and the actor covariate \code{x}.
#' }
#'
#' @template erpm-initergmterm-args
#'
#' @export
InitErgmTerm.cov_diff <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_diff"

  # Debug helpers:
  # - dbg: logical flag, controlled by an R option;
  # - dbgcat(): emits prefixed debug messages when dbg is TRUE.
  dbg    <- isTRUE(getOption("ERPM.cov_diff.debug", TRUE))
  dbgcat <- function(...) if (dbg) cat("[cov_diff][DEBUG]", ..., "\n", sep = "")

  # Run standard ERGM term checks and parse user arguments:
  # - enforce bipartite network;
  # - accept 'cov', 'clique_size', and a normalization argument
  #   ('normalize', 'normalized', or 'norm') with flexible types;
  # - let \pkg{ergm} handle generic validations (missing args, etc.).
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov",
                      "clique_size",
                      "normalize",
                      "normalized",
                      "norm"),
    vartypes      = c("character,numeric,logical,vector",
                      "numeric,integer",
                      "logical,character,numeric",
                      "logical,character,numeric",
                      "logical,character,numeric"),
    defaultvalues = list(NULL,
                         2,
                         NULL,
                         NULL,
                         NULL),
    required      = c(TRUE,
                      FALSE,
                      FALSE,
                      FALSE,
                      FALSE)
  )

  # ----- 1) Actor-mode size n1 -----------------------------------------------
  # n1 is the number of actors, retrieved from the bipartite network attribute.
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L) stop(termname, ": strictly non-bipartite network.")
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
    stop(termname, ": length of the covariate < n1.")

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
    stop(termname, ": the covariate must be coercible to numeric.")

  if (anyNA(cov_vec))
    stop(termname, ": NA values are not allowed in the actor-mode covariate.")

  # ----- 4) Subset size 'clique_size' = k >= 2 -------------------------------
  # 'clique_size' is interpreted as the subset size k used in the definition
  # of the statistic. It must be a single finite numeric value, rounded to
  # an integer and required to be at least 2.
  k_raw <- a$clique_size
  if (length(k_raw) != 1L || !is.numeric(k_raw))
    stop(termname, ": 'clique_size' must be a numeric scalar.")
  k <- as.integer(round(k_raw))
  if (!is.finite(k) || k < 2L)
    stop(termname, ": 'clique_size' must be an integer >= 2.")
  dbgcat("clique_size k = ", k)

  # ----- 5) Normalization mode (raw / by_group / global) ---------------------
  # The normalization argument can be provided under the names:
  # - normalize
  # - normalized
  # - norm
  #
  # Its semantics:
  # - NULL / FALSE (or missing): raw statistic (no normalisation);
  # - TRUE or "by_group"       : per-group average over k-subsets;
  # - "global"                 : per-group average per actor;
  # - numeric 0,1,2            : raw, by-group, global.
  norm_raw <- a$normalize
  if (is.null(norm_raw)) norm_raw <- a$normalized
  if (is.null(norm_raw)) norm_raw <- a$norm

  if (is.null(norm_raw)) {
    normalized <- "none"
  } else if (is.logical(norm_raw)) {
    normalized <- if (isTRUE(norm_raw)) "by_group" else "none"
  } else if (is.character(norm_raw) && length(norm_raw) == 1L) {
    normalized <- match.arg(tolower(norm_raw), c("by_group", "global"))
  } else if (is.numeric(norm_raw) && length(norm_raw) == 1L) {
    if (norm_raw == 0) {
      normalized <- "none"
    } else if (norm_raw == 1) {
      normalized <- "by_group"
    } else if (norm_raw == 2) {
      normalized <- "global"
    } else {
      stop(termname, ": numeric 'normalize' must be 0 (raw), 1 ('by_group'), or 2 ('global').")
    }
  } else {
    stop(termname, ": 'normalize'/'normalized'/'norm' must be numeric or 'by_group'/'global'.")
  }

  norm_mode <- switch(normalized,
                      none     = 0L,
                      by_group = 1L,
                      global   = 2L)
  norm_label <- normalized

  dbgcat("normalized = ", norm_label, " (mode=", norm_mode, ")")

  # ----- 6) Coefficient name --------------------------------------------------
  # The coefficient name encodes the covariate label, the subset size k and
  # the normalization mode, for interpretability in model summaries.
  suffix_norm <- switch(norm_label,
                        none     = "",
                        by_group = "_bygrp",
                        global   = "_glob")
  coef.name <- sprintf(
    "cov_diff[%s]_k%d%s",
    cov_label, k, suffix_norm
  )
  dbgcat("coef.name = ", coef.name)

  # ----- 7) Build INPUT_PARAM for the C layer --------------------------------
  # Layout of INPUT_PARAM:
  #   [1]   = n1            (actor-mode size)
  #   [2]   = k             (subset size)
  #   [3]   = norm_mode     (0 raw, 1 by-group, 2 global)
  #   [4..] = cov_vec[1:n1] (numeric covariate on the actor mode)
  inputs <- c(
    as.double(n1),
    as.double(k),
    as.double(norm_mode),
    as.double(cov_vec)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " k=", k, " norm_mode=", norm_mode,
         " | cov[1:6]=", paste(utils::head(signif(cov_vec, 5L), 6L), collapse = ","))

  # ----- 8) Return ERGM term specification -----------------------------------
  # IMPORTANT:
  # - `d_func = TRUE` tells ergm to call the multi-toggle (D_) changestat entrypoint.
  # - Without it, ergm assumes a one-toggle C_ changestat and may call with the
  #   wrong signature (=> segfault) if only a D_ function is compiled.
  list(
    name         = "cov_diff",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, k, norm_mode, x[n1]
    dependence   = TRUE,
    d_func       = TRUE,        # <-- REQUIRED for D_CHANGESTAT_FN (multi-toggle)
    minval       = 0,
    maxval       = Inf,
    emptynwstats = 0
  )
}