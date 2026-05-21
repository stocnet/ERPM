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
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' \code{dyadcov} is an ERGM term for bipartite networks that aggregates a
#' numeric dyadic covariate \eqn{Z = (z_{ij})} over cliques of actors within
#' each group.
#'
#' The network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side of the bipartite graph).
#' }
#'
#' For a group g, let A(g) be the set of adjacent actors. For a fixed
#' \eqn{k \ge 2}, define:
#' \deqn{
#'   S_g^{(k)}(Z)
#'     = \sum_{C \in C_k(g)} \prod_{i<j,\, i,j \in C} (z_{ij} + z_{ji}).
#' }
#'
#' Three variants are supported, controlled by \code{normalize} (aliases
#' \code{normalized} and \code{norm} for backward compatibility):
#' \itemize{
#'   \item raw sum (no normalisation);
#'   \item global normalisation \code{"global"}: factor \eqn{1 / n_g} per group;
#'   \item by-group normalisation \code{"by_group"}: factor
#'         \eqn{1 / \binom{n_g}{k}} per group.
#' }
#'
#' @details
#' The dyadic covariate \eqn{Z} is defined on the actor mode (size \eqn{n1})
#' and is read in column-major order (R convention). Z may be non-symmetric:
#' for each unordered pair \{i,j\}, the term uses (z_ij + z_ji).
#'
#' IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
#' - This term must be safe under multi-toggle MCMC proposals (lists of toggles).
#' - The compiled change-statistic is implemented using the D_CHANGESTAT_FN API
#'   (multi-toggle) as `d_dyadcov`.
#' - Therefore, the R initializer MUST return `d_func = TRUE`.
#'   If you forget this, ergm will call the function using the one-toggle
#'   signature, which is a signature mismatch and can crash R.
#'
#' @template erpm-initergmterm-args
#'
#' @export
InitErgmTerm.dyadcov <- function(nw, arglist, ..., version = packageVersion("ergm")) {
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
  # Actor-mode size (n1) from bipartite attribute
  # ---------------------------------------------------------------------------
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L)
    stop(termname, ": strictly bipartite network required (attribut %n% 'bipartite' manquant ou invalide).")
  dbgcat("n1 = ", n1)

  # ---------------------------------------------------------------------------
  # Retrieve the dyadic matrix (n1 x n1)
  # ---------------------------------------------------------------------------
  dyad_raw   <- a$dyadcov
  dyad_label <- NULL

  if (is.character(dyad_raw) && length(dyad_raw) == 1L) {
    dyad_name <- dyad_raw

    dyad_mat <- nw %n% dyad_name
    src <- "network attribute"

    if (is.null(dyad_mat)) {
      dyads_list <- nw %n% "dyads"
      if (is.list(dyads_list) && !is.null(dyads_list[[dyad_name]])) {
        dyad_mat <- dyads_list[[dyad_name]]
        src <- "nw %n% 'dyads' list"
      }
    }

    if (is.null(dyad_mat)) {
      stop(
        termname, ": dyad not found : ", sQuote(dyad_name),
        " (looked up in nw %n% ", sQuote(dyad_name),
        " and in nw %n% 'dyads'[[", sQuote(dyad_name), "]])."
      )
    }

    dyad_label <- if (identical(src, "nw %n% 'dyads' list")) paste0("dyads$", dyad_name) else dyad_name
    dbgcat("dyadcov source = ", src, " ", sQuote(dyad_label))
  } else {
    dyad_mat   <- dyad_raw
    dyad_label <- "dyadcov"
    dbgcat("dyadcov source = literal matrix")
  }

  if (!is.matrix(dyad_mat))
    stop(termname, ": 'dyadcov' must be a matrix or the name of a network-level attribute.")

  nr <- nrow(dyad_mat)
  nc <- ncol(dyad_mat)

  if (nr < n1 || nc < n1)
    stop(termname, ": dyadic matrix dimensions (", nr, "x", nc,
         ") insufficient for n1 = ", n1, ".")

  if (nr > n1 || nc > n1) {
    dyad_mat <- dyad_mat[seq_len(n1), seq_len(n1), drop = FALSE]
    dbgcat("dyadcov truncated to ", n1, "x", n1)
  }

  if (!is.numeric(dyad_mat))
    stop(termname, ": dyadic matrix should be numeric.")
  if (anyNA(dyad_mat))
    stop(termname, ": NA values are not allowed in the dyadic matrix.")

  # Optional symmetry diagnostics (non-blocking)
  tol <- 1e-8
  max_asym <- max(abs(dyad_mat - t(dyad_mat)))
  if (dbg && is.finite(max_asym) && max_asym > tol) {
    dbgcat("warning : dyadcov matrix is not symmetric, max |Z - t(Z)| = ",
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
      stop(termname, ": 'clique_size' must be numeric or integer.")
    k <- as.integer(round(k_raw[1L]))
  }
  if (!is.finite(k) || k < 2L)
    stop(termname, ": 'clique_size' must be an integer >= 2.")
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
    if (isTRUE(norm_raw[1L])) {
      norm_mode  <- 1L
      norm_label <- "global"
    } else {
      norm_mode  <- 0L
      norm_label <- "none"
    }
  } else if (is.character(norm_raw)) {
    val <- match.arg(norm_raw[1L], c("none",
                                    "global", "by_group",
                                    "size", "cliques"))
    if (val == "none") {
      norm_mode  <- 0L
      norm_label <- "none"
    } else if (val %in% c("global", "size")) {
      norm_mode  <- 1L
      norm_label <- "global"
    } else {
      norm_mode  <- 2L
      norm_label <- "by_group"
    }
  } else {
    stop(termname, ": 'normalize' must be logical or character. ",
         "(\"none\", \"global\", \"by_group\").")
  }

  dbgcat("normalized mode = ", norm_label, " (code=", norm_mode, ")")

  # ---------------------------------------------------------------------------
  # Coefficient name
  # ---------------------------------------------------------------------------
  base_label <- sprintf("dyadcov[%s]_k%d", dyad_label, k)
  coef.name  <- switch(
    as.character(norm_mode),
    "0" = base_label,
    "1" = paste0(base_label, "_global"),
    "2" = paste0(base_label, "_bygrp"),
    base_label
  )
  dbgcat("coef.name = ", coef.name)

  # ---------------------------------------------------------------------------
  # Build INPUT_PARAM for the C side (unchanged layout)
  # ---------------------------------------------------------------------------
  inputs <- c(
    as.double(n1),
    as.double(k),
    as.double(norm_mode),
    as.double(dyad_mat)  # column-major
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " k=", k, " norm_mode=", norm_mode,
         " | Z[1:6]=", paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L),
                             collapse = ","))

  # ---------------------------------------------------------------------------
  # Standard ERGM term initialization return value
  # ---------------------------------------------------------------------------
  # IMPORTANT: d_func = TRUE is REQUIRED because the compiled changestat is D_.
  list(
    name         = "dyadcov",
    coef.names   = coef.name,
    inputs       = inputs,
    dependence   = TRUE,
    d_func       = TRUE,   # <-- REQUIRED for D_CHANGESTAT_FN (multi-toggle)
    minval       = -Inf,
    maxval       = Inf,
    emptynwstats = 0
  )
}