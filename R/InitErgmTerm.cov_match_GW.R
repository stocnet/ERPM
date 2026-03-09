# ==============================================================================
# File    : R/InitErgmTerm.cov_match_GW.R
# Term    : cov_match_GW(cov, lambda = 2, category = NULL,
#                        normalized = c("none","by_group","global"))
# Project : ERPM / ERGM extensions
# ==============================================================================
# Statistic (informal summary)
# ------------------------------------------------------------------------------
# Non-normalized:
#   S_GW(B; c, l)       = sum_g sum_r l * (1 - r_l^{ n_{g,r} })
# Targeted category (category = k):
#   S_GW^{(k)}(B; c, l) = sum_g l * (1 - r_l^{ n_{g,k} })
# By-group:
#   sum_g [ Num(g) / Den(g) ] with
#     Num(g) = sum_r l(1 - r_l^{n_{g,r}})    (or l(1 - r_l^{n_{g,k}}) when targeted)
#     Den(g) = l(1 - r_l^{n_g})
# Global:
#   [ sum_g Num(g) ] / [ l(1 - r_l^{N_A}) ]
# where:
#   - N_A is the number of actors in the actor mode,
#   - n_{g,r} is the number of actors of category r in group g,
#   - r_l = (l - 1) / l with l > 1.
#
# IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
# ------------------------------------------------------------------------------
# This term MUST support multi-toggle moves (swap/split/merge decomposed into
# several membership toggles).
#
# Therefore:
#   - The compiled changestat MUST be implemented as:
#       D_CHANGESTAT_FN(d_cov_match_GW)
#   - On the R side, we MUST advertise this to ergm by returning:
#       d_func = TRUE
#
# If you forget `d_func = TRUE`, ergm will try to call the changestat as a
# one-toggle C_CHANGESTAT_FN with the wrong signature (=> crash/segfault).
#
# Symbol naming convention:
#   - Implement ONLY `d_cov_match_GW` (D-signature).
#   - Avoid exporting a `c_cov_match_GW` symbol with a D-signature: ergm may
#     resolve it as the one-toggle entrypoint and crash.
# ==============================================================================

#' ERGM term: cov_match_GW (geometrically weighted monochromatic cliques)
#' @name InitErgmTerm.cov_match_GW
#' @aliases cov_match_GW
#' @note InitErgmTerm.cov_match_GW.R
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' \code{cov_match_GW} is an ERGM term for bipartite networks that applies a
#' geometrically weighted transform to the monochromatic clique counts produced
#' by \code{cov_match}, based on a categorical actor-level covariate.
#'
#' The network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side of the bipartite graph).
#' }
#'
#' For each group in the group mode and each category of the actor covariate, the
#' term replaces raw counts of monochromatic cliques by a geometrically weighted
#' function driven by \eqn{\lambda > 1}, with ratio
#' \deqn{
#'   r_\lambda = \frac{\lambda - 1}{\lambda} \in (0, 1).
#' }
#'
#' Several normalization variants are supported (\code{"none"}, \code{"by_group"},
#' \code{"global"}), as well as an optional targeted category.
#'
#' @details
#' (Documentation inchangée : voir version précédente du fichier.)
#'
#' @keywords ERGM term bipartite categorical covariate geometrically weighted
#' @md
#'
#' @export
InitErgmTerm.cov_match_GW <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_match_GW"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  # Global option:
  #   options(ERPM.cov_match_GW.debug = TRUE/FALSE)
  # When TRUE, the initializer prints diagnostic messages to the console.
  dbg    <- isTRUE(getOption("ERPM.cov_match_GW.debug", TRUE))
  dbgcat <- function(...) if (dbg) cat("[cov_match_GW][DEBUG]", ..., "\n", sep = "")

  dbgcat("InitErgmTerm.cov_match_GW called with args: ",
         paste(names(arglist), collapse = ", "))

  # ---------------------------------------------------------------------------
  # Base ERGM term validation and argument parsing
  # ---------------------------------------------------------------------------
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov","lambda","category","normalized"),
    vartypes      = c("character","numeric","character","logical,character"),
    defaultvalues = list(NULL,        2,        NULL,       "none"),
    required      = c(TRUE,           FALSE,    FALSE,      FALSE)
  )

  # ---------------------------------------------------------------------------
  # Strict bipartite guard and actor-mode size
  #   - n1 is the size of the actor mode
  # ---------------------------------------------------------------------------
  n1 <- tryCatch(nw %n% "bipartite", error = function(e) NA_integer_)
  if (!is.numeric(n1) || !is.finite(n1) || n1 <= 0) {
    ergm_Init_stop(
      sQuote(termname),
      ": non-bipartite network or missing/invalid %n% 'bipartite' attribute."
    )
  }
  n1 <- as.integer(n1)
  dbgcat("bipartite attribute (n1) =", n1)

  # ---------------------------------------------------------------------------
  # Normalization mode: map "none" / "by_group" / "global" to an integer flag
  # ---------------------------------------------------------------------------
  normalized <- a$normalized
  if (is.logical(normalized)) normalized <- if (isTRUE(normalized)) "by_group" else "none"
  normalized <- match.arg(tolower(as.character(normalized)), c("none","by_group","global"))
  norm_mode  <- switch(normalized, none = 0L, by_group = 1L, global = 2L)
  dbgcat("normalized =", normalized, " -> norm_mode =", norm_mode)

  # ---------------------------------------------------------------------------
  # Lambda handling:
  #   - default lambda = 2
  #   - require all lambda > 1 and finite
  #   - deduplicate to get K distinct values
  # ---------------------------------------------------------------------------
  lambdas <- as.double(a$lambda)
  if (!length(lambdas)) lambdas <- 2
  if (any(!is.finite(lambdas)) || any(lambdas <= 1)) {
    ergm_Init_stop(sQuote(termname), ": 'lambda' must be > 1 (numeric, finite).")
  }
  lambdas <- as.double(unique(lambdas))
  K <- length(lambdas)
  dbgcat("lambda (unique) =", paste(format(lambdas), collapse = ", "), " | K =", K)

  # ---------------------------------------------------------------------------
  # Categorical actor attribute and targeted category handling
  #   - covname identifies an actor-level attribute
  #   - transform to factor and encode as integer codes 1..R (0 for NA)
  #   - category, if provided, is ensured to be among the levels
  # ---------------------------------------------------------------------------
  covname <- a$cov
  if (!(is.character(covname) && length(covname) == 1L)) {
    ergm_Init_stop(sQuote(termname), ": 'cov' must be the name of an actor attribute (factor/character).")
  }

  # Indices for the actor mode (here simply 1..n1)
  ia <- seq_len(n1)

  # Retrieve the actor-level covariate values
  vals <- network::get.vertex.attribute(nw, covname)
  if (is.null(vals)) {
    ergm_Init_stop(sQuote(termname), ": nonexistent attribute : ", sQuote(covname), ".")
  }

  # Coerce to factor and restrict to actors
  f <- as.factor(vals[ia])

  # Targeted category handling:
  #   - if category is not NULL and not in levels, extend levels so that the
  #     targeted category exists with zero frequency
  category <- a$category
  if (!is.null(category) && !(category %in% levels(f))) {
    dbgcat("category not in levels -> extending levels with: ", category)
    levels(f) <- c(levels(f), category)
  }

  # Encode categorical values to integer codes:
  #   - 1..R for valid levels
  #   - NA mapped to 0 (ignored in C code)
  z <- as.integer(f)
  z[!is.finite(z)] <- 0L

  # Targeted category code:
  #   - 0 if no category is targeted
  #   - otherwise, the integer level index of the targeted category
  kappa_code <- if (is.null(category)) 0L else as.integer(match(category, levels(f)))
  has_kappa  <- as.integer(kappa_code > 0L)

  # Label for coefficient names:
  cov_label  <- if (is.null(category)) covname else paste0(covname, "==", category)

  dbgcat("cov =", covname,
         " | levels =", paste(levels(f), collapse = ","),
         " | has_kappa =", has_kappa,
         " | kappa_code =", kappa_code)

  # ---------------------------------------------------------------------------
  # Pack INPUT_PARAM for the C change-statistic (multi-toggle)
  #   Layout:
  #     [0] = n1
  #     [1] = K
  #     [2] = norm_mode
  #     [3] = has_kappa
  #     [4] = kappa_code
  #     [5 .. 5+K-1]      = lambdas
  #     [5+K .. 5+K+n1-1] = z[1..n1]
  # ---------------------------------------------------------------------------
  inputs <- c(
    as.double(n1),
    as.double(K),
    as.double(norm_mode),
    as.double(has_kappa),
    as.double(kappa_code),
    lambdas,
    as.double(z)
  )

  dbgcat("inputs length =", length(inputs), " (= 5 + K + n1)")

  # ---------------------------------------------------------------------------
  # Coefficient names
  # ---------------------------------------------------------------------------
  suffix_norm <- switch(normalized,
                        none     = "",
                        by_group = "_bygrp",
                        global   = "_glob")

  # Format lambda values compactly for coefficient names
  fmt_lambda  <- function(x) sub("\\.?0+$", "", format(x, trim = TRUE))

  coef.names <- paste0(
    "cov_match_GW[", cov_label, "]_l",
    vapply(lambdas, fmt_lambda, ""),
    suffix_norm
  )

  dbgcat("coef.names =", paste(coef.names, collapse = " | "))

  # ---------------------------------------------------------------------------
  # Standard ERGM term specification
  # ---------------------------------------------------------------------------
  # IMPORTANT:
  # - `d_func = TRUE` tells ergm to call the multi-toggle (D_) changestat entrypoint.
  # - The compiled symbol must be `d_cov_match_GW`.
  list(
    name         = "cov_match_GW",
    coef.names   = coef.names,  # length = K
    inputs       = inputs,
    dependence   = TRUE,
    d_func       = TRUE,        # <-- REQUIRED for D_CHANGESTAT_FN (multi-toggle)
    emptynwstats = rep(0, K)
  )
}