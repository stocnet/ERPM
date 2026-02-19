# ==============================================================================
# File    : R/InitErgmTerm.cov_match.R
# Term    : cov_match(cov, clique_size = 2, category = NULL,
#                     normalized = c("none","by_group","global"))
# Project : ERPM / ERGM extensions
# ==============================================================================
# Statistic (informal summary):
#   Non-normalized:
#     S_k(B; c)        = sum_g sum_r C(n_{g,r}, k)
#   Targeted category (category = κ):
#     S_k^{(κ)}(B; c)  = sum_g C(n_{g,κ}, k)
#   by_group:
#     sum_g [ S_k(g) / C(n_g, k) ]   (or C(n_{g,κ}, k) / C(n_g, k) when targeted)
#   global:
#     sum_g [ S_k(g) / n_g ]         (or C(n_{g,κ}, k) / n_g when targeted),
#   with groups of size n_g = 0 contributing 0.
# ==============================================================================
#
# IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
# - The compiled changestat is now a D_ entrypoint: D_CHANGESTAT_FN(d_cov_match).
# - Therefore this initializer MUST return `d_func = TRUE`.
# - If you forget d_func=TRUE, ergm will call it as a one-toggle C_ function,
#   leading to signature mismatch and typically a segfault.
#
# Debug initializer output:
#   options(erpm.debug.cov_match_init = TRUE)
#
# Debug C output:
#   set DEBUG_COV_MATCH=1 in src/changestat_cov_match.c and recompile.
# ==============================================================================

#' ERGM term: cov_match (monochromatic cliques by actor covariate)
#'
#' @name InitErgmTerm.cov_match
#' @aliases cov_match
#' @note InitErgmTerm.cov_match.R
#'
#' @description
#' \code{cov_match} is an ERGM term for bipartite networks that counts
#' monochromatic cliques of actors within each group, based on a categorical
#' actor-level covariate.
#'
#' The term is vectorized in \code{clique_size}: several values of \eqn{k} can
#' be specified, yielding one statistic per \eqn{k}.
#'
#' @details
#' See the long description already present in your version; the key change is:
#'
#' IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
#' - The C implementation is multi-toggle (D_ entrypoint).
#' - We MUST advertise this to ergm by returning `d_func = TRUE`.
#'
#' @export
InitErgmTerm.cov_match <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_match"

  # -------------------------------------------------------------------------
  # Optional debug flag for this initializer
  #   options(erpm.debug.cov_match_init = TRUE)
  # -------------------------------------------------------------------------
  DEBUG <- isTRUE(getOption("erpm.debug.cov_match_init", TRUE))

  # -------------------------------------------------------------------------
  # Base ERGM term validation and argument parsing
  # -------------------------------------------------------------------------
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov","clique_size","category","normalized"),
    vartypes      = c("numeric,character","numeric","character","logical,character"),
    defaultvalues = list(NULL,            2,            NULL,       "none"),
    required      = c(TRUE,               FALSE,        FALSE,      FALSE)
  )

  # -------------------------------------------------------------------------
  # 0) Bipartite guard and actor-mode size
  # -------------------------------------------------------------------------
  n1 <- tryCatch(as.integer(nw %n% "bipartite"), error = function(e) NA_integer_)
  if (!is.finite(n1) || n1 <= 0L) {
    ergm_Init_stop(sQuote(termname), ": non-bipartite network or missing/invalid %n% 'bipartite' attribute.")
  }

  # -------------------------------------------------------------------------
  # Debug: verify ERPM wrapper metadata are attached (pure diagnostics)
  # -------------------------------------------------------------------------
  if (DEBUG) {
    .safe_get_n_attr <- function(nw, key) {
      tryCatch(nw %n% key, error = function(e) NULL)
    }
    nodes_meta <- .safe_get_n_attr(nw, "nodes")
    dyads_meta <- .safe_get_n_attr(nw, "dyads")

    cat(sprintf("[Init:%s] debug=TRUE\n", termname))
    cat(sprintf("[Init:%s] network size=%d | bipartite n1=%d\n",
                termname, network::network.size(nw), n1))
    cat(sprintf("[Init:%s] nw %%n%% \"nodes\" : %s\n", termname,
                if (is.null(nodes_meta)) "ABSENT" else paste0("PRESENT (", paste(class(nodes_meta), collapse="/"), ")")))
    cat(sprintf("[Init:%s] nw %%n%% \"dyads\" : %s\n", termname,
                if (is.null(dyads_meta)) "ABSENT" else paste0("PRESENT (", paste(class(dyads_meta), collapse="/"), ")")))
  }

  # -------------------------------------------------------------------------
  # 1) Normalize the 'normalized' argument to an internal mode flag
  # -------------------------------------------------------------------------
  normalized <- a$normalized
  if (is.logical(normalized)) {
    normalized <- if (isTRUE(normalized)) "by_group" else "none"
  }
  normalized <- match.arg(tolower(as.character(normalized)), c("none","by_group","global"))
  norm_mode  <- switch(normalized, none = 0L, by_group = 1L, global = 2L)

  # -------------------------------------------------------------------------
  # 2) Normalize clique sizes k (clique_size)
  # -------------------------------------------------------------------------
  ks <- as.integer(round(a$clique_size))
  if (length(ks) < 1L || any(!is.finite(ks)) || any(ks < 1L)) {
    ergm_Init_stop(sQuote(termname), ": 'clique_size' must contain finite integers >= 1.")
  }

  .allow_k1_nn <- isTRUE(getOption("ERPM.allow.k1.nonnormalized", FALSE))
  if (any(ks == 1L) && (normalized %in% c("none","global")) && !.allow_k1_nn) {
    ergm_Init_stop(
      sQuote(termname),
      ": cov_match(..., clique_size=1) with normalized='", normalized,
      "' is constant in the partition setting. ",
      "Use k>=2, normalized='by_group', or offset(...). ",
      "To force: options(ERPM.allow.k1.nonnormalized=TRUE)."
    )
  }

  ks <- sort(unique(ks))
  K  <- length(ks)

  # -------------------------------------------------------------------------
  # 3) Build actor-level category codes (z) and targeted category info
  # -------------------------------------------------------------------------
  cov      <- a$cov
  category <- a$category

  get_actor_codes <- function(nw, cov, category = NULL, n1, termname, DEBUG = FALSE) {
    ia <- seq_len(n1)

    if (is.character(cov) && length(cov) == 1L) {
      vals <- network::get.vertex.attribute(nw, cov)
      if (is.null(vals)) {
        ergm_Init_stop(sQuote(termname), ": missing vertex attribute: ", sQuote(cov), ".")
      }

      x <- vals[ia]
      if (is.numeric(x)) {
        ergm_Init_stop(sQuote(termname), ": 'cov_match' requires a categorical covariate (factor/character), not numeric.")
      }

      f <- as.factor(x)
      if (!is.null(category)) {
        category <- as.character(category)[1L]
        if (!(category %in% levels(f))) levels(f) <- c(levels(f), category)
      }

      z <- as.integer(f)
      z[is.na(z)] <- 0L
      levs <- levels(f)

      kappa_code <- 0L
      cov_label  <- cov
      if (!is.null(category)) {
        kappa_code <- as.integer(match(category, levs))
        cov_label  <- paste0(cov, "==", category)
      }

      if (DEBUG) {
        cat(sprintf("[Init:%s] cov attribute=%s | actor N=%d | NA actors=%d\n",
                    termname, cov, n1, sum(z == 0L)))
        cat(sprintf("[Init:%s] levels (%d): %s\n",
                    termname, length(levs), paste(levs, collapse = ", ")))
        if (!is.null(category)) {
          cat(sprintf("[Init:%s] targeted category=%s | kappa_code=%d\n",
                      termname, category, kappa_code))
        }
      }

      return(list(
        z          = as.double(z),
        kappa_code = as.double(kappa_code),
        cov_label  = cov_label,
        levels     = levs
      ))
    }

    if (is.numeric(cov)) {
      ergm_Init_stop(sQuote(termname), ": 'cov_match' requires a categorical covariate (factor/character), not a numeric vector.")
    }
    if (length(cov) < n1) {
      ergm_Init_stop(sQuote(termname), ": length(cov) < |A| = ", n1, ".")
    }

    x <- cov[ia]
    f <- as.factor(x)

    if (!is.null(category)) {
      category <- as.character(category)[1L]
      if (!(category %in% levels(f))) levels(f) <- c(levels(f), category)
    }

    z <- as.integer(f)
    z[is.na(z)] <- 0L
    levs <- levels(f)

    kappa_code <- 0L
    cov_label  <- "cov"
    if (!is.null(category)) {
      kappa_code <- as.integer(match(category, levs))
      cov_label  <- paste0("cov==", category)
    }

    if (DEBUG) {
      cat(sprintf("[Init:%s] cov vector | actor N=%d | NA actors=%d\n",
                  termname, n1, sum(z == 0L)))
      cat(sprintf("[Init:%s] levels (%d): %s\n",
                  termname, length(levs), paste(levs, collapse = ", ")))
      if (!is.null(category)) {
        cat(sprintf("[Init:%s] targeted category=%s | kappa_code=%d\n",
                    termname, category, kappa_code))
      }
    }

    return(list(
      z          = as.double(z),
      kappa_code = as.double(kappa_code),
      cov_label  = cov_label,
      levels     = levs
    ))
  }

  ax <- get_actor_codes(nw, cov = cov, category = category, n1 = n1, termname = termname, DEBUG = DEBUG)
  z_codes    <- ax$z
  kappa_code <- ax$kappa_code
  cov_label  <- ax$cov_label

  has_kappa <- as.double(as.integer(kappa_code > 0))

  if (DEBUG) {
    cat(sprintf("[Init:%s] normalized=%s (mode=%d) | K=%d | ks={%s}\n",
                termname, normalized, norm_mode, K, paste(ks, collapse=",")))
  }

  # -------------------------------------------------------------------------
  # 4) Build INPUT_PARAM vector for the C change-statistic (D_ entrypoint)
  # -------------------------------------------------------------------------
  inputs <- c(
    as.double(n1),
    as.double(K),
    as.double(norm_mode),
    has_kappa,
    as.double(kappa_code),
    as.double(ks),
    z_codes
  )

  # -------------------------------------------------------------------------
  # 5) Coefficient names
  # -------------------------------------------------------------------------
  suffix_norm <- switch(normalized,
                        none     = "",
                        by_group = "_bygrp",
                        global   = "_glob")
  coef.names  <- paste0("cov_match[", cov_label, "]_k", ks, suffix_norm)

  # -------------------------------------------------------------------------
  # 6) Standard ERGM term specification
  # -------------------------------------------------------------------------
  # IMPORTANT:
  # - name stays "cov_match" (term name in formulas).
  # - d_func=TRUE tells ergm to call D_CHANGESTAT_FN(d_cov_match).
  list(
    name         = "cov_match",
    coef.names   = coef.names,
    inputs       = inputs,
    dependence   = TRUE,
    d_func       = TRUE,          # <-- REQUIRED (multi-toggle C entrypoint)
    emptynwstats = rep(0, K)
  )
}