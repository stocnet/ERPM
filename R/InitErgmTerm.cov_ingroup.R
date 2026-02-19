# ==============================================================================
# File    : R/InitErgmTerm.cov_ingroup.R
# Term    : cov_ingroup(cov, size = NULL, category = NULL)
# Project : ERPM / ERGM extensions
# ==============================================================================
# Informal definition of the statistic
# ------------------------------------------------------------------------------
# Let:
#   - A be the actor mode (|A| = n1 = nw %n% "bipartite");
#   - G be the group mode (complementary side of the bipartite graph);
#   - B be the actor-group incidence matrix;
#   - x_i be a numeric value attached to actor i in A;
#   - n_g be the size of group g (number of adjacent actors);
#   - S be a set of admissible group sizes (encoded by `size`).
#
# The statistic is:
#   T(B; x, S) = sum_g [ n_g * (sum_{i in g} x_i) * 1[n_g in S] ].
#
# When `category` is provided and `cov` is categorical, x_i is replaced by
# an indicator 1[c_i == category] so that the term becomes a weighted count
# of actors in the targeted category inside each group, scaled by group size.
#
# ------------------------------------------------------------------------------
# IMPORTANT (multi-toggle / D_CHANGESTAT_FN)
# ------------------------------------------------------------------------------
# This term MUST support multi-toggle proposals (swap/split/merge) represented
# by a list of toggles inside ergm's MCMC.
#
# Therefore:
#   - The compiled change-statistic MUST be implemented with D_CHANGESTAT_FN
#     (multi-toggle signature).
#   - On the R side, we MUST return d_func = TRUE so that ergm calls the D_
#     entrypoint (and does NOT attempt to call a one-toggle C_ entrypoint).
#
# Compiled symbol naming convention:
#   - Implement the C function as `d_cov_ingroup` via D_CHANGESTAT_FN(d_cov_ingroup).
#   - Do NOT expose a `c_cov_ingroup` symbol with a D signature (ergm may resolve
#     it as one-toggle and crash).
# ==============================================================================

#' ERGM term: cov_ingroup (group-size weighted covariate sums)
#' @name InitErgmTerm.cov_ingroup
#' @aliases cov_ingroup
#' @note InitErgmTerm.cov_ingroup.R
#'
#' @description
#' \code{cov_ingroup} is an ERGM term for bipartite networks that aggregates
#' actor-level covariates within groups, with an optional filter on group sizes.
#' The network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side of the bipartite graph).
#' }
#'
#' For each group in the group mode, let \eqn{n_g} be the group size (number of
#' adjacent actors) and \eqn{x_i} a numeric covariate value attached to each
#' actor \eqn{i} in the actor mode. The term computes a size-weighted sum of
#' within-group covariate totals:
#' \deqn{
#'   T(B; x, S) = \sum_{g \in G} n_g \left(\sum_{i \in g} x_i\right)\,
#'   \mathbf{1}[n_g \in S],
#' }
#' where \eqn{S} is the set of admissible group sizes derived from the
#' \code{size} argument. If \code{size} is \code{NULL}, all group sizes are
#' included and \eqn{\mathbf{1}[n_g \in S] \equiv 1}.
#'
#' When \code{category} is supplied and \code{cov} is categorical, the actor
#' covariate is replaced by the indicator \eqn{x_i = 1[c_i = \kappa]}, where
#' \eqn{\kappa} is the targeted category. In that case, the term becomes:
#' \deqn{
#'   T(B; \kappa, S) =
#'   \sum_{g \in G} n_g \left(\sum_{i \in g} \mathbf{1}[c_i = \kappa]\right)
#'   \mathbf{1}[n_g \in S],
#' }
#' i.e. a group-size-weighted measure of how many actors with category
#' \eqn{\kappa} are present inside each group.
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic under the name
#' \code{d_cov_ingroup} (multi-toggle form). The R initializer:
#' \itemize{
#'   \item enforces that the network is bipartite and retrieves the actor-mode
#'         size from \code{nw \%n\% "bipartite"};
#'   \item builds a numeric vector \eqn{x \in \mathbb{R}^{n_1}} of covariate
#'         values on actors, either from a vertex attribute or a direct vector;
#'   \item normalizes and validates the size filter \code{size}, converting it
#'         into a sorted set \eqn{S} of positive integers;
#'   \item packs \eqn{n_1}, the size filter and \eqn{x} into a compact
#'         \code{INPUT_PARAM} vector consumed by the C code.
#'   \item sets \code{d_func = TRUE} to advertise multi-toggle support to ergm.
#' }
#'
#' On each toggle list (possibly with several toggles affecting the same group),
#' the C change-statistic processes toggles sequentially, temporarily applying
#' intermediate toggles so degrees and neighbour sums are consistent, then
#' undoes them before returning.
#'
#' @section INPUT_PARAM layout (C side):
#' The numeric vector passed to \code{d_cov_ingroup} has the following layout:
#'
#' \preformatted{
#'   INPUT_PARAM = c(
#'     n1,          # actor-mode size |A|
#'     L,           # number of sizes in S (length(size))
#'     sizes[1:L],  # admissible group sizes S (possibly L = 0 => all sizes)
#'     x[1:n1]      # covariate values on actors (numeric)
#'   )
#' }
#'
#' @export
InitErgmTerm.cov_ingroup <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_ingroup"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  # Global option:
  #   options(ERPM.cov_ingroup.debug = TRUE/FALSE)
  dbg    <- isTRUE(getOption("ERPM.cov_ingroup.debug", TRUE))
  dbgcat <- function(...) if (dbg) cat("[cov_ingroup][DEBUG]", ..., "\n", sep = "")

  dbgcat("InitErgmTerm.cov_ingroup called with args: ",
         paste(names(arglist), collapse = ", "))

  # ---------------------------------------------------------------------------
  # Base ERGM term validation and argument parsing
  #   - enforce bipartite network
  #   - accept `cov`, optional `size`, optional `category`
  # ---------------------------------------------------------------------------
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov",               "size",   "category"),
    vartypes      = c("numeric,character", "numeric","character"),
    defaultvalues = list(NULL,             NULL,     NULL),
    required      = c(TRUE,                FALSE,    FALSE)
  )

  # ---------------------------------------------------------------------------
  # Actor-mode size (n1)
  #   - strict bipartite guard via nw %n% "bipartite"
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
  # Build the actor-level covariate vector x (length n1)
  #   - handle both attribute name and direct numeric vector
  #   - handle optional category for categorical attributes
  # ---------------------------------------------------------------------------
  cov      <- a$cov
  category <- a$category

  get_actor_cov <- function(nw, cov, category = NULL) {
    n1 <- as.integer(nw %n% "bipartite")

    # Indices for actor mode: by convention, actors are 1..n1.
    ia <- seq_len(n1)

    vn <- network::network.vertex.names(nw)
    if (length(vn) >= n1) {
      ia_guess <- which(!grepl("^G\\d+$", vn))
      if (length(ia_guess) == n1) ia <- ia_guess
    }

    # Case 1: cov is the name of a vertex attribute
    if (is.character(cov) && length(cov) == 1L) {
      vals <- network::get.vertex.attribute(nw, cov)
      if (is.null(vals)) {
        ergm_Init_stop(sQuote(termname), ": attribut inexistant: ", sQuote(cov), ".")
      }
      x <- vals[ia]

      # Category => indicator x_i = 1[c_i == category]
      if (!is.null(category)) {
        xb <- as.integer(as.character(x) == category)
        xb[is.na(xb)] <- 0L
        return(list(
          x         = as.double(xb),
          cov_label = paste0(cov, "==", category)
        ))
      }

      # Numeric covariate (coerce + strict finite)
      x_num <- suppressWarnings(as.numeric(x))
      if (any(!is.finite(x_num))) {
        bad <- which(!is.finite(x_num))[1]
        ergm_Init_stop(
          sQuote(termname),
          ": numeric actor attribute contains NA/NaN/Inf. ",
          "Example: vertex=", if (length(vn) >= ia[bad]) vn[ia[bad]] else ia[bad],
          ", value=", as.character(x[bad])
        )
      }
      return(list(
        x         = as.double(x_num),
        cov_label = cov
      ))
    }

    # Case 2: direct numeric vector
    x_num <- suppressWarnings(as.numeric(cov))
    if (any(!is.finite(x_num))) {
      ergm_Init_stop(sQuote(termname), ": vector 'cov' contains NA/NaN/Inf.")
    }
    if (length(x_num) < n1) {
      ergm_Init_stop(sQuote(termname), ": length(cov) < |A| = ", n1, ".")
    }
    if (!is.null(category)) {
      ergm_Init_stop(sQuote(termname), ": 'category' does not apply when 'cov' is a direct numeric vector.")
    }

    list(
      x         = as.double(x_num[seq_len(n1)]),
      cov_label = "cov"
    )
  }

  ax <- get_actor_cov(nw, cov = cov, category = category)
  x         <- ax$x
  cov_label <- ax$cov_label

  dbgcat("cov_label =", cov_label, " | x[1:5] = ",
         paste(utils::head(x, 5), collapse = ","))

  # ---------------------------------------------------------------------------
  # Normalize the size filter S (argument 'size')
  #   - NULL or empty => all group sizes (encoded as L = 0)
  #   - otherwise: distinct positive integers, sorted
  # ---------------------------------------------------------------------------
  S <- a$size
  if (is.null(S) || length(S) == 0L) {
    sizes <- integer(0)
  } else {
    S <- unique(as.integer(S))
    if (any(!is.finite(S)) || any(S < 1L)) {
      ergm_Init_stop(sQuote(termname), ": 'size' must contain integers >= 1.")
    }
    sizes <- sort(S)
  }
  L <- length(sizes)

  dbgcat("sizes filter L =", L, " | sizes = ",
         if (L) paste(sizes, collapse = ",") else "<ALL>")

  # ---------------------------------------------------------------------------
  # Build INPUT_PARAM for the C change-statistic
  #   Layout:
  #     [1]           = n1
  #     [2]           = L
  #     [3..(2+L)]    = sizes
  #     [3+L.. ]      = x[1..n1]
  # ---------------------------------------------------------------------------
  inputs <- c(
    as.double(n1),
    as.double(L),
    as.double(sizes),
    as.double(x)
  )

  # ---------------------------------------------------------------------------
  # Coefficient name
  # ---------------------------------------------------------------------------
  pretty_sizes <- if (L == 0L) "all" else paste0("S{", paste(sizes, collapse = ","), "}")
  coef.names   <- paste0("cov_ingroup[", cov_label, "]_", pretty_sizes)

  dbgcat("coef.names =", coef.names, " | inputs length =", length(inputs))

  # ---------------------------------------------------------------------------
  # Standard ERGM term specification
  # ---------------------------------------------------------------------------
  # IMPORTANT:
  # - d_func = TRUE tells ergm to call the multi-toggle (D_) changestat entrypoint
  #   implemented as D_CHANGESTAT_FN(d_cov_ingroup).
  list(
    name         = "cov_ingroup",
    coef.names   = coef.names,
    inputs       = inputs,
    dependence   = TRUE,
    d_func       = TRUE,  # <-- REQUIRED (multi-toggle)
    emptynwstats = 0
  )
}