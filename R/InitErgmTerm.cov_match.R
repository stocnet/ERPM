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
# where:
#   - n_{g,r} is the number of actors of category r in group g,
#   - n_g     is the size of group g.
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
#' actor-level covariate. The network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side of the bipartite graph).
#' }
#'
#' For each group in the group mode and each category of the actor covariate,
#' the term considers all subsets of actors in that group who share the same
#' category and counts cliques of size \eqn{k}. Several variants of
#' normalization are supported:
#' \itemize{
#'   \item \code{normalized = "none"}: raw counts of monochromatic cliques;
#'   \item \code{normalized = "by_group"}: per-group normalization by
#'         \eqn{\binom{n_g}{k}}, where \eqn{n_g} is group size;
#'   \item \code{normalized = "global"}: per-group normalization by
#'         \eqn{n_g}, i.e. each group contribution is divided by its size
#'         and the resulting terms are summed over groups.
#' }
#' An optional targeted category \eqn{\kappa} focuses the statistic on cliques
#' whose actors all share that specific category.
#'
#' @details
#' Let:
#' \itemize{
#'   \item \eqn{A} be the set of actor-mode nodes, with \eqn{|A| = N_A};
#'   \item \eqn{G} be the set of group-mode nodes;
#'   \item \eqn{B} be the actor–group incidence (bipartite) matrix;
#'   \item \eqn{c : A \to \{1,\dots,R\}} be a categorical covariate assigning a
#'         category \eqn{r} to each actor;
#'   \item \eqn{n_{g,r}} be the number of actors of category \eqn{r} attached to
#'         group \eqn{g};
#'   \item \eqn{n_g = \sum_r n_{g,r}} be the size of group \eqn{g}.
#' }
#'
#' For a fixed clique size \eqn{k \ge 1}, define:
#'
#' \code{
#'   S_k(B; c)       = sum_g sum_r C(n_{g,r}, k),
#' }
#'
#' i.e. the total number of size-\eqn{k} monochromatic subsets of actors within
#' all groups.
#'
#' If a targeted category \eqn{\kappa} is specified, define:
#'
#' \code{
#'   S_k^{(κ)}(B; c) = sum_g C(n_{g,κ}, k),
#' }
#'
#' i.e. the count restricted to actors whose category equals \eqn{\kappa}.
#'
#' The three normalization modes correspond to:
#' \itemize{
#'   \item \code{"none"}:
#'     \deqn{
#'       T_k(B; c) =
#'       \begin{cases}
#'         S_k(B; c)             & \text{if no category is targeted}, \\
#'         S_k^{(\kappa)}(B; c)  & \text{if category } \kappa \text{ is targeted};
#'       \end{cases}
#'     }
#'   \item \code{"by_group"}: for each group \eqn{g} we form the ratio
#'         \eqn{\frac{S_k(g)}{\binom{n_g}{k}}} (or
#'         \eqn{\frac{C(n_{g,\kappa}, k)}{\binom{n_g}{k}}} when a category is
#'         targeted), and sum these ratios over groups;
#'   \item \code{"global"}: for each group \eqn{g} we form
#'     \deqn{
#'       \frac{S_k(g)}{n_g}
#'       \quad\text{or}\quad
#'       \frac{C(n_{g,\kappa}, k)}{n_g}
#'     }
#'     (with the convention that groups of size \eqn{n_g = 0} contribute 0),
#'     and sum these contributions over groups:
#'     \deqn{
#'       T_k^{\text{global}}(B; c) =
#'       \sum_{g \in G}
#'       \begin{cases}
#'         \dfrac{S_k(g)}{n_g}            & \text{if no category is targeted}, \\
#'         \dfrac{C(n_{g,\kappa}, k)}{n_g}& \text{if category } \kappa \text{ is targeted}.
#'       \end{cases}
#'     }
#' }
#'
#' The term is vectorized in \code{clique_size}: several values of \eqn{k} can
#' be specified, yielding one statistic per \eqn{k}.
#'
#' The term is implemented as a native ERGM C change-statistic
#' \code{c_cov_match}. The R initializer below:
#' \itemize{
#'   \item enforces that the network is bipartite and retrieves the actor-mode
#'         size from \code{nw \%n\% "bipartite"};
#'   \item validates that the actor covariate is categorical (factor/character)
#'         and encodes it into integer codes \eqn{0, 1, \dots, R};
#'   \item handles the optional targeted category by mapping it to its level
#'         index \eqn{\kappa};
#'   \item maps the normalization choice to an internal integer flag;
#'   \item builds a compact \code{INPUT_PARAM} vector encoding the actor-mode
#'         size, the clique sizes, normalization mode, and covariate codes.
#' }
#'
#' @section INPUT_PARAM layout (C side):
#' The numeric input vector passed to \code{c_cov_match} is:
#'
#' \code{
#'   INPUT_PARAM = c(
#'     n1,          # actor-mode size |A|
#'     K,           # number of distinct clique sizes
#'     norm_mode,   # 0=none, 1=by_group, 2=global
#'     has_kappa,   # 0/1: whether a targeted category is used
#'     kappa_code,  # level index of the targeted category (0 if none)
#'     ks[1:K],     # vector of clique sizes k >= 1
#'     z[1:n1]      # actor covariate codes (0=missing / undefined, 1..R for levels)
#'   )
#' }
#'
#' On each toggle of an actor–group edge, the C code recomputes the local
#' contribution for the affected group and updates the statistic accordingly,
#' respecting the chosen normalization and targeted category.
#'
#' @section Arguments:
#' The initializer is called internally by \pkg{ergm} and should not be invoked
#' directly by users. The user-facing term is:
#'
#' \code{
#'   cov_match(cov,
#'             clique_size = 2,
#'             category    = NULL,
#'             normalized  = c("none","by_group","global"))
#' }
#'
#' The term arguments are passed to this initializer through \code{arglist}
#' with expected components:
#' \itemize{
#'   \item \code{cov}: character (vertex attribute name) or a factor/character vector;
#'   \item \code{clique_size}: integer(s) \eqn{k \ge 1};
#'   \item \code{category}: optional targeted category (character or \code{NULL});
#'   \item \code{normalized}: \code{"none"}, \code{"by_group"}, or \code{"global"}
#'         (logical values are accepted as shorthand).
#' }
#'
#' @param nw A \pkg{network} object.
#' @param arglist A named list of term arguments (see \sQuote{Arguments}).
#' @param ... Passed through by \pkg{ergm}; not used.
#' @param version ERGM API version; not used.
#'
#' @return
#' A standard \pkg{ergm} term specification list with components:
#' \itemize{
#'   \item \code{name}         = \code{"cov_match"};
#'   \item \code{coef.names}   = coefficient names encoding the covariate label,
#'         clique size, and normalization mode;
#'   \item \code{inputs}       = the \code{INPUT_PARAM} numeric vector described
#'         above;
#'   \item \code{dependence}   = \code{TRUE};
#'   \item \code{emptynwstats} = \code{0}.
#' }
#'
#' @note
#' \itemize{
#'   \item The network must be bipartite and interpreted as actors versus groups.
#'         The actor mode is identified by \code{nw \%n\% "bipartite"} and must
#'         be a strictly positive integer.
#'   \item The covariate must be categorical (factor or character). Numeric
#'         vectors are rejected fail-fast because \code{cov_match} relies on
#'         category frequencies, not continuous values.
#'   \item By default, \code{clique_size = 1} with \code{normalized = "none"} or
#'         \code{"global"} is disallowed, because the statistic is then constant
#'         under edge toggles in the typical ERPM partition setting. This behavior
#'         can be overridden (for advanced use) by setting:
#'         \code{options(ERPM.allow.k1.nonnormalized = TRUE)}.
#'   \item Debug logging for the initializer can be enabled via:
#'         \code{options(erpm.debug.cov_match_init = TRUE)}. When enabled, the
#'         initializer prints diagnostic information about actor-mode size,
#'         clique sizes, normalization mode, and covariate level mapping.
#'         It also verifies whether ERPM wrapper metadata are attached:
#'         \code{nw \%n\% "nodes"} and \code{nw \%n\% "dyads"}.
#' }
#'
#' @examples
#' \dontrun{
#'   library(network)
#'   library(ergm)
#'
#'   # -----------------------------------------------------------------------
#'   # Build a small bipartite network: 5 actors, 2 groups
#'   # -----------------------------------------------------------------------
#'   n_actors <- 5
#'   n_groups <- 2
#'   n_total  <- n_actors + n_groups
#'
#'   adj <- matrix(0, n_total, n_total)
#'
#'   # Actors = 1..5, Groups = 6..7
#'   # Group 6: actors 1, 2, 3
#'   adj[1, 6] <- adj[6, 1] <- 1
#'   adj[2, 6] <- adj[6, 2] <- 1
#'   adj[3, 6] <- adj[6, 3] <- 1
#'   # Group 7: actors 3, 4, 5
#'   adj[3, 7] <- adj[7, 3] <- 1
#'   adj[4, 7] <- adj[7, 4] <- 1
#'   adj[5, 7] <- adj[7, 5] <- 1
#'
#'   nw <- network(adj, directed = FALSE, matrix.type = "adjacency")
#'   nw %n% "bipartite" <- n_actors  # actor-mode size
#'
#'   # Actor covariate: two categories "A" / "B"
#'   cov_vals <- c("A", "A", "B", "B", "A")
#'   set.vertex.attribute(nw, "grp", c(cov_vals, rep(NA, n_groups)))
#'
#'   # -----------------------------------------------------------------------
#'   # Example 1: raw counts of k=2 monochromatic cliques (pairs)
#'   # -----------------------------------------------------------------------
#'   summary(
#'     nw ~ cov_match("grp", clique_size = 2, normalized = "none"),
#'     constraints = ~ b1part
#'   )
#'
#'   # -----------------------------------------------------------------------
#'   # Example 2: k=2 with by-group normalization
#'   # -----------------------------------------------------------------------
#'   summary(
#'     nw ~ cov_match("grp", clique_size = 2, normalized = "by_group"),
#'     constraints = ~ b1part
#'   )
#'
#'   # -----------------------------------------------------------------------
#'   # Example 3: targeted category, k=3, global normalization (by group size)
#'   # -----------------------------------------------------------------------
#'   summary(
#'     nw ~ cov_match("grp", clique_size = 3,
#'                    category   = "A",
#'                    normalized = "global"),
#'     constraints = ~ b1part
#'   )
#'
#'   # Example ERGM fit with two clique sizes
#'   fit <- ergm(
#'     nw ~ cov_match("grp", clique_size = c(2, 3), normalized = "by_group"),
#'     constraints = ~ b1part
#'   )
#'   summary(fit)
#'
#'   # Example call through the ERPM wrapper (network already bipartite)
#'   # erpm(nw ~ cov_match("grp", clique_size = 2, normalized = "by_group"))
#' }
#'
#' @section Tests:
#' Self-tests for \code{cov_match} (not shown here) typically:
#' \itemize{
#'   \item construct small bipartite networks with an explicit partition of
#'         actors into groups and a known categorical covariate;
#'   \item compute the reference counts of monochromatic cliques of size
#'         \eqn{k} directly in R, both with and without a targeted category,
#'         and under each normalization mode ("none", "by_group", "global");
#'   \item compare these reference values to
#'         \code{summary(nw ~ cov_match(...), constraints = ~ b1part)};
#'   \item check that toggling a single actor–group edge modifies the statistic
#'         by a difference that matches the local recomputation of clique counts
#'         in the affected group, in agreement with the C change-statistic
#'         \code{c_cov_match}.
#' }
#'
#' @keywords ERGM term bipartite categorical covariate monochromatic cliques
#' @md
#'
#' @export
InitErgmTerm.cov_match <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_match"

  # -------------------------------------------------------------------------
  # Optional debug flag for this initializer
  #   options(erpm.debug.cov_match_init = TRUE)
  # will trigger additional console output.
  # -------------------------------------------------------------------------
  DEBUG <- isTRUE(getOption("erpm.debug.cov_match_init", FALSE))

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
  if (!is.finite(n1) || n1 <= 0L)
    ergm_Init_stop(sQuote(termname), ": non-bipartite network or missing/invalid %n% 'bipartite' attribute.")

  # -------------------------------------------------------------------------
  # Debug: verify ERPM wrapper metadata are attached
  #   - build_bipartite_from_inputs() stores dyads in: nw %n% "dyads"
  #   - nodes may be stored in: nw %n% "nodes"
  # This does not change semantics. It only prints diagnostics when enabled.
  # -------------------------------------------------------------------------
  if (DEBUG) {
    .safe_get_n_attr <- function(nw, key) {
      tryCatch(nw %n% key, error = function(e) NULL)
    }

    nodes_meta <- .safe_get_n_attr(nw, "nodes")
    dyads_meta <- .safe_get_n_attr(nw, "dyads")

    cat(sprintf("[Init:%s] debug=TRUE\n", termname))
    cat(sprintf("[Init:%s] network size=%d | bipartite n1=%d\n", termname, network::network.size(nw), n1))

    if (is.null(nodes_meta)) {
      cat(sprintf("[Init:%s] nw %%n%% \"nodes\" : ABSENT\n", termname))
    } else {
      cls <- paste(class(nodes_meta), collapse = "/")
      cat(sprintf("[Init:%s] nw %%n%% \"nodes\" : PRESENT (class=%s)\n", termname, cls))
      if (is.data.frame(nodes_meta)) {
        cat(sprintf("[Init:%s] nodes: nrow=%d ncol=%d\n", termname, nrow(nodes_meta), ncol(nodes_meta)))
        cat(sprintf("[Init:%s] nodes colnames: %s\n", termname, paste(colnames(nodes_meta), collapse = ", ")))
        cat(sprintf("[Init:%s] nodes head:\n", termname))
        print(utils::head(nodes_meta, 3))
      }
    }

    if (is.null(dyads_meta)) {
      cat(sprintf("[Init:%s] nw %%n%% \"dyads\" : ABSENT\n", termname))
    } else {
      cls <- paste(class(dyads_meta), collapse = "/")
      cat(sprintf("[Init:%s] nw %%n%% \"dyads\" : PRESENT (class=%s)\n", termname, cls))
      if (is.data.frame(dyads_meta)) {
        cat(sprintf("[Init:%s] dyads: nrow=%d ncol=%d\n", termname, nrow(dyads_meta), ncol(dyads_meta)))
        cat(sprintf("[Init:%s] dyads colnames: %s\n", termname, paste(colnames(dyads_meta), collapse = ", ")))
        cat(sprintf("[Init:%s] dyads head:\n", termname))
        print(utils::head(dyads_meta, 3))
      }
    }
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
  if (length(ks) < 1L || any(!is.finite(ks)) || any(ks < 1L))
    ergm_Init_stop(sQuote(termname), ": 'clique_size' must contain finite integers >= 1.")

  # Option to explicitly allow k = 1 with non-normalized/global modes
  .allow_k1_nn <- isTRUE(getOption("ERPM.allow.k1.nonnormalized", FALSE))

  # Forbid k = 1 when normalized is "none" or "global" unless explicitly overridden
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
      # Case 1: 'cov' is the name of a vertex attribute
      vals <- network::get.vertex.attribute(nw, cov)
      if (is.null(vals))
        ergm_Init_stop(sQuote(termname), ": missing vertex attribute: ", sQuote(cov), ".")

      # Only actor-mode values are used.
      x <- vals[ia]

      if (is.numeric(x))
        ergm_Init_stop(sQuote(termname), ": 'cov_match' requires a categorical covariate (factor/character), not numeric.")

      f <- as.factor(x)

      # If a targeted category is absent, add it to levels so match() is defined.
      if (!is.null(category)) {
        category <- as.character(category)[1L]
        if (!(category %in% levels(f))) levels(f) <- c(levels(f), category)
      }

      z <- as.integer(f)            # 1..R or NA
      z[is.na(z)] <- 0L            # 0 = "absent/undefined"
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

    # Case 2: direct vector supplied
    if (is.numeric(cov)) {
      ergm_Init_stop(sQuote(termname), ": 'cov_match' requires a categorical covariate (factor/character), not a numeric vector.")
    }

    if (length(cov) < n1)
      ergm_Init_stop(sQuote(termname), ": length(cov) < |A| = ", n1, ".")

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
  levs       <- ax$levels

  has_kappa <- as.double(as.integer(kappa_code > 0))

  # -------------------------------------------------------------------------
  # Optional local debug output (continued)
  # -------------------------------------------------------------------------
  if (DEBUG) {
    cat(sprintf("[Init:%s] normalized=%s (mode=%d) | K=%d | ks={%s}\n",
                termname, normalized, norm_mode, K, paste(ks, collapse=",")))
  }

  # -------------------------------------------------------------------------
  # 4) Build INPUT_PARAM vector for the C change-statistic
  # -------------------------------------------------------------------------
  # Layout:
  #   [0]          n1
  #   [1]          K
  #   [2]          norm_mode  (0 none, 1 by_group, 2 global)
  #   [3]          has_kappa  (0/1)
  #   [4]          kappa_code (0 if no targeted category)
  #   [5 .. 5+K-1] ks (clique sizes)
  #   [5+K .. ]    z[1..n1] (actor covariate codes)
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
  # Example patterns:
  #   cov_match[sex]_k2
  #   cov_match[sex==F]_k3_bygrp
  #   cov_match[group]_k4_glob
  suffix_norm <- switch(normalized,
                        none     = "",
                        by_group = "_bygrp",
                        global   = "_glob")
  coef.names  <- paste0("cov_match[", cov_label, "]_k", ks, suffix_norm)

  # -------------------------------------------------------------------------
  # 6) Standard ERGM term specification
  # -------------------------------------------------------------------------
  list(
    name         = "cov_match",         # must match C_CHANGESTAT_FN(c_cov_match)
    coef.names   = coef.names,          # length = K
    inputs       = inputs,              # as described above
    dependence   = TRUE,
    emptynwstats = rep(0, K)
  )
}
