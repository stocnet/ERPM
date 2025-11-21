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
#     (1 / C(N_A, k)) * S_k(·)
# where:
#   - N_A is the number of actors in the actor mode,
#   - n_{g,r} is the number of actors of category r in group g.
# ==============================================================================

#' ERGM term: cov_match (monochromatic cliques by actor covariate)
#'
#' @file InitErgmTerm.cov_match.R
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
#'   \item \code{normalized = "global"}: global normalization by
#'         \eqn{\binom{N_A}{k}}, where \eqn{N_A} is the size of the actor mode.
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
#'         S_k(B; c)       & \text{if no category is targeted}, \\
#'         S_k^{(\kappa)}(B; c) & \text{if category } \kappa \text{ is targeted};
#'       \end{cases}
#'     }
#'   \item \code{"by_group"}: for each group \eqn{g} we form the ratio
#'         \eqn{\frac{S_k(g)}{\binom{n_g}{k}}} (or
#'         \eqn{\frac{C(n_{g,\kappa}, k)}{\binom{n_g}{k}}} when a category is
#'         targeted), and sum these ratios over groups;
#'   \item \code{"global"}: we normalize by \eqn{\binom{N_A}{k}}:
#'     \deqn{
#'       T_k^{\text{global}}(B; c) =
#'       \frac{1}{\binom{N_A}{k}}
#'       \times
#'       \begin{cases}
#'         S_k(B; c)       & \text{if no category is targeted}, \\
#'         S_k^{(\kappa)}(B; c) & \text{if category } \kappa \text{ is targeted}.
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
#' @code{
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
#' The initializer is called internally by {ergm} and should not be invoked
#' directly by users. The user-facing term is:
#'
#' @code{
#'   cov_match(cov,
#'             clique_size = 2,
#'             category    = NULL,
#'             normalized  = c("none","by_group","global"))
#' }
#'
#' @param cov character|factor|numeric
#'   Either:
#'   \itemize{
#'     \item the name of an actor-level vertex attribute (factor or character)
#'           defined on all actors in the actor mode; or
#'     \item a vector of length at least \code{|A|} (number of actors), which
#'           is interpreted as an actor-level covariate.
#'   }
#'   In both cases, the covariate is coerced to a factor and then encoded as
#'   integer codes \eqn{1,\dots,R}; \code{NA} values are mapped to 0 (meaning
#'   "absent/undefined" and ignored in clique counts). For numeric covariates,
#'   \code{cov_match} is not meaningful and the initializer fails fast with an
#'   error.
#'
#' @param clique_size integer|numeric
#'   One or several clique sizes \eqn{k \ge 1}. The term is vectorized in
#'   \code{clique_size}, so each distinct value of \eqn{k} yields one statistic
#'   and one coefficient. Values are rounded to integers and must be finite and
#'   at least 1.
#'
#' @param category character|NULL
#'   Optional targeted category. If \code{NULL}, all categories contribute to
#'   the statistic. If a character string, the initializer ensures that the
#'   category appears in the factor levels; if it does not, the level is added
#'   with zero frequency so that the resulting statistic is structurally zero
#'   without error.
#'
#' @param normalized character|logical
#'   Normalization mode, one of:
#'   \itemize{
#'     \item \code{"none"}: raw counts of monochromatic cliques;
#'     \item \code{"by_group"}: per-group normalization by \eqn{\binom{n_g}{k}};
#'     \item \code{"global"}: normalization by \eqn{\binom{N_A}{k}}.
#'   }
#'   Logical values are supported as shorthand:
#'   \code{TRUE} is equivalent to \code{"by_group"} and \code{FALSE} to
#'   \code{"none"}.
#'
#' @return
#' A standard {ergm} term specification list with components:
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
#'         under edge toggles and uninformative for ERGM fitting. This behavior
#'         can be overridden (for advanced use) by setting:
#'         \code{options(ERPM.allow.k1.nonnormalized = TRUE)}.
#'   \item Debug logging for the initializer can be enabled via:
#'         \code{options(erpm.debug.cov_match_init = TRUE)}. When enabled, the
#'         initializer prints diagnostic information about actor-mode size,
#'         clique sizes, normalization mode, and covariate level mapping.
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
#'   # Example 3: targeted category, k=3, global normalization
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
#' @test
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
    options(erpm.debug.cov_match_init = FALSE)

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
    n1 <- tryCatch(nw %n% "bipartite", error = function(e) NA_integer_)
    if (!is.numeric(n1) || !is.finite(n1) || n1 <= 0)
        ergm_Init_stop(sQuote(termname), ": réseau non biparti ou attribut %n% 'bipartite' manquant/invalide.")

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
    ks <- as.integer(a$clique_size)
    if (any(!is.finite(ks)) || any(ks < 1L))
        ergm_Init_abort(sQuote(termname), ": 'clique_size' doit contenir des entiers >= 1.")

    # Option to explicitly allow k = 1 with non-normalized/global modes
    .allow_k1_nn <- isTRUE(getOption("ERPM.allow.k1.nonnormalized", FALSE))

    # Forbid k = 1 when normalized is "none" or "global" unless explicitly overridden
    if (any(ks == 1L) && (normalized %in% c("none","global")) && !.allow_k1_nn) {
        ergm_Init_abort(
            sQuote(termname),
            ": cov_match(..., clique_size=1) avec normalized='", normalized,
            "' est constant pour les fits ERGM. ",
            "Utilisez k>=2, normalized='by_group', ou offset(...). ",
            "Pour forcer: options(ERPM.allow.k1.nonnormalized=TRUE)."
        )
    }

    ks <- unique(ks)
    K  <- length(ks)

    # -------------------------------------------------------------------------
    # 3) Build actor-level category codes (z) and targeted category info
    # -------------------------------------------------------------------------
    cov      <- a$cov
    category <- a$category

    get_actor_codes <- function(nw, cov, category = NULL) {
        # Determine actor indices: by default 1..n1, but try to detect actor
        # labels when vertex names suggest a pattern "G<id>" for groups.
        n1 <- as.integer(nw %n% "bipartite")
        ia <- seq_len(n1)
        vn <- network::network.vertex.names(nw)
        if (length(vn) >= n1) {
            ia_guess <- which(!grepl("^G\\d+$", vn))
            if (length(ia_guess) == n1) ia <- ia_guess
        }

        if (is.character(cov) && length(cov) == 1L) {
            # Case 1: 'cov' is the name of a vertex attribute
            vals <- network::get.vertex.attribute(nw, cov)
            if (is.null(vals))
                ergm_Init_stop(sQuote(termname), ": attribut inexistant: ", sQuote(cov), ".")
            x <- vals[ia]

            # Coerce to factor for clean category coding
            f <- as.factor(x)
            z <- as.integer(f)     # codes 1..R, NA -> NA
            z[!is.finite(z)] <- 0L # 0 = "absent/undefined"

            kappa_code <- 0L
            cov_label  <- cov
            if (!is.null(category)) {
                # If the targeted category is absent, add it to the levels
                # so that it has zero frequency without error.
                if (!(category %in% levels(f))) {
                    levels(f) <- c(levels(f), category)
                    # Recompute codes with updated levels
                    z <- as.integer(f)
                    z[!is.finite(z)] <- 0L
                }
                kappa_code <- as.integer(match(category, levels(f)))
                cov_label  <- paste0(cov, "==", category)
            }

            # Return codes, targeted category code, and levels (for debug)
            return(list(
                z          = as.double(z),
                kappa_code = as.double(kappa_code),
                cov_label  = cov_label,
                levels     = levels(f)
            ))

        } else {
            # Case 2: direct numeric vector
            # For cov_match, a purely numeric covariate is not meaningful,
            # so we fail fast after basic checks.
            x_num <- suppressWarnings(as.numeric(cov))
            if (any(!is.finite(x_num)))
                ergm_Init_stop(sQuote(termname), ": vecteur 'cov' contient NA/NaN/Inf.")
            if (length(x_num) < n1)
                ergm_Init_stop(sQuote(termname), ": longueur(cov) < |A| = ", n1, ".")
            ergm_Init_stop(sQuote(termname), ": 'cov_match' requiert un attribut catégoriel (factor/character).")
        }
    }

    ax <- get_actor_codes(nw, cov = cov, category = category)
    z_codes    <- ax$z          # double*, codes 0 (NA) or 1..R
    kappa_code <- ax$kappa_code # 0 if no targeted category
    cov_label  <- ax$cov_label
    levs       <- ax$levels %||% character(0)

    has_kappa <- as.double(as.integer(kappa_code > 0))

    # -------------------------------------------------------------------------
    # Optional local debug output
    #   Enable with: options(erpm.debug.cov_match_init = TRUE)
    # -------------------------------------------------------------------------
    DEBUG_COV_MATCH_INIT <- isTRUE(getOption("erpm.debug.cov_match_init", FALSE)) || FALSE
    if (DEBUG_COV_MATCH_INIT) {
        cat(sprintf("[Init:%s] n1=%d | normalized=%s (mode=%d) | K=%d | ks={%s}\n",
                    termname, n1, normalized, norm_mode, K, paste(ks, collapse=",")))
        cat(sprintf("[Init:%s] cov=%s | has_kappa=%s | kappa_code=%s\n",
                    termname, deparse(substitute(cov)), as.logical(has_kappa), as.integer(kappa_code)))
        if (length(levs)) {
            map <- paste(sprintf("%d:%s", seq_along(levs), levs), collapse=", ")
            cat(sprintf("[Init:%s] levels map: %s\n", termname, map))
            zi <- as.integer(z_codes)
            zi[!is.finite(zi)] <- 0L
            tab <- as.integer(table(factor(zi, levels = 0:length(levs))))
            cat(sprintf("[Init:%s] freq codes (0=NA): {%s}\n",
                        termname, paste(tab, collapse=", ")))
        }
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
        name         = "cov_match",   # must match C_CHANGESTAT_FN(c_cov_match)
        coef.names   = coef.names,    # length = K
        inputs       = inputs,        # as described above
        dependence   = TRUE,
        emptynwstats = 0
    )
}