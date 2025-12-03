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
#   - B be the actor–group incidence matrix;
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
#' \code{c_cov_ingroup}. The R initializer:
#' \itemize{
#'   \item enforces that the network is bipartite and retrieves the actor-mode
#'         size from \code{nw \%n\% "bipartite"};
#'   \item builds a numeric vector \eqn{x \in \mathbb{R}^{n_1}} of covariate
#'         values on actors, either from a vertex attribute or a direct vector;
#'   \item normalizes and validates the size filter \code{size}, converting it
#'         into a sorted set \eqn{S} of positive integers;
#'   \item packs \eqn{n_1}, the size filter and \eqn{x} into a compact
#'         \code{INPUT_PARAM} vector consumed by the C code.
#' }
#'
#' On each toggle of an actor–group edge, the C change-statistic recomputes the
#' contribution of the affected group and updates the statistic in \eqn{O(n_g)}
#' time for that group.
#'
#' @section Mathematical definition:
#' Let:
#' \itemize{
#'   \item \eqn{A} denote the set of actor-mode nodes with \eqn{|A| = n_1};
#'   \item \eqn{G} denote the set of group-mode nodes;
#'   \item \eqn{B} be the bipartite adjacency between actors and groups;
#'   \item \eqn{x : A \to \mathbb{R}} be a numeric covariate on actors;
#'   \item \eqn{S \subseteq \mathbb{N}} be a set of admissible group sizes.
#' }
#' For each group \eqn{g \in G}, define:
#' \itemize{
#'   \item \eqn{n_g = \sum_{i \in A} B_{i,g}} the group size (number of actors);
#'   \item \eqn{X_g = \sum_{i \in A} B_{i,g} x_i} the sum of covariate values
#'         of actors in group \eqn{g}.
#' }
#' The statistic is:
#' \deqn{
#'   T(B; x, S) = \sum_{g \in G} n_g X_g \mathbf{1}[n_g \in S].
#' }
#' When a categorical attribute \eqn{c_i} and a category \eqn{\kappa} are used,
#' we set \eqn{x_i = \mathbf{1}[c_i = \kappa]} and obtain:
#' \deqn{
#'   T(B; \kappa, S) =
#'   \sum_{g \in G} n_g \left(\sum_{i \in g} \mathbf{1}[c_i = \kappa]\right)
#'   \mathbf{1}[n_g \in S].
#' }
#'
#' @section INPUT_PARAM layout (C side):
#' The numeric vector passed to \code{c_cov_ingroup} has the following layout:
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
#' The C code reconstructs:
#'
#' \itemize{
#'   \item the actor-mode size \eqn{n_1};
#'   \item the set \eqn{S} of admissible group sizes;
#'   \item the numeric covariate vector \eqn{x};
#' }
#' and then recomputes the local changes to \eqn{T(B; x, S)} when edges between
#' actors and groups are toggled.
#'
#' @section Usage:
#' The user-facing term is:
#'
#' \preformatted{
#'   cov_ingroup(cov,
#'               size     = NULL,
#'               category = NULL)
#' }
#'
#' Typical usage in an ERGM formula:
#'
#' \preformatted{
#'   # Bipartite network with actor mode A and group mode G
#'   summary(nw ~ cov_ingroup("age") + b1part)
#'
#'   # Restrict to groups with sizes 3, 4 or 5
#'   summary(nw ~ cov_ingroup("age", size = 3:5) + b1part)
#'
#'   # Target a specific category of a categorical covariate
#'   summary(nw ~ cov_ingroup("gender", category = "F") + b1part)
#' }
#'
#' When using the ERPM wrapper, the same term can be used either on a bipartite
#' network or directly on a partition representation:
#'
#' \preformatted{
#'   erpm(nw ~ cov_ingroup("age"))
#'   erpm(partition ~ cov_ingroup("gender", category = "F"))
#' }
#'
#' @param cov character|numeric  
#'   Either:
#'   \itemize{
#'     \item the name of an actor-level vertex attribute (numeric or
#'           categorical) to be evaluated on the actor mode; or
#'     \item a numeric vector of length at least \eqn{|A| = n_1} giving covariate
#'           values directly for actors.
#'   }
#'   When \code{category} is \code{NULL} and \code{cov} refers to a vertex
#'   attribute, the attribute is coerced to numeric and used as is. When
#'   \code{category} is not \code{NULL} and \code{cov} refers to a vertex
#'   attribute, an indicator is constructed for the targeted category
#'   (\code{1[cov == category]} on the actor mode). If \code{cov} is a numeric
#'   vector, \code{category} must be \code{NULL}.
#'
#' @param size integer|numeric|NULL  
#'   Optional set \eqn{S} of admissible group sizes. If \code{NULL} or empty, all
#'   group sizes are included. Otherwise, \code{size} is converted to a sorted
#'   vector of distinct positive integers; only groups with size in this set
#'   contribute to the statistic.
#'
#' @param category character|NULL  
#'   Optional targeted category when \code{cov} denotes a categorical actor
#'   attribute. If provided, the actor covariate is replaced by the indicator
#'   \eqn{1[c_i = \text{category}]}. If \code{cov} is given as a numeric vector,
#'   \code{category} must remain \code{NULL}.
#'
#' @return
#' A standard {ergm} term specification list with components:
#' \itemize{
#'   \item \code{name}         = \code{"cov_ingroup"};
#'   \item \code{coef.names}   = a single coefficient name encoding the
#'         covariate label and the size filter;
#'   \item \code{inputs}       = numeric vector \code{INPUT_PARAM} as described
#'         above;
#'   \item \code{dependence}   = \code{TRUE};
#'   \item \code{emptynwstats} = \code{0}.
#' }
#'
#' @note
#' \itemize{
#'   \item The network must be bipartite and interpreted as actors versus
#'         groups. The actor mode size is taken from \code{nw \%n\% "bipartite"}
#'         and must be a strictly positive finite integer.
#'   \item When \code{cov} is given as an attribute name and \code{category} is
#'         \code{NULL}, the attribute is coerced to numeric; non-finite values
#'         (\code{NA}, \code{NaN}, \code{Inf}) are rejected in a fail-fast
#'         manner with a descriptive error.
#'   \item When \code{category} is provided, the underlying attribute is
#'         compared to the category at the R level, and missing values are
#'         treated as non-matching (indicator 0).
#' }
#'
#' @examples
#' \dontrun{
#'   library(network)
#'   library(ergm)
#'
#'   # -----------------------------------------------------------------------
#'   # Build a small bipartite network: 4 actors, 2 groups
#'   # -----------------------------------------------------------------------
#'   n_actors <- 4
#'   n_groups <- 2
#'   n_total  <- n_actors + n_groups
#'
#'   adj <- matrix(0, n_total, n_total)
#'
#'   # Actors = 1..4, Groups = 5..6
#'   # Group 5: actors 1, 2
#'   adj[1, 5] <- adj[5, 1] <- 1
#'   adj[2, 5] <- adj[5, 2] <- 1
#'   # Group 6: actors 2, 3, 4
#'   adj[2, 6] <- adj[6, 2] <- 1
#'   adj[3, 6] <- adj[6, 3] <- 1
#'   adj[4, 6] <- adj[6, 4] <- 1
#'
#'   nw <- network(adj, directed = FALSE, matrix.type = "adjacency")
#'   nw \%n\% "bipartite" <- n_actors  # actor-mode size
#'
#'   # Numeric actor covariate (e.g. age)
#'   age <- c(25, 30, 28, 40)
#'   set.vertex.attribute(nw, "age", c(age, rep(NA_real_, n_groups)))
#'
#'   # Categorical actor covariate (e.g. gender)
#'   gender <- c("F", "M", "F", "M")
#'   set.vertex.attribute(nw, "gender", c(gender, rep(NA_character_, n_groups)))
#'
#'   # -----------------------------------------------------------------------
#'   # Example 1: numeric covariate, all group sizes
#'   # -----------------------------------------------------------------------
#'   summary(
#'     nw ~ cov_ingroup("age") + b1part
#'   )
#'
#'   # -----------------------------------------------------------------------
#'   # Example 2: numeric covariate, filter groups with size in {2, 3}
#'   # -----------------------------------------------------------------------
#'   summary(
#'     nw ~ cov_ingroup("age", size = c(2, 3)) + b1part
#'   )
#'
#'   # -----------------------------------------------------------------------
#'   # Example 3: categorical covariate, targeted category "F"
#'   # -----------------------------------------------------------------------
#'   summary(
#'     nw ~ cov_ingroup("gender", category = "F") + b1part
#'   )
#'
#'   # Example ERGM fit
#'   fit <- ergm(
#'     nw ~ cov_ingroup("age", size = c(2, 3)) + b1part
#'   )
#'   summary(fit)
#' }
#'
#' @section Tests:
#' Self-tests for \code{cov_ingroup} (not shown here) typically:
#' \itemize{
#'   \item build small bipartite networks with a known partition of actors into
#'         groups and known covariate values on the actor mode;
#'   \item compute, in pure R, the reference value
#'         \eqn{T(B; x, S) = \sum_g n_g (\sum_{i \in g} x_i)\mathbf{1}[n_g \in S]};
#'   \item compare these reference values to
#'         \code{summary(nw ~ cov_ingroup(...), constraints = ~ b1part)};
#'   \item verify that toggling a single actor–group tie changes the statistic
#'         by the local increment obtained by recomputing the contribution of
#'         the affected group only, as implemented in the C change-statistic
#'         \code{c_cov_ingroup}.
#' }
#'
#' @keywords ERGM term bipartite covariate groups
#' @md
#'
#' @export
InitErgmTerm.cov_ingroup <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_ingroup"

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
  if (!is.numeric(n1) || !is.finite(n1) || n1 <= 0)
    ergm_Init_stop(sQuote(termname), ": réseau non biparti ou attribut %n% 'bipartite' manquant/invalide.")

  # ---------------------------------------------------------------------------
  # Build the actor-level covariate vector x (length n1)
  #   - handle both attribute name and direct numeric vector
  #   - handle optional category for categorical attributes
  # ---------------------------------------------------------------------------
  cov      <- a$cov
  category <- a$category

  get_actor_cov <- function(nw, cov, category = NULL) {
    n1 <- as.integer(nw %n% "bipartite")

    # Indices for the actor mode:
    # by convention from the builder: actors are listed first (1..n1).
    # If vertex names suggest an actor/group split, use that as a best effort.
    ia <- seq_len(n1)
    vn <- network::network.vertex.names(nw)
    if (length(vn) >= n1) {
      ia_guess <- which(!grepl("^G\\d+$", vn))
      if (length(ia_guess) == n1) ia <- ia_guess
    }

    # Case 1: cov is the name of a vertex attribute
    if (is.character(cov) && length(cov) == 1L) {
      vals <- network::get.vertex.attribute(nw, cov)
      if (is.null(vals))
        ergm_Init_stop(sQuote(termname), ": attribut inexistant: ", sQuote(cov), ".")
      x <- vals[ia]

      # If a category is provided, build an indicator covariate x_i = 1[c_i == category]
      if (!is.null(category)) {
        xb <- as.integer(as.character(x) == category)
        xb[is.na(xb)] <- 0L
        return(list(
          x         = as.double(xb),
          cov_label = paste0(cov, "==", category)
        ))
      } else {
        # Otherwise, interpret the attribute as numeric
        x_num <- suppressWarnings(as.numeric(x))
        if (any(!is.finite(x_num))) {
          bad <- which(!is.finite(x_num))[1]
          ergm_Init_stop(
            sQuote(termname),
            ": numeric actor attribute contains NA/NaN/Inf. ",
            "Example: vertex=", vn[ia[bad]], ", value=", as.character(x[bad])
          )
        }
        return(list(
          x         = as.double(x_num),
          cov_label = cov
        ))
      }

    } else {
      # Case 2: direct numeric vector for actor covariate
      x_num <- suppressWarnings(as.numeric(cov))
      if (any(!is.finite(x_num)))
        ergm_Init_stop(sQuote(termname), ": vector 'cov' contains NA/NaN/Inf.")
      if (length(x_num) < n1)
        ergm_Init_stop(sQuote(termname), ": length(cov) < |A| = ", n1, ".")
      if (!is.null(category))
        ergm_Init_stop(sQuote(termname), ": 'category' does not apply when 'cov' is a direct numeric vector.")
      return(list(
        x         = as.double(x_num[seq_len(n1)]),
        cov_label = "cov"
      ))
    }
  }

  ax <- get_actor_cov(nw, cov = cov, category = category)
  x         <- ax$x
  cov_label <- ax$cov_label

  # ---------------------------------------------------------------------------
  # Normalize the size filter S (argument 'size')
  #   - NULL or empty => all group sizes (encoded as L = 0)
  #   - otherwise: distinct positive integers, sorted
  # ---------------------------------------------------------------------------
  S <- a$size
  if (is.null(S) || length(S) == 0L) {
    sizes <- integer(0)  # S = all group sizes
  } else {
    S <- unique(as.integer(S))
    if (any(!is.finite(S)) || any(S < 1L))
      ergm_Init_stop(sQuote(termname), ": 'size' must contain integers >= 1.")
    sizes <- sort(S)
  }

  # ---------------------------------------------------------------------------
  # Build INPUT_PARAM for the C change-statistic
  #   Layout (1-based indexing):
  #     [1]           = n1
  #     [2]           = L = length(sizes)
  #     [3..(2+L)]    = sizes (possibly L = 0)
  #     [3+L.. ]      = x[1..n1]
  # ---------------------------------------------------------------------------
  L <- length(sizes)
  inputs <- c(
    as.double(n1),
    as.double(L),
    as.double(sizes),
    as.double(x)
  )

  # ---------------------------------------------------------------------------
  # Coefficient name:
  #   - encodes the covariate label and the size filter S
  # ---------------------------------------------------------------------------
  pretty_sizes <- if (L == 0L) "all" else paste0("S{", paste(sizes, collapse = ","), "}")
  coef.names   <- paste0("cov_ingroup[", cov_label, "]_", pretty_sizes)

  # ---------------------------------------------------------------------------
  # Standard ERGM term specification
  #   - name must match C_CHANGESTAT_FN(c_cov_ingroup)
  # ---------------------------------------------------------------------------
  list(
    name         = "cov_ingroup",   # must match C_CHANGESTAT_FN(c_cov_ingroup)
    coef.names   = coef.names,
    inputs       = inputs,          # n1, L, sizes[L], x[n1]
    dependence   = TRUE,
    emptynwstats = 0
  )
}
