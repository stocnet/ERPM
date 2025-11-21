#' ERGM term: cov_diff (range over k-actor subsets within groups)
#'
#' @file InitErgmTerm.cov_diff.R
#'
#' @description
#' \code{cov_diff} is an ERGM term for bipartite networks that measures, for each
#' group, the dispersion of a numeric actor covariate over all \eqn{k}-actor
#' subsets within that group. The bipartite network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw %n% "bipartite"});
#'   \item a \emph{group mode} (the complementary side, representing groups).
#' }
#'
#' Each actor in the actor mode carries a numeric covariate value \eqn{x_i}.
#' For a fixed integer \eqn{k \ge 2}, and for each group, the term considers all
#' \eqn{k}-actor subsets of that group and computes, on each subset, the
#' max–min range of the covariate. The statistic can be used either in its raw
#' form (sum over all subsets in all groups) or in a by-group normalized form
#' (average range per group), controlled by the \code{normalized} argument.
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic, declared in the
#' compiled code under a symbol compatible with \code{name = "cov_diff"}. The R
#' initializer below:
#' \itemize{
#'   \item enforces that the network is bipartite via \code{nw %n% "bipartite"};
#'   \item extracts an actor-level covariate from a vertex attribute or from a
#'         literal vector, coercing it to numeric and failing fast on \code{NA};
#'   \item normalizes the \code{clique_size} argument into a single integer
#'         \eqn{k \ge 2};
#'   \item interprets \code{normalized} as either:
#'         \itemize{
#'           \item \code{FALSE}/\code{"none"}: raw sum over all \eqn{k}-subsets;
#'           \item \code{TRUE}/\code{"by_group"}: group-wise average over
#'                 \eqn{k}-subsets;
#'         }
#'   \item packs the covariate and scalar parameters into a compact
#'         \code{INPUT_PARAM} vector for the C layer.
#' }
#'
#' The ERGM infrastructure will call the C change-statistic whenever a toggle
#' affects an edge between an actor and a group. The C code recomputes the
#' local contribution of the affected group(s) to the \eqn{k}-subset range sum.
#'
#' @section Mathematical definition:
#' Let:
#' \itemize{
#'   \item \eqn{A} be the set of actor-mode nodes;
#'   \item \eqn{G} be the set of group-mode nodes;
#'   \item \eqn{B} the bipartite adjacency between actors and groups;
#'   \item \eqn{x_i \in \mathbb{R}} the numeric covariate value of actor
#'         \eqn{i \in A};
#'   \item \eqn{A_g = \{ i \in A : B_{ig}=1\}} the set of actors in group
#'         \eqn{g \in G};
#'   \item \eqn{n_g = |A_g|} the size of group \eqn{g};
#'   \item \eqn{\mathcal{C}_k(g)} the family of all \eqn{k}-element subsets of
#'         \eqn{A_g}.
#' }
#'
#' For each subset \eqn{S \in \mathcal{C}_k(g)}, define its range as
#' \deqn{
#'   R(S) = \max_{i \in S} x_i - \min_{i \in S} x_i.
#' }
#' The \emph{raw} statistic is
#' \deqn{
#'   T_k^{\text{raw}}(B;x)
#'   =
#'   \sum_{g \in G}
#'   \sum_{S \in \mathcal{C}_k(g)} R(S),
#' }
#' where groups with \eqn{n_g < k} contribute zero (they have no
#' \eqn{k}-subsets). The \emph{by-group normalized} statistic divides each group
#' contribution by the number of \eqn{k}-subsets in that group, i.e.
#' \deqn{
#'   T_k^{\text{by-group}}(B;x)
#'   =
#'   \sum_{g \in G}
#'     \mathbf{1}[\,n_g \ge k\,]
#'     \frac{1}{\binom{n_g}{k}}
#'     \sum_{S \in \mathcal{C}_k(g)} R(S).
#' }
#'
#' The argument \code{normalized} selects between these two variants.
#'
#' @section Usage:
#' Typical usage with {ergm} on a bipartite network \code{nw}:
#' \preformatted{
#'   # Raw (non-normalized) cov_diff on actor covariate "x_attr" with k = 2
#'   summary(nw ~ cov_diff(cov = "x_attr", clique_size = 2, normalized = FALSE))
#'
#'   # By-group normalized variant
#'   summary(nw ~ cov_diff(cov = "x_attr", clique_size = 3, normalized = TRUE))
#' }
#'
#' When using the ERPM wrapper on a partition-based workflow, the term can be
#' invoked indirectly as:
#' \preformatted{
#'   erpm(partition ~ cov_diff(cov = "x_attr", clique_size = 2, normalized = TRUE))
#' }
#' provided that the wrapper translates the partition into a bipartite network
#' with a consistent actor mode and group mode.
#'
#' @note
#' The network must be strictly bipartite:
#' \itemize{
#'   \item the actor mode size is given by \code{nw %n% "bipartite"} and must be
#'         a positive integer;
#'   \item the group mode is the complementary set of nodes;
#'   \item the covariate is read on the actor mode only.
#' }
#'
#' The covariate is coerced to numeric if possible (logical and integer types
#' are accepted), and any \code{NA} on the actor mode is rejected in a
#' fail-fast manner. The \code{clique_size} argument must be a single finite
#' numeric value that is rounded to an integer \eqn{k \ge 2}. The
#' \code{normalized} argument accepts logical values, the character strings
#' \code{"none"} and \code{"by_group"}, or a scalar numeric flag interpreted as
#' zero (raw) versus non-zero (by-group).
#'
#' Internally, \code{cov_diff} passes its parameters and covariate to the C
#' layer through \code{INPUT_PARAM} with layout
#' \code{c(n1, k, norm_flag, x[1:n1])}, where:
#' \itemize{
#'   \item \code{n1} is the actor-mode size;
#'   \item \code{k} is the subset size \eqn{k};
#'   \item \code{norm_flag} equals \code{0} for raw, \code{1} for by-group;
#'   \item \code{x} is the numeric actor covariate restricted to the actor mode.
#' }
#'
#' @examples
#' \dontrun{
#'   library(network)
#'   library(ergm)
#'
#'   # Build a small bipartite network: 4 actors, 2 groups
#'   n_actors <- 4
#'   n_groups <- 2
#'   n_total  <- n_actors + n_groups
#'
#'   # Adjacency: actors 1..4, groups 5..6
#'   adj <- matrix(0, n_total, n_total)
#'
#'   # Group 5 has actors 1,2,3
#'   adj[1, 5] <- adj[5, 1] <- 1
#'   adj[2, 5] <- adj[5, 2] <- 1
#'   adj[3, 5] <- adj[5, 3] <- 1
#'
#'   # Group 6 has actors 2,4
#'   adj[2, 6] <- adj[6, 2] <- 1
#'   adj[4, 6] <- adj[6, 4] <- 1
#'
#'   nw <- network(adj, directed = FALSE, matrix.type = "adjacency")
#'   nw %n% "bipartite" <- n_actors  # actor mode size
#'
#'   # Numeric covariate on the actor mode
#'   x_attr <- c(1.0, 3.0, 5.0, 10.0)
#'   set.vertex.attribute(nw, "x_attr", c(x_attr, rep(NA_real_, n_groups)))
#'
#'   # Example: k = 2 (all 2-actor subsets inside each group)
#'   summary(nw ~ cov_diff(cov = "x_attr", clique_size = 2))
#'
#'   # By-group normalized variant
#'   summary(nw ~ cov_diff(cov = "x_attr", clique_size = 2, normalized = TRUE))
#'
#'   # Fit a simple ERGM including the term
#'   fit <- ergm(nw ~ cov_diff(cov = "x_attr", clique_size = 2))
#'   summary(fit)
#' }
#'
#' @test
#' Self-tests for \code{cov_diff} construct small bipartite networks with
#' known group memberships and numeric covariates, and compare:
#' \itemize{
#'   \item the ERGM summary \code{summary(nw ~ cov_diff(...))};
#'   \item direct evaluations of
#'         \eqn{\sum_g \sum_{S \in \mathcal{C}_k(g)} (\max_{i \in S} x_i - \min_{i \in S} x_i)}
#'         (raw case);
#'   \item direct evaluations of
#'         \eqn{\sum_g \mathbf{1}[n_g \ge k] \binom{n_g}{k}^{-1}
#'               \sum_{S \in \mathcal{C}_k(g)} (\max_{i \in S} x_i - \min_{i \in S} x_i)}
#'         (by-group normalized case).
#' }
#' Additional checks verify that toggling an actor–group tie updates the
#' statistic by the expected local change in the range over all \eqn{k}-subsets
#' within the affected group, including cases where the group size crosses
#' the threshold \eqn{n_g = k}.
#'
#' @keywords ERGM term bipartite groups covariate range subsets
#' @md
#'
#' @export
InitErgmTerm.cov_diff <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_diff"

  # Debug helpers:
  # - dbg: logical flag, controlled by an R option;
  # - dbgcat(): emits prefixed debug messages when dbg is TRUE.
  dbg    <- isTRUE(getOption("ERPM.cov_diff.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[cov_diff][DEBUG]", ..., "\n", sep = "")

  # Run standard ERGM term checks and parse user arguments:
  # - enforce bipartite network;
  # - accept 'cov', 'clique_size', 'normalized' with flexible types;
  # - let {ergm} handle generic validations (missing args, etc.).
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov",                              "clique_size",            "normalized"),
    vartypes      = c("character,numeric,logical,vector", "numeric,integer",        "logical,character,numeric"),
    defaultvalues = list(NULL,                             2,                         FALSE),
    required      = c(TRUE,                                FALSE,                     FALSE)
  )

  # ----- 1) Actor-mode size n1 -----------------------------------------------
  # n1 is the number of actors, retrieved from the bipartite network attribute.
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L) stop(termname, ": réseau non biparti strict.")
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
    stop(termname, ": longueur de la covariée < n1.")

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
    stop(termname, ": la covariée doit être convertissable en numérique.")

  if (anyNA(cov_vec))
    stop(termname, ": NA non autorisé dans la covariée du mode acteurs.")

  # ----- 4) Subset size 'clique_size' = k >= 2 -------------------------------
  # 'clique_size' is interpreted as the subset size k used in the definition
  # of the statistic. It must be a single finite numeric value, rounded to
  # an integer and required to be at least 2.
  k_raw <- a$clique_size
  if (length(k_raw) != 1L || !is.numeric(k_raw))
    stop(termname, ": 'clique_size' doit être un scalaire numérique.")
  k <- as.integer(round(k_raw))
  if (!is.finite(k) || k < 2L)
    stop(termname, ": 'clique_size' doit être un entier >= 2.")
  dbgcat("clique_size k = ", k)

  # ----- 5) Normalization flag (raw vs by_group) -----------------------------
  # The 'normalized' argument can take several forms:
  # - NULL: treated as non-normalized (raw);
  # - logical: TRUE => by-group, FALSE => raw;
  # - character: "none" or "by_group";
  # - numeric: 0 => raw, non-zero => by-group.
  norm_raw   <- a$normalized
  norm_flag  <- 0L
  norm_label <- "raw"

  if (is.null(norm_raw)) {
    norm_flag  <- 0L
    norm_label <- "raw"
  } else if (is.logical(norm_raw)) {
    norm_flag  <- if (isTRUE(norm_raw)) 1L else 0L
    norm_label <- if (norm_flag == 1L) "by_group" else "raw"
  } else if (is.character(norm_raw) && length(norm_raw) == 1L) {
    nr <- match.arg(norm_raw, c("none", "by_group"))
    norm_flag  <- if (nr == "by_group") 1L else 0L
    norm_label <- nr
  } else if (is.numeric(norm_raw) && length(norm_raw) == 1L) {
    norm_flag  <- if (norm_raw != 0) 1L else 0L
    norm_label <- if (norm_flag == 1L) "by_group" else "raw"
  } else {
    stop(termname, ": 'normalized' doit être logique, numérique ou 'none'/'by_group'.")
  }

  dbgcat("normalized = ", norm_label, " (flag=", norm_flag, ")")

  # ----- 6) Coefficient name --------------------------------------------------
  # The coefficient name encodes the covariate label, the subset size k and
  # whether a by-group normalization is in use, for interpretability in model
  # summaries.
  coef.name <- sprintf(
    "cov_diff[%s]_k%d%s",
    cov_label, k,
    if (norm_flag == 1L) "_norm" else ""
  )
  dbgcat("coef.name = ", coef.name)

  # ----- 7) Build INPUT_PARAM for the C layer --------------------------------
  # Layout of INPUT_PARAM:
  #   [1]   = n1            (actor-mode size)
  #   [2]   = k             (subset size)
  #   [3]   = norm_flag     (0 raw, 1 by-group)
  #   [4..] = cov_vec[1:n1] (numeric covariate on the actor mode)
  inputs <- c(
    as.double(n1),
    as.double(k),
    as.double(norm_flag),
    as.double(cov_vec)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " k=", k, " norm_flag=", norm_flag,
         " | cov[1:6]=", paste(utils::head(signif(cov_vec, 5L), 6L), collapse = ","))

  # ----- 8) Return ERGM term specification -----------------------------------
  list(
    name         = "cov_diff",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, k, norm_flag, x[n1]
    dependence   = TRUE,
    minval       = 0,
    maxval       = Inf,
    emptynwstats = 0
  )
}