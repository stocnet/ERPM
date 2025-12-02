#' ERGM term: cov_fullmatch (group-level unanimity)
#'
#' @note InitErgmTerm.cov_fullmatch.R
#'
#' @description
#' \code{cov_fullmatch} is an ERGM term for bipartite networks that counts
#' groups whose actors are unanimously homogeneous on a categorical covariate.
#' The bipartite network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side that represents groups).
#' }
#'
#' Each actor in the actor mode carries a covariate value, which is mapped
#' internally to integer categories \eqn{1,\dots,K}. For each group node, we
#' consider the set of adjacent actors (its group of members) and check whether
#' all actors in that group share the same category. An optional size filter
#' restricts the set of group sizes that contribute to the statistic, and an
#' optional \code{category} argument targets a single covariate level.
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic, declared in the
#' compiled code under a symbol compatible with \code{name = "cov_fullmatch"}.
#' The R initializer below:
#' \itemize{
#'   \item enforces that the network is bipartite via \code{nw \%n\% "bipartite"};
#'   \item extracts and encodes a covariate on the actor mode as integer
#'         categories \eqn{1,\dots,K}, failing fast on \code{NA};
#'   \item normalizes the optional \code{size} filter into a sorted set
#'         of positive integers;
#'   \item encodes an optional targeted category into an integer code
#'         \eqn{\texttt{target} \in \{0,\dots,K\}} (with \code{0} = "no target");
#'   \item packs a compact \code{INPUT_PARAM} vector for the C layer.
#' }
#'
#' The ERGM infrastructure will call the C change-statistic whenever a toggle
#' affects an edge between an actor and a group. The C code reconstructs the
#' affected group memberships from the bipartite structure and updates the
#' count of unanimously homogeneous groups accordingly.
#'
#' @section Mathematical definition:
#' Let:
#' \itemize{
#'   \item \eqn{A} be the set of actor-mode nodes;
#'   \item \eqn{G} be the set of group-mode nodes;
#'   \item \eqn{B} the bipartite adjacency between actors and groups;
#'   \item \eqn{c_i \in \{1,\dots,K\}} be the covariate category of actor
#'         \eqn{i \in A};
#'   \item \eqn{A_g = \{ i \in A : B_{ig}=1\}} the set of actors in group
#'         \eqn{g \in G};
#'   \item \eqn{n_g = |A_g|} the size of group \eqn{g};
#'   \item \eqn{S} the set of allowed group sizes, derived from the
#'         \code{size} argument (\eqn{S = \mathbb{N}^*} if \code{size = NULL});
#'   \item \eqn{\kappa \in \{1,\dots,K\}} an optional targeted category
#'         (or "no target" if not provided).
#' }
#'
#' In the untargeted case (no \code{category}), the statistic is
#' \deqn{
#'   T(B;c)
#'   =
#'   \sum_{g \in G}
#'     \mathbf{1}[\,n_g \in S\,]
#'     \mathbf{1}\big[\exists r \in \{1,\dots,K\} :
#'         \forall i \in A_g,\; c_i = r\big],
#' }
#' i.e. the number of groups whose actors are unanimously in the same category,
#' restricted to group sizes in \eqn{S}.
#'
#' In the targeted case, when a category \eqn{\kappa} is specified, the statistic
#' becomes
#' \deqn{
#'   T_\kappa(B;c)
#'   =
#'   \sum_{g \in G}
#'     \mathbf{1}[\,n_g \in S\,]
#'     \mathbf{1}\big[\forall i \in A_g,\; c_i = \kappa\big],
#' }
#' i.e. the number of groups that are unanimously in the targeted category
#' \eqn{\kappa} and whose size belongs to \eqn{S}.
#'
#' @section Usage:
#' Typical usage with {ergm} on a bipartite network \code{nw}:
#' \preformatted{
#'   # Count all groups that are unanimously homogeneous on 'cov_attr',
#'   # with no size filter:
#'   summary(nw ~ cov_fullmatch(cov = "cov_attr"))
#'
#'   # Restrict to groups of size 3 or 4:
#'   summary(nw ~ cov_fullmatch(cov = "cov_attr", size = c(3, 4)))
#'
#'   # Only count groups unanimously in category "A":
#'   summary(nw ~ cov_fullmatch(cov = "cov_attr", category = "A"))
#' }
#'
#' When using the ERPM wrapper on a partition-based workflow, the term can be
#' invoked indirectly via:
#' \preformatted{
#'   erpm(partition ~ cov_fullmatch(cov = "cov_attr", size = c(2, 3)))
#' }
#' provided that the wrapper translates the partition into a bipartite
#' representation with a consistent actor mode and group mode.
#'
#' @note
#' The network must be strictly bipartite:
#' \itemize{
#'   \item the actor mode size is identified by \code{nw \%n\% "bipartite"} and
#'         must be a positive integer;
#'   \item the group mode is the complementary set of nodes;
#'   \item the term assumes that the bipartite structure encodes actor–group
#'         membership in a way consistent with the ERPM builder.
#' }
#' The covariate is read on the actor mode only and must contain no \code{NA}
#' values after restriction to the actor mode. Any \code{NA} is rejected in a
#' fail-fast manner. The optional \code{size} filter must contain strictly
#' positive integers if provided.
#'
#' Internally, \code{cov_fullmatch} encodes the covariate as a factor and then
#' maps its levels to integer codes \eqn{1,\dots,K}. The C change-statistic
#' receives these codes together with the size filter and, optionally, a
#' targeted category.
#'
#' @examples
#' \dontrun{
#'   library(network)
#'   library(ergm)
#'
#'   # Build a small bipartite network: 4 actors, 3 groups
#'   n_actors <- 4
#'   n_groups <- 3
#'   n_total  <- n_actors + n_groups
#'
#'   # Adjacency: actors 1..4, groups 5..7
#'   adj <- matrix(0, n_total, n_total)
#'   # Group 5 has actors 1 and 2
#'   adj[1, 5] <- adj[5, 1] <- 1
#'   adj[2, 5] <- adj[5, 2] <- 1
#'   # Group 6 has actors 3 and 4
#'   adj[3, 6] <- adj[6, 3] <- 1
#'   adj[4, 6] <- adj[6, 4] <- 1
#'   # Group 7 is empty
#'
#'   nw <- network(adj, directed = FALSE, matrix.type = "adjacency")
#'   nw \%n\% "bipartite" <- n_actors  # actor mode size
#'
#'   # Add a categorical covariate on the actor mode
#'   cov_attr <- c("A", "A", "B", "B")
#'   set.vertex.attribute(nw, "cov_attr", c(cov_attr, rep(NA, n_groups)))
#'
#'   # Group 5 is unanimously "A", group 6 is unanimously "B"
#'   summary(nw ~ cov_fullmatch(cov = "cov_attr"))
#'
#'   # Only count groups unanimously "A"
#'   summary(nw ~ cov_fullmatch(cov = "cov_attr", category = "A"))
#'
#'   # Only count groups of size exactly 2
#'   summary(nw ~ cov_fullmatch(cov = "cov_attr", size = 2))
#'
#'   # Fit a simple ERGM including the term
#'   fit <- ergm(nw ~ cov_fullmatch(cov = "cov_attr"))
#'   summary(fit)
#' }
#'
#' @section Tests:
#' Self-tests for \code{cov_fullmatch} construct small bipartite networks with
#' known group memberships and covariate patterns, and compare:
#' \itemize{
#'   \item the ERGM summary \code{summary(nw ~ cov_fullmatch(...))};
#'   \item direct counts of groups that are unanimously in a given category,
#'         with and without size filters.
#' }
#' Additional checks verify that toggling actor–group ties updates the statistic
#' by the expected local increment, including edge cases such as empty groups,
#' groups of size 1, and groups whose covariate composition switches between
#' homogeneous and heterogeneous.
#'
#' @keywords ERGM term bipartite groups covariate homogeneity
#' @md
#'
#' @export
InitErgmTerm.cov_fullmatch <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_fullmatch"

  # Debug helpers (controlled by an R option):
  # - dbg: boolean flag;
  # - dbgcat(): prefixed debug output when dbg is TRUE.
  dbg <- isTRUE(getOption("ERPM.cov_fullmatch.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[cov_fullmatch][DEBUG]", ..., "\n", sep = "")

  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov",                              "size",             "category"),
    vartypes      = c("character,numeric,logical,vector", "numeric,integer",  "character,numeric,logical"),
    defaultvalues = list(NULL,                             NULL,               NULL),
    required      = c(TRUE,                                FALSE,              FALSE)
  )

  # ----- 1) Actor-mode size n1 -----------------------------------------------
  # n1 is the number of actors, as stored in the bipartite network attribute.
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L) stop(termname, ": réseau non biparti strict.")
  dbgcat("n1 = ", n1)

  # ----- 2) Extract the actor covariate vector (length >= n1) -----------------
  # The covariate can be given as:
  # - the name of a vertex attribute (preferred);
  # - a literal vector. In all cases we keep only the first n1 entries for actors.
  cov_raw <- a$cov
  if (is.character(cov_raw) && length(cov_raw) == 1L) {
    cov_vec   <- network::get.vertex.attribute(nw, cov_raw)
    cov_label <- cov_raw
    if (is.null(cov_vec)) stop(termname, ": attribut inexistant: ", sQuote(cov_raw), ".")
    dbgcat("cov source = vertex attribute ", sQuote(cov_label))
  } else {
    cov_vec   <- cov_raw
    cov_label <- "cov"
    dbgcat("cov source = vector literal")
  }
  if (length(cov_vec) < n1) stop(termname, ": longueur de la covariée < n1.")
  cov_vec <- cov_vec[seq_len(n1)]
  dbgcat("cov length = ", length(cov_vec), " | head = ", paste(utils::head(as.character(cov_vec), 6L), collapse = ","))

  # ----- 3) Fail-fast on missing values ---------------------------------------
  # Any NA on the actor mode is rejected, since the change-statistic expects
  # well-defined categories for all actors.
  if (anyNA(cov_vec)) stop(termname, ": NA non autorisé dans la covariée du mode acteurs.")

  # ----- 4) Normalize to integer category codes 1..K --------------------------
  # We map the actor covariate to a factor, then to integer codes:
  # - levels_f: distinct covariate values;
  # - K: number of categories;
  # - cats: integer codes in {1..K} for each actor.
  if (is.logical(cov_vec)) cov_vec <- as.integer(cov_vec)
  f        <- factor(cov_vec)           # no NA at this stage
  levels_f <- levels(f)
  K        <- length(levels_f)
  if (K == 0L) stop(termname, ": aucune modalité valide trouvée.")
  cats     <- as.integer(f)             # 1..K
  dbgcat("K = ", K, " | levels = {", paste(levels_f, collapse = ","), "}")

  # ----- 5) Targeted category -> 'target' code (0 if not provided) ------------
  # If a single 'category' is given and matches one of the factor levels, we
  # encode it as an integer code in 1..K. Otherwise, target=0 means "no target".
  has_category <- {
    x <- a$category
    !is.null(x) && length(x) == 1L && !is.na(x) && nzchar(as.character(x))
  }
  dbgcat("has_category = ", has_category)

  target <- 0L
  if (has_category) {
    cat_val <- a$category
    if (is.logical(cat_val)) cat_val <- as.integer(cat_val)
    ix <- match(as.character(cat_val), levels_f, nomatch = 0L)
    if (ix == 0L) {
      warning(termname, ": 'category' non trouvée dans les modalités; cible ignorée.", call. = FALSE)
      dbgcat("category ", sQuote(as.character(a$category)), " not found -> target=0 (ignored)")
    } else {
      target <- as.integer(ix)
      cov_label <- paste0(cov_label, "==", levels_f[target])
      dbgcat("category target = ", target, " -> label suffix = ", levels_f[target])
    }
  }

  # ----- 6) Size filter S (argument 'size') -----------------------------------
  # The 'size' argument defines the set S of allowed group sizes:
  # - NULL => all group sizes are allowed (L=0, empty size vector);
  # - numeric => converted to a sorted, unique vector of positive integers.
  sizes <- a$size
  if (is.null(sizes)) {
    L <- 0L
    sizes_vec <- numeric(0)
    dbgcat("size filter = <ALL> (L=0)")
  } else {
    if (!is.numeric(sizes))
      stop("cov_fullmatch: 'size' doit être numérique.")
    if (length(sizes) == 0L)
      stop("cov_fullmatch: 'size' vide (integer(0)) interdit. Utilisez NULL pour toutes tailles.")
    sizes <- as.integer(round(sizes))
    if (any(!is.finite(sizes)) || any(sizes <= 0L))
      stop("cov_fullmatch: 'size' doit contenir des entiers positifs.")
    sizes <- sort(unique(sizes))
    L <- length(sizes)
    sizes_vec <- as.double(sizes)
    dbgcat("size filter = {", paste(sizes, collapse = ","), "} (L=", L, ")")
  }

  # ----- 7) Coefficient name (kept backward compatible) -----------------------
  # We encode both the covariate label and the size filter in the coefficient
  # name, keeping the existing convention for compatibility.
  size_label <- if (L > 0L) sprintf("_S{%s}", paste(sizes, collapse = ",")) else "_all"
  coef.name  <- sprintf("cov_fullmatch[%s]%s", cov_label, size_label)
  dbgcat("coef.name = ", coef.name)

  # ----- 8) Build INPUT_PARAM for the C layer ---------------------------------
  # Layout of INPUT_PARAM:
  #   [1]   = n1                 (actor-mode size)
  #   [2]   = L                  (number of size values)
  #   [3..] = sizes_vec          (L entries, may be empty)
  #   [...] = K                  (number of categories)
  #   [...] = target             (0 if no target, or category code in 1..K)
  #   [...] = cats[1..n1]        (category code for each actor)
  inputs <- c(
    as.double(n1),
    as.double(L),
    sizes_vec,
    as.double(K),
    as.double(target),
    as.double(cats)
  )
  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " L=", L, " K=", K, " target=", target, " | cats[1:6]=",
         paste(utils::head(as.integer(cats), 6L), collapse = ","))

  # ----- 9) Return the ERGM term specification --------------------------------
  # This structure is what {ergm} expects from an InitErgmTerm.* initializer.
  list(
    name         = "cov_fullmatch",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, L, sizes[L], K, target, cats[n1]
    dependence   = TRUE,
    minval       = 0,
    maxval       = Inf,
    emptynwstats = 0
  )
}
