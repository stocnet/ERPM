#' ERPM wrapper: translate ERPM formulas to {ergm} calls and optionally fit
#' @name erpm_wrapper
#' @note erpm_wrapper.R
#' @description
#' This module provides a thin wrapper around {ergm} to:
#' \enumerate{
#'   \item Build a bipartite network from a partition, with optional node and dyadic inputs.
#'   \item Translate ERPM RHS terms (e.g., \code{groups}, \code{cov_match}, \code{cliques})
#'         into {ergm} terms, with optional encapsulations (\code{Proj1}, \code{B}).
#'   \item Compose a standard call \code{ergm(nw ~ <translated RHS>, constraints = ~ b1part, ...)}.
#'   \item Either return the call (dry-run) or evaluate it and return the fitted model.
#' }
#'
#' The file also includes:
#' \itemize{
#'   \item argument normalizers for ERPM-specific terms such as \code{groups(...)}; 
#'   \item internal validators for node data and dyadic covariates;
#'   \item a builder for bipartite networks from partitions;
#'   \item a single-term translator used by the ERPM→{ergm} translation layer.
#' }
#'
#' All user-facing console messages are in English for consistency with {ergm}.
#'
#' @note
#' This wrapper assumes that the user has attached the {ergm} and {network} packages
#' (or that they are available in the search path) and that the \code{b1part}
#' constraint is meaningful for the constructed network.
#'
#' @examples
#' \dontrun{
#'   # Basic usage from a partition vector
#'   partition <- c(1, 1, 2, 2, 3)
#'   fit <- erpm(partition ~ groups + cliques(3), estimate = "MLE")
#'
#'   # Dry-run: only inspect the translated {ergm} call
#'   call_only <- erpm(partition ~ groups(2), eval_call = FALSE)
#'   print(call_only)
#' }
#'
#' @note
#' The main wrapper is exercised in self-tests and MWEs under \code{scripts/test}
#' by comparing ERPM-based fits to direct {ergm} calls on constructed bipartite networks.
#'
#' @keywords ERPM ERGM wrapper bipartite translation

if (!exists(".__erpm_wrapper_loaded", envir = .GlobalEnv)) {

  # ============================================================================
  # 1. Small generic helpers
  # ============================================================================

  #' Null-coalescing operator
  #'
  #' Return `a` if it is not `NULL`, otherwise return `b`.
  #' Used throughout the wrapper to avoid nested `if` / `else` for default values.
  #'
  #' @param a any R object, possibly `NULL`
  #' @param b any R object used as fallback
  #' @return `a` if not `NULL`, otherwise `b`
  #' @examples
  #' 1 %||% 2
  #' NULL %||% 2
  #'
  #' @note
  #' This is a small utility for readability; it does not introduce any new
  #' semantics beyond a standard `if (is.null(a)) b else a`.
  #'
  #' @note
  #' See tests where wrapper options are normalized using `%||%`.
  #' @keywords internal
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  #' One-line deparser for logging
  #'
  #' Convert an object to a single-line string suitable for compact logging
  #' and error messages.
  #'
  #' @param x an R object
  #' @return character scalar
  #' @examples
  #' .oneline(quote(a + b + c))
  #'
  #' @note
  #' Width cutoff is fixed at 500L to reduce line breaks when deparsing formulas.
  #'
  #' @note
  #' Used in logging inside \code{erpm()} for user-facing messages.
  #' @keywords internal
  .oneline <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")

  #' Remove all whitespace from a string
  #'
  #' Replace every run of whitespace (spaces, tabs, newlines) by nothing.
  #' This is used to compactify formulas and calls for logging.
  #'
  #' @param s character scalar
  #' @return character scalar without whitespace
  #' @examples
  #' .tight(" a  +  b \n c ")
  #'
  #' @note
  #' Unlike a simple trimming, this removes \emph{all} whitespace tokens.
  #' This matches the current behavior of the ERPM logging format.
  #'
  #' @note
  #' Checked implicitly by tests that compare log messages in error paths.
  #' @keywords internal
  .tight <- function(s) gsub("\\s+", "", s)

  #' Minimal bipartite network test
  #'
  #' Check if an object looks like a bipartite {network} by testing
  #' the presence and finiteness of the `%n% "bipartite"` attribute.
  #'
  #' @param x any R object
  #' @return logical, `TRUE` if `x` is a bipartite network
  #' @examples
  #' \dontrun{
  #'   is_bipartite_network(nw)
  #' }
  #'
  #' @note
  #' This helper is currently not used inside this file but kept for potential
  #' future sanity checks on user-supplied networks.
  #'
  #' @note
  #' May be exercised in higher-level tests that inspect pre-built networks.
  #' @keywords internal
  # is_bipartite_network <- function(x) {
  #   inherits(x, "network") &&
  #     !is.null(x %n% "bipartite") &&
  #     is.finite(x %n% "bipartite")
  # }

  # ============================================================================
  # 2. Node / dyad validation helpers
  # ============================================================================

  #' Guess a label column in a node data frame
  #'
  #' Search for a column to be used as node labels. Preference order:
  #' \code{prefer}, then a small set of common aliases.
  #'
  #' @param nodes a data.frame of node attributes
  #' @param prefer character scalar, preferred label column name
  #' @return character scalar column name or `NULL` if none found
  #' @examples
  #' df <- data.frame(id = 1:3, name = c("A","B","C"))
  #' .get_label_col(df)
  #'
  #' @note
  #' The aliases include "label", "nom", "name", and "id" to be robust
  #' to typical naming patterns in user data.
  #'
  #' @note
  #' `build_bipartite_from_inputs()` indirectly tests this helper by using
  #' various node input configurations.
  #' @keywords internal
  .get_label_col <- function(nodes, prefer = "label") {
    stopifnot(is.data.frame(nodes))
    nms <- trimws(names(nodes))
    aliases <- c(prefer, "label", "nom", "name", "id")
    # Avoid duplicate "label" entries when prefer = "label"
    hit <- intersect(aliases, nms)
    if (length(hit)) hit[1L] else NULL
  }

  #' Validate a node data frame
  #'
  #' Requirements:
  #' \itemize{
  #'   \item `nodes` must be a data.frame;
  #'   \item if a label-like column exists, it must not contain duplicates.
  #' }
  #'
  #' @param nodes a data.frame of node attributes
  #' @return invisibly `TRUE` on success; otherwise raises an error
  #' @examples
  #' df <- data.frame(label = c("A","B","C"))
  #' .check_nodes_df(df)
  #'
  #' @note
  #' The label column is chosen using \code{.get_label_col()} and may be NULL
  #' if no reasonable candidate is found.
  #'
  #' @note
  #' Used by \code{build_bipartite_from_inputs()} when validating node inputs.
  #' @keywords internal
  .check_nodes_df <- function(nodes) {
    stopifnot(is.data.frame(nodes))
    lab <- .get_label_col(nodes)  # may be NULL if no obvious label column
    if (!is.null(lab) && anyDuplicated(nodes[[lab]]))
      stop(sprintf("nodes$%s contains duplicates.", lab))
    invisible(TRUE)
  }

  #' Validate dyadic n×n matrices
  #'
  #' Ensure that each entry of a list of dyadic matrices is:
  #' \itemize{
  #'   \item a matrix;
  #'   \item square n×n;
  #'   \item if row/colnames are present, they match the actor label order.
  #' }
  #'
  #' These matrices are later stored on the network as `%n%` attributes.
  #'
  #' @param dyads named list of matrices
  #' @param n integer number of actors (mode 1)
  #' @param labels character vector of actor labels, used for name checks
  #' @return invisibly `TRUE` on success; otherwise raises an error
  #' @examples
  #' M <- diag(3)
  #' .check_dyads(list(ex = M), n = 3, labels = c("A","B","C"))
  #'
  #' @note
  #' This helper does not attach the matrices; it only validates their shape
  #' and naming. Attachment is performed in \code{build_bipartite_from_inputs()}.
  #'
  #' @note
  #' Tested indirectly by ERPM self-tests using dyadic covariates.
  #' @keywords internal
  .check_dyads <- function(dyads, n, labels) {
    if (length(dyads) == 0L) return(invisible(TRUE))
    stopifnot(is.list(dyads))
    for (nm in names(dyads)) {
      M <- dyads[[nm]]
      stopifnot(is.matrix(M), nrow(M) == n, ncol(M) == n)
      # If names are present, enforce exact matching to actor label order.
      if (!is.null(rownames(M)) && !is.null(colnames(M))) {
        if (!identical(rownames(M), labels) || !identical(colnames(M), labels))
          stop(sprintf("dyads['%s']: row/colnames must match the actor order.", nm))
      }
    }
    invisible(TRUE)
  }

  # ============================================================================
  # 3. Bipartite builder and helper
  # ============================================================================

  #' Helper for bipartite builder: stop with usage message
  #'
  #' Centralize the error and usage message for
  #' \code{build_bipartite_from_inputs()}, so that all validation failures
  #' provide the same minimal reproducible example.
  #'
  #' @param msg character scalar, specific error message
  #' @return this function always calls [base::stop()]
  #' @examples
  #' \dontrun{
  #'   .stop_build_bipartite("partition must be non-empty")
  #' }
  #'
  #' @note
  #' This function is not exported. It is used internally by
  #' \code{build_bipartite_from_inputs()}.
  #'
  #' @note
  #' The content of the usage message is not asserted, but calls from
  #' invalid inputs should raise an error including the provided \code{msg}.
  #' @keywords internal
  .stop_build_bipartite <- function(msg) {
    usage <- paste(
      "[ERPM] build_bipartite_from_inputs usage:",
      "",
      "  build_bipartite_from_inputs(partition, nodes = NULL, dyads = list())",
      "",
      "  # Minimal example:",
      "  partition <- c(1, 1, 2, 2, 2, 3)",
      "  nodes <- data.frame(",
      "    label  = c('A','B','C','D','E','F'),",
      "    gender = c(1, 1, 2, 1, 2, 2),",
      "    age    = c(20, 22, 25, 30, 30, 31)",
      "  )",
      "  friendship <- matrix(c(",
      "    0,1,1,1,0,0,",
      "    1,0,0,0,1,0,",
      "    1,0,0,0,1,0,",
      "    1,0,0,0,0,0,",
      "    0,1,1,0,0,1,",
      "    0,0,0,0,1,0",
      "  ), 6, 6, byrow = TRUE)",
      "  dyads <- list(friendship = friendship)",
      "",
      "  built <- build_bipartite_from_inputs(partition, nodes = nodes, dyads = dyads)",
      "  nw    <- built$network",
      sep = "\n"
    )
    stop(paste0(msg, "\n\n", usage), call. = FALSE)
  }

  #' Build a padded bipartite network from a partition
  #'
  #' Construct a bipartite {network} object from:
  #' \itemize{
  #'   \item a partition vector of length \eqn{n} (mode-1 actors),
  #'   \item optional node attributes in a data.frame,
  #'   \item optional dyadic n×n matrices, attached as `%n%` attributes.
  #' }
  #'
  #' The function pads the number of groups to \eqn{G = n} to maximize ERGM
  #' variability (unused groups become empty mode-2 vertices).
  #'
  #' @param partition integer vector of group ids, length n, values in 1..G
  #' @param nodes optional data.frame with a label column (auto-detected)
  #' @param dyads optional named list of n×n matrices to attach as `%n%` attributes
  #' @return list with components:
  #'   \itemize{
  #'     \item \code{network}: the bipartite {network} object;
  #'     \item \code{partition}: the original partition vector;
  #'     \item \code{actor_labels}: character vector of actor labels;
  #'     \item \code{group_labels}: character vector of group labels.
  #'   }
  #' @examples
  #' \dontrun{
  #'   partition <- c(1, 1, 2, 2, 3)
  #'   nodes <- data.frame(
  #'     label = paste0("A", 1:5),
  #'     age   = c(20, 21, 22, 23, 24)
  #'   )
  #'   friendship <- matrix(0, 5, 5)
  #'   built <- build_bipartite_from_inputs(partition, nodes = nodes,
  #'                                       dyads = list(friendship = friendship))
  #'   nw <- built$network
  #' }
  #'
  #' @note
  #' This builder is used by \code{erpm()} when the LHS of the formula is
  #' a partition vector rather than a pre-built network.
  #' The bipartite attribute is set to \code{n} so that \code{constraints = ~ b1part}
  #' is meaningful for {ergm}.
  #'
  #' @note
  #' Self-tests construct bipartite networks from partitions and compare their
  #' summaries to reference partitions.
  #' @keywords internal
  build_bipartite_from_inputs <- function(partition = NULL,
                                          nodes     = NULL,
                                          dyads     = list()) {
    # -- Checks on partition ----------------------------------------------------
    if (is.null(partition) || !is.atomic(partition) || length(partition) < 1L) {
      .stop_build_bipartite("partition must be a non-empty atomic vector.")
    }
    if (any(!is.finite(partition)) || any(partition < 1)) {
      .stop_build_bipartite("partition must contain finite positive integers.")
    }
    n <- length(partition)

    # -- Nodes / labels ---------------------------------------------------------
    if (is.null(nodes)) {
      labels <- sprintf("A%d", seq_len(n))
      nodes  <- data.frame(label = labels, stringsAsFactors = FALSE)
    } else {
      # Wrap .check_nodes_df errors in a more informative message
      ok_nodes <- try(.check_nodes_df(nodes), silent = TRUE)
      if (inherits(ok_nodes, "try-error")) {
        .stop_build_bipartite(
          paste0("Invalid `nodes` data.frame: ", conditionMessage(attr(ok_nodes, "condition")))
        )
      }
      if (nrow(nodes) != n) {
        .stop_build_bipartite("nrow(nodes) must equal length(partition).")
      }
      labcol <- .get_label_col(nodes)
      if (is.null(labcol)) {
        labels <- sprintf("A%d", seq_len(n))
        nodes$label <- labels
      } else {
        labels <- as.character(nodes[[labcol]])
        if (anyNA(labels) || any(!nzchar(labels)))
          .stop_build_bipartite("empty/NA labels are not allowed in `nodes`.")
        if (anyDuplicated(labels))
          .stop_build_bipartite("duplicate labels are not allowed in `nodes`.")
        if (!("label" %in% names(nodes))) nodes$label <- labels
      }
    }

    # -- Padded groups: G = n ---------------------------------------------------
    G_obs <- max(partition, na.rm = TRUE)
    if (G_obs > n) {
      .stop_build_bipartite("max(partition) cannot exceed n when padding with G = n.")
    }
    G <- n  # force as many groups as actors

    # Validate dyadic matrices
    ok_dyads <- try(.check_dyads(dyads, n, labels), silent = TRUE)
    if (inherits(ok_dyads, "try-error")) {
      .stop_build_bipartite(
        paste0("Invalid `dyads` list: ", conditionMessage(attr(ok_dyads, "condition")))
      )
    }

    # -- Names and adjacency ----------------------------------------------------
    g_names <- sprintf("G%d", seq_len(G))
    all_v   <- c(labels, g_names)

    adj <- matrix(0L, n + G, n + G, dimnames = list(all_v, all_v))
    # actor indices: 1..n ; group indices: n + (1..G)
    idx_actor <- seq_len(n)
    idx_group <- n + as.integer(partition)

    # Add bipartite edges for observed membership only; extra groups remain empty.
    adj[cbind(idx_actor, idx_group)] <- 1L
    adj[cbind(idx_group, idx_actor)] <- 1L

    # -- Build network object ---------------------------------------------------
    nw <- network::network(adj, directed = FALSE, matrix.type = "adjacency")
    network::set.network.attribute(nw, "bipartite", n)
    network::set.vertex.attribute(nw, "vertex.names", all_v)

    # Push actor attributes; pad groups with NA
    if (ncol(nodes) > 1L) {
      for (a in setdiff(names(nodes), c("label"))) {
        vals <- nodes[[a]]
        network::set.vertex.attribute(nw, a, c(vals, rep(NA, G)))
      }
    }

    # Attach dyadic n×n matrices as %n% attributes, with enforced dimnames order
    if (length(dyads)) {
      for (nm in names(dyads)) {
        M <- dyads[[nm]]
        dimnames(M) <- list(labels, labels)
        nw %n% nm <- M
      }
    }

    list(
      network      = nw,
      partition    = partition,
      actor_labels = labels,
      group_labels = g_names
    )
  }

  # ============================================================================
  # 4. ERPM term normalizers and splitters
  # ============================================================================

  #' Normalize `groups(...)` arguments into a (from, to) interval
  #'
  #' Supported ERPM signatures:
  #' \itemize{
  #'   \item \code{groups}
  #'   \item \code{groups(k)} or \code{groups(size = k)}
  #'   \item \code{groups(from = ..., to = ...)}
  #' }
  #'
  #' Rules:
  #' \itemize{
  #'   \item \code{from} must be a finite integer \eqn{\ge 0};
  #'   \item \code{to} must be a finite integer strictly greater than `from` or `Inf`;
  #'   \item `Inf` is allowed for `to`, but not for `from`.
  #' }
  #'
  #' @param args_list list of arguments as captured from the call
  #' @return list with elements \code{from} (integer) and \code{to} (integer or `quote(Inf)`)
  #' @examples
  #' .normalize_groups_args(list())
  #' .normalize_groups_args(list(3))
  #' .normalize_groups_args(list(from = 2, to = 5))
  #'
  #' @note
  #' This normalizer is used exclusively by \code{.translate_one_term()} for
  #' the special-case translation of \code{groups(...)} into \code{b2degrange}.
  #'
  #' @note
  #' Covered indirectly by tests that exercise ERPM formulas with \code{groups()}.
  #' @keywords internal
  .normalize_groups_args <- function(args_list) {
    nm <- names(args_list)

    # -- helpers ---------------------------------------------------------------
    .deparse1 <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
    .eval_num_scalar <- function(x, env = parent.frame()) {
      if (is.numeric(x) && length(x) == 1L) return(x)
      vx <- try(eval(x, envir = env), silent = TRUE)
      if (inherits(vx, "try-error")) return(NULL)
      if (is.numeric(vx) && length(vx) == 1L) return(vx)
      NULL
    }
    .as_from_int <- function(x, env = parent.frame()) {
      if (is.symbol(x) && identical(x, as.name("Inf")))
        stop("groups(from,to): 'from' cannot be Inf.")
      v <- .eval_num_scalar(x, env)
      if (is.null(v) || !is.finite(v))
        stop(sprintf("groups(from): must be a finite integer >= 0. Got: %s", .deparse1(x)))
      iv <- as.integer(round(v))
      if (!isTRUE(all.equal(v, iv))) stop(sprintf("groups(from): integer required. Got: %s", format(v)))
      if (iv < 0L) stop("groups(from): must be >= 0.")
      iv
    }
    .as_to_val <- function(x, env = parent.frame()) {
      # Accept both symbol Inf and numeric Inf.
      if ((is.symbol(x) && identical(x, as.name("Inf"))) ||
          (is.numeric(x) && length(x) == 1L && is.infinite(x))) {
        return(quote(Inf))
      }
      v <- .eval_num_scalar(x, env)
      if (is.null(v) || !is.finite(v))
        stop(sprintf("groups(to): must be a finite integer or Inf. Got: %s", .deparse1(x)))
      iv <- as.integer(round(v))
      if (!isTRUE(all.equal(v, iv))) stop(sprintf("groups(to): integer or Inf required. Got: %s", format(v)))
      iv
    }
    # -------------------------------------------------------------------------

    # No args: groups ≡ [1, Inf)
    if (length(args_list) == 0L) return(list(from = 1L, to = quote(Inf)))

    # Alias support: groups(size = k)
    if (!is.null(nm) && "size" %in% nm && !("from" %in% nm) && !("to" %in% nm)) {
      args_list <- list(args_list[["size"]]); nm <- NULL
    }

    # Single positional: groups(k) ≡ [k, k+1)
    if (length(args_list) == 1L && (is.null(nm) || isTRUE(nm[1L] == ""))) {
      k <- .as_from_int(args_list[[1L]], parent.frame())
      return(list(from = k, to = k + 1L))
    }

    # Named pair: groups(from=..., to=...)
    if (!is.null(nm) && all(c("from","to") %in% nm)) {
      from <- .as_from_int(args_list[["from"]], parent.frame())
      to   <- .as_to_val  (args_list[["to"]],   parent.frame())

      # Resolve numeric to for check; keep quote(Inf) if Inf.
      to_num <- if (is.language(to)) Inf else to
      if (!(is.infinite(to_num) || (is.finite(to_num) && to_num > from))) {
        to_str <- if (is.language(to)) "Inf" else as.character(to_num)
        stop(sprintf("groups(from,to): requires 'from' < 'to'. Got: from=%d, to=%s", from, to_str))
      }
      return(list(from = from, to = to))
    }

    stop("groups(): use `groups`, `groups(k)`/`groups(size=k)`, or `groups(from=..,to=..)`. ")
  }

  #' Split a sum-of-terms RHS recursively
  #'
  #' Decompose a RHS expression of the form
  #' \code{edges + triangles + nodematch("grp")} into a list of call objects,
  #' each representing a single term.
  #'
  #' @param expr a language object representing the RHS expression
  #' @return list of call objects (each element is either a call or an atomic term)
  #' @examples
  #' expr <- quote(edges + triangles + nodematch("grp"))
  #' .split_sum_terms(expr)
  #'
  #' @note
  #' This is a syntactic splitter. It does not interpret the meaning of the
  #' underlying terms; it only follows the binary \code{+} structure.
  #'
  #' @note
  #' Used by \code{erpm()} to translate each term independently.
  #' @keywords internal
  .split_sum_terms <- function(expr) {
    if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
      return(c(.split_sum_terms(expr[[2]]),
               .split_sum_terms(expr[[3]])))
    }
    list(expr)
  }

  # ============================================================================
  # 5. Single-term translator ERPM → {ergm}
  # ============================================================================

  #' Translate a single ERPM term into an {ergm} term
  #'
  #' This function takes a single RHS term from an ERPM formula and returns
  #' a corresponding {ergm} call. It supports:
  #' \itemize{
  #'   \item \code{groups(...)} → \code{b2degrange(from, to)};
  #'   \item \code{cliques(...)}: normalize aliases and keep as ERPM term name;
  #'   \item generic renames via \code{rename_map} and optional wrapping inside
  #'         \code{Proj1(~ term)} or \code{B(~ term, form = "nonzero")}.
  #' }
  #'
  #' Bare symbols are turned into zero-argument calls to unify processing.
  #'
  #' @param term_call a call object representing a single RHS term
  #' @param rename_map named character vector mapping old → new term names
  #' @param wrap_proj1 character vector of term names to wrap in `Proj1(~ term)`
  #' @param wrap_B character vector of term names to wrap in `B(~ term, form="nonzero")`
  #' @return a translated call suitable for {ergm}
  #' @examples
  #' .translate_one_term(quote(groups(3)), rename_map = c())
  #' .translate_one_term(quote(cliques(k = 3, normalized = TRUE)), rename_map = c())
  #'
  #' @note
  #' In the current version, \code{rename_map}, \code{wrap_proj1}, and \code{wrap_B}
  #' are typically empty, but they are kept for forward compatibility with other
  #' ERPM terms that may be translated into standard ERGM effects.
  #'
  #' @note
  #' Indirectly tested by ERPM self-tests that use `groups()` and `cliques()`
  #' in ERPM formulas.
  #' @keywords internal
  .translate_one_term <- function(term_call,
                                  rename_map,
                                  wrap_proj1 = character(),
                                  wrap_B     = character()) {
    # Turn bare symbols into zero-arg calls to unify processing.
    if (is.symbol(term_call)) term_call <- as.call(list(term_call))
    if (!is.call(term_call)) return(term_call)

    fun_sym   <- term_call[[1L]]
    args_list <- as.list(term_call)[-1L]

    # Special case: `groups(...)`
    if (is.symbol(fun_sym) && identical(fun_sym, as.name("groups"))) {
      gt <- .normalize_groups_args(args_list)
      return(as.call(list(as.name("b2degrange"), from = gt$from, to = gt$to)))
    }

    # Special case: `cliques(...)` (normalize aliases, keep as ERPM term name)
    #
    # Important:
    # - Do NOT recompute or override `normalized` here (we forward the user value).
    # - Let `InitErgmTerm.cliques()` handle defaults and validation for `k`
    #   and `normalized`.
    if (is.symbol(fun_sym) && identical(fun_sym, as.name("cliques"))) {
      # If there are no arguments, or only positional ones, we can just pass
      # the call through: InitErgmTerm.cliques() already supports:
      #   cliques()
      #   cliques(3)
      #   cliques(k = 3)
      #   cliques(clique_size = 3)
      if (length(args_list) == 0L) {
        return(term_call)
      }

      al <- as.pairlist(args_list)

      # If the user used the alias `clique_size` but not `k`, rename it to `k`
      # so that the underlying ERGM term sees a canonical argument name.
      if (!is.null(names(al))) {
        if (!is.null(al$clique_size) && is.null(al$k)) {
          al$k <- al$clique_size
          al$clique_size <- NULL
        }
      }

      # Rebuild the call with possibly renamed arguments, but without touching
      # `normalized` (or any other argument).
      return(as.call(c(as.name("cliques"), as.list(al))))
    }

    # Generic rename, then optional wrappers.
    fname <- if (is.symbol(fun_sym)) as.character(fun_sym) else deparse(fun_sym)[1L]
    fname <- if (fname %in% names(rename_map)) rename_map[[fname]] else fname
    out <- as.call(c(as.name(fname), args_list))
    if (fname %in% wrap_B)     out <- call("B",     call("~", out), form = "nonzero")
    if (fname %in% wrap_proj1) out <- call("Proj1", call("~", out))
    out
  }

  # ============================================================================
  # 6. Main ERPM → {ergm} wrapper
  # ============================================================================

  #' ERPM main wrapper: translate and optionally fit with {ergm}
  #'
  #' This function:
  #' \enumerate{
  #'   \item Interprets the LHS of a formula as either a partition vector
  #'         or a pre-built bipartite {network};
  #'   \item Builds a bipartite network from a partition if needed;
  #'   \item Splits and translates ERPM RHS terms into {ergm} terms;
  #'   \item Constructs an \code{ergm()} call with constraint \code{~ b1part};
  #'   \item Either returns the unevaluated call (dry-run) or evaluates it.
  #' }
  #'
  #' @param formula `lhs ~ <ERPM terms>`, where `lhs` is either a partition
  #'        vector or an already built bipartite `network` object.
  #' @param eval_call logical. If TRUE evaluate the {ergm} call, else return it.
  #' @param verbose logical. If TRUE print wrapper translation log and options.
  #'        When explicitly set in the `erpm()` call, its value is also forwarded
  #'        to \code{ergm()} as the \code{verbose=} argument.
  #' @param estimate character or NULL.
  #'        If non-NULL, its value is forwarded to \code{ergm()} as the
  #'        \code{estimate=} argument, after a light normalization that maps
  #'        `"MCMLE"` to `"MLE"`. If NULL, no \code{estimate} argument is
  #'        supplied and \code{ergm()} uses its own default.
  #' @param eval.loglik logical or NULL.
  #'        If non-NULL, its value is forwarded to \code{ergm()} as the
  #'        \code{eval.loglik=} argument. If NULL, no \code{eval.loglik}
  #'        argument is supplied and \code{ergm()} uses its own default
  #'        (possibly depending on \code{estimate} and \code{control}).
  #' @param control list or `control.ergm` object, or NULL.
  #'        If NULL, `ergm()` receives no `control` argument and uses its defaults.
  #'        When a control object with an `init` vector of incompatible length
  #'        is passed, the wrapper drops `init` to let {ergm} recompute its
  #'        defaults.
  #' @param timeout numeric seconds; if set, evaluation runs under
  #'        `R.utils::withTimeout`.
  #' @param nodes optional data.frame for actor attributes and labels.
  #'        Used only when the LHS is a partition vector.
  #' @param dyads optional list of n×n matrices to attach to `%n%`.
  #'        Used only when the LHS is a partition vector.
  #' @return If `eval_call=TRUE`, an {ergm} fit. Otherwise, the unevaluated call.
  #' @examples
  #' \dontrun{
  #'   # Translate only
  #'   partition <- c(1, 1, 2, 2, 3)
  #'   call_only <- erpm(partition ~ groups(2), eval_call = FALSE)
  #'   print(call_only)
  #'
  #'   # Fit immediately with explicit estimate
  #'   fit <- erpm(partition ~ groups + cliques(3), estimate = "MLE")
  #' }
  #'
  #' @note
  #' The wrapper always imposes the \code{b1part} constraint and assumes a bipartite
  #' representation of the partition. The internal translation layer can later
  #' be extended by populating \code{effect_rename_map}, \code{wrap_with_proj1},
  #' and \code{wrap_with_B}.
  #'
  #' @note
  #' The behavior of \code{erpm()} is exercised in self-tests that compare
  #' ERPM-based fits to direct {ergm} calls on constructed bipartite networks.
  #' @export
  erpm <- function(formula, 
                   eval_call   = TRUE, 
                   verbose     = TRUE,
                   estimate    = NULL,
                   eval.loglik = NULL,
                   control     = NULL,
                   timeout     = NULL,
                   nodes       = NULL,
                   dyads       = list()) {

    # TRUE si l'argument verbose n'a PAS été fourni explicitement dans l'appel
    verbose_arg_missing <- missing(verbose)

    # --- 0) Normalize options ---------------------------------------------------
    if (!is.null(estimate)) {
      # Keep the set small and map legacy "MCMLE" to "MLE" so that
      # users can keep older scripts without diverging from {ergm}.
      estimate <- match.arg(estimate, c("MLE", "CD", "MPLE", "MCMLE"))
      if (identical(estimate, "MCMLE")) estimate <- "MLE"
    }
    # Do not fabricate a default eval.loglik here:
    # - if eval.loglik is NULL, we omit the argument and let ergm() decide;
    # - if non-NULL, we pass it through unchanged.

    # Guard: require a formula with a valid LHS
    if (!inherits(formula, "formula"))
      stop("Expected a `lhs ~ ...` formula with lhs = partition OR bipartite network.")

    # Compact string version of the original formula for logs and error reports.
    user_formula_str <- .tight(.oneline(formula))

    # --- 1) Inspect LHS type ----------------------------------------------------
    env0     <- environment(formula) %||% parent.frame()
    lhs      <- formula[[2]]
    rhs_expr <- formula[[3]]

    # Evaluate LHS safely in the formula env
    lhs_val <- tryCatch(eval(lhs, envir = env0), error = function(e) e)

    # --- 2) Case A: LHS is a partition -> build the bipartite network ----------
    if (!(inherits(lhs_val, "error")) &&
        is.atomic(lhs_val) && !inherits(lhs_val, "network")) {

      built <- build_bipartite_from_inputs(partition = lhs_val, nodes = nodes, dyads = dyads)
      nw2   <- built$network

      # Rebind the RHS over `nw` in a new formula environment to avoid global leaks.
      new_formula <- as.formula(bquote(nw ~ .(rhs_expr)))
      environment(new_formula) <- list2env(list(nw = nw2), parent = env0)

    } else {
      # --- 3) Case B: LHS is already a network -> do NOT rebuild ----------------
      if (!(inherits(lhs_val, "error")) && inherits(lhs_val, "network")) {
        bip <- tryCatch(lhs_val %n% "bipartite", error = function(e) NULL)
        if (is.null(bip) || is.na(bip))
          stop("LHS network is not bipartite or missing `%n% 'bipartite'` attribute.")
      }
      new_formula <- formula
      if (!(inherits(lhs_val, "error"))) environment(new_formula) <- env0
    }

    # --- 4) Translate RHS terms -------------------------------------------------
    orig_formula <- new_formula
    rhs_expr     <- new_formula[[3]]

    # Rename ERPM terms to {ergm} terms
    effect_rename_map <- c()
    wrap_with_proj1   <- c()
    wrap_with_B       <- c()

    # Translate and validate shapes immediately to fail fast on malformed inputs.
    safe_translate <- function() {
      rhs_terms <- .split_sum_terms(rhs_expr)
      translated_terms <- lapply(
        rhs_terms,
        .translate_one_term,
        rename_map = effect_rename_map,
        wrap_proj1 = wrap_with_proj1,
        wrap_B     = wrap_with_B
      )

      # Validate `b2degrange` intervals after translation.
      validate_groups_shape <- function(term_call) {
        if (is.symbol(term_call)) term_call <- as.call(list(term_call))
        if (!is.call(term_call)) return(term_call)
        f <- term_call[[1L]]
        if (is.symbol(f) && identical(f, as.name("b2degrange"))) {
          al <- as.pairlist(as.list(term_call)[-1L])
          env_eval <- environment(new_formula)
          from <- eval(al$from, envir = env_eval)
          to   <- if (is.null(al$to)) Inf else eval(al$to, envir = env_eval)
          if (is.finite(from) && from < 0L)
            stop(sprintf("groups(from|k): 'from'/'k' must be >= 0. Got: %d", from))
          if (!(is.infinite(to) || (is.finite(to) && to > from))) {
            to_str <- if (is.infinite(to)) "Inf" else as.character(to)
            stop(sprintf("groups(from,to): requires 'from' < 'to'. Got: from=%d, to=%s",
                        from, to_str))
          }
        }
        term_call
      }
      
      lapply(translated_terms, validate_groups_shape)
    }

    translated_terms <- tryCatch(
      safe_translate(),
      error = function(e) {
        if (isTRUE(verbose)) {
          cat(sprintf("[ERPM] call initial : erpm(%s)\n", user_formula_str))
          cat("\t error during translation: ", conditionMessage(e), "\n", sep = "")
        }
        stop(e)
      }
    )

    # Note: no artificial upper bound related to the number of groups.
    validate_groups_shape <- function(term_call) {
      if (is.symbol(term_call)) term_call <- as.call(list(term_call))
      if (!is.call(term_call)) return(term_call)
      f <- term_call[[1L]]
      if (is.symbol(f) && identical(f, as.name("b2degrange"))) {
        al <- as.pairlist(as.list(term_call)[-1L])
        env_eval <- environment(new_formula)
        from <- eval(al$from, envir = env_eval)
        to   <- if (is.null(al$to)) Inf else eval(al$to, envir = env_eval)
        if (is.finite(from) && from < 0L) stop(sprintf("groups(from|k): 'from'/'k' must be >= 0. Got: %d", from))
        if (!(is.infinite(to) || (is.finite(to) && to > from)))
          stop("groups(from,to): requires 'from' < 'to'.")
      }
      term_call
    }
    translated_terms <- lapply(translated_terms, validate_groups_shape)

    # Recompose the RHS from translated pieces.
    new_rhs <- if (length(translated_terms) == 1L) translated_terms[[1L]]
               else Reduce(function(x, y) call("+", x, y), translated_terms)

    # Inject translated RHS back into the formula.
    env <- environment(new_formula)
    new_formula[[3]] <- new_rhs
    if (!is.null(env)) environment(new_formula) <- env

    # --- 5) Constraints and control --------------------------------------------
    constraint_expression <- as.formula(~ b1part)

    ctrl_sym  <- NULL
    eval_env  <- if (!is.null(env)) env else parent.frame()

    # Build control only if user supplied a non-NULL control argument
    if (!is.null(control)) {
      ctrl <- if (inherits(control, "control.ergm")) control
              else do.call(ergm::control.ergm, as.list(control))

      # Guard: drop `init` when its length does not match the number of stats.
      # This keeps the API tolerant while letting {ergm} handle the usual cases.
      k <- length(summary(new_formula, constraints = ~ b1part))
      if (!is.null(ctrl$init) && length(ctrl$init) != k) ctrl$init <- NULL

      ctrl_sym  <- as.name(".ctrl_erpm")
      assign(as.character(ctrl_sym), ctrl, envir = eval_env)
    }

    # --- 6) Build the ergm(...) call -------------------------------------------
    call_args <- list(
      as.name("ergm"),
      new_formula,
      constraints = constraint_expression
    )
    if (!is.null(estimate))    call_args$estimate    <- estimate
    if (!is.null(eval.loglik)) call_args$eval.loglik <- eval.loglik
    if (!is.null(ctrl_sym))    call_args$control     <- ctrl_sym
    # Propagation rétroactive : si verbose a été explicitement fourni à erpm(),
    # on le transmet aussi à ergm(verbose = <valeur>).
    if (!verbose_arg_missing)  call_args$verbose     <- verbose

    ergm_call <- as.call(call_args)

    # --- 7) Compact logging -----------------------------------------------------
    if (isTRUE(verbose)) {
      final_fun <- as.call(list(as.name("ergm"), new_formula))
      final_str <- .tight(.oneline(final_fun))
      cat(sprintf("[ERPM] call initial : erpm(%s) -> call final : %s\n",
                  user_formula_str, final_str))
      cat("\t Constraints:  ~ b1part\n")
      cat("\t Options: estimate=",
          if (is.null(estimate)) "NULL" else estimate,
          ", eval.loglik=",
          if (is.null(eval.loglik)) "NULL" else eval.loglik,
          ", control=",
          if (is.null(ctrl_sym)) "NULL" else as.character(ctrl_sym),
          "\n", sep = "")
    }

    # --- 8) Evaluate or return the call ----------------------------------------
    if (!isTRUE(eval_call)) {
      if (isTRUE(verbose)) {
        cat("\t dry-run ergm call : ",
            paste(deparse(ergm_call, width.cutoff = 500L), collapse = " "),
            "\n", sep = "")
      }
      return(ergm_call)
    }

    # Evaluate with optional timeout, trap errors to reformat as ERPM errors.
    res <- try({
      if (is.null(timeout)) {
        eval(ergm_call, envir = eval_env)
      } else {
        R.utils::withTimeout(
          eval(ergm_call, envir = eval_env),
          timeout   = as.numeric(timeout),
          onTimeout = "silent"
        )
      }
    }, silent = TRUE)

    if (inherits(res, "try-error")) {
      msg <- paste0(
        "[ERPM ERROR]\n",
        "  user call : erpm(", user_formula_str, ")\n",
        "  ergm call : ", paste(deparse(ergm_call, width.cutoff = 500L), collapse = " "), "\n",
        "  message   : ", conditionMessage(attr(res, "condition"))
      )
      stop(msg, call. = FALSE)
    }
    return(res)
  }

  # ============================================================================
  # 7. Export functions to the global environment
  # ============================================================================

  assign("build_bipartite_from_inputs", build_bipartite_from_inputs,  envir = .GlobalEnv)
  assign(".normalize_groups_args",      .normalize_groups_args,       envir = .GlobalEnv)
  assign(".split_sum_terms",            .split_sum_terms,             envir = .GlobalEnv)
  assign(".translate_one_term",         .translate_one_term,          envir = .GlobalEnv)
  assign("erpm",                        erpm,                         envir = .GlobalEnv)
  assign(".__erpm_wrapper_loaded",      TRUE,                         envir = .GlobalEnv)
}