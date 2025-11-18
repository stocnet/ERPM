#' ERPM wrapper: translate ERPM formulas to {ergm} calls and optionally fit
#'
#' @file erpm_wrapper.R
#' @description
#' This module provides a thin wrapper to:
#' 1) Build a bipartite network from a partition, with optional node and dyadic inputs.
#' 2) Translate ERPM RHS terms (e.g., `groups`, `cov_match`, `cliques`) into {ergm} terms,
#'    with optional encapsulations (`Proj1`, `B`).
#' 3) Compose a standard `ergm(nw ~ <translated RHS>, constraints=~b1part, ...)` call.
#' 4) Either return the call (dry-run) or evaluate it and return the fitted model.
#'
#' The file also includes argument normalizers and internal validation helpers.
#' All user-facing console messages are in English for consistency with {ergm}.
#'
#' @keywords ERPM ERGM wrapper bipartite translation
#' @md

if (!exists(".__erpm_wrapper_loaded", envir = .GlobalEnv)) {

  # ============================================================================
  # Internal helpers
  # ============================================================================

  # Inline coalescing operator: return `a` if not NULL, else `b`.
  # Used to avoid nested if/else for option defaults.
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # One-line deparser for stable logging and error messages.
  .oneline <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")

  # Compacts whitespace for logs: "a   b\nc" -> "a b c".
  .tight   <- function(s) gsub("\\s+", "", s)

  #' Minimal bipartite test
  #'
  #' @param x any R object
  #' @return TRUE if `x` is a network with a finite `bipartite` attribute
  #' @keywords internal
  is_bipartite_network <- function(x) {
    inherits(x, "network") && !is.null(x %n% "bipartite") && is.finite(x %n% "bipartite")
  }

  # Try to locate a column to be used as node labels.
  # Preference order: `prefer`, then common aliases.
  .get_label_col <- function(nodes, prefer = "label") {
    stopifnot(is.data.frame(nodes))
    nms <- trimws(names(nodes))
    aliases <- c(prefer, "label", "nom", "name", "id")
    hit <- intersect(aliases, nms) # Avoid duplicate "label" entries when prefer = "label"
    if (length(hit)) hit[1L] else NULL
  }

  # Validate a nodes data.frame:
  # - must be a data.frame
  # - when a label-like column exists, it cannot contain duplicates
  .check_nodes_df <- function(nodes) {
    stopifnot(is.data.frame(nodes))
    lab <- .get_label_col(nodes)           # may be NULL if no obvious label column
    if (!is.null(lab) && anyDuplicated(nodes[[lab]]))
      stop(sprintf("nodes$%s contains duplicates.", lab))
    invisible(TRUE)
  }

  #' Validate dyadic n×n matrices
  #'
  #' Ensures matrices are square, of size n, and optionally that row/colnames
  #' match the actor label order. Stored on the network as `%n%` attributes.
  #'
  #' @param dyads named list of matrices
  #' @param n integer number of actors (mode 1)
  #' @param labels character vector of actor labels, used for name checks
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

  #' Normalize `groups(...)` into `b2degrange(from, to)`
  #'
  #' Accepted ERPM signatures:
  #'   groups
  #'   groups(k) / groups(size = k)
  #'   groups(from = ..., to = ...)
  #' Rules:
  #'   - from: integer >= 0
  #'   - to: integer > from  OR Inf
  #' Returns list(from = <int>, to = <int|quote(Inf)>)
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
  #' Turns `edges + triangles + nodematch("grp")` into a list of calls.
  #'
  #' @param expr a language object representing the RHS expression
  #' @return list of call objects
  #' @keywords internal
  .split_sum_terms <- function(expr) {
    if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
      return(c(.split_sum_terms(expr[[2]]),
               .split_sum_terms(expr[[3]])))
    }
    list(expr)
  }

  # Helper interne : stop message + usage bloc
  .stop_build_bipartite <- function(msg) {
    usage <- paste(
      "[ERPM] build_bipartite_from_inputs usage:",
      "",
      "  build_bipartite_from_inputs(partition, nodes = NULL, dyads = list())",
      "",
      "  # Exemple minimal :",
      "  partition <- c(1, 1, 2, 2, 2, 3)",
      "  nodes <- data.frame(",
      "    label  = c('A','B','C','D','E','F'),",
      "    gender = c(1,1,2,1,2,2),",
      "    age    = c(20,22,25,30,30,31)",
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

  #' Build a padded bipartite: G = n groups to maximize ERGM variability.
  #' Requires helper functions present in the wrapper:
  #'   .check_nodes_df(), .get_label_col(), .check_dyads()  
  #' 
  #' Actors are mode 1 (size n), groups are mode 2 (size G).
  #' Adjacency matrix is built symmetrically (undirected bipartite).
  #'
  #' @param partition integer vector of group ids, length n, values in 1..G
  #' @param nodes optional data.frame with a label column (auto-detected)
  #' @param dyads optional named list of n×n matrices to attach as `%n%` attributes
  #' @return list(network=..., partition=..., actor_labels=..., group_labels=...)
  #' @keywords internal 
  build_bipartite_from_inputs <- function(partition = NULL,
                                          nodes     = NULL,
                                          dyads     = list()) {
    # -- Checks on partition -------------------------------------------------------
    if (is.null(partition) || !is.atomic(partition) || length(partition) < 1L) {
      .stop_build_bipartite("partition must be a non-empty atomic vector.")
    }
    if (any(!is.finite(partition)) || any(partition < 1)) {
      .stop_build_bipartite("partition must contain finite positive integers.")
    }
    n <- length(partition)

    # -- Nodes / labels ------------------------------------------------------------
    if (is.null(nodes)) {
      labels <- sprintf("A%d", seq_len(n))
      nodes  <- data.frame(label = labels, stringsAsFactors = FALSE)
    } else {
      # erreurs de .check_nodes_df enveloppées
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

    # -- Padded groups: G = n ------------------------------------------------------
    G_obs <- max(partition, na.rm = TRUE)
    if (G_obs > n) {
      .stop_build_bipartite("max(partition) cannot exceed n when padding with G = n.")
    }
    G <- n  # force as many groups as actors

    # contrôle des matrices dyadiques
    ok_dyads <- try(.check_dyads(dyads, n, labels), silent = TRUE)
    if (inherits(ok_dyads, "try-error")) {
      .stop_build_bipartite(
        paste0("Invalid `dyads` list: ", conditionMessage(attr(ok_dyads, "condition")))
      )
    }

    # -- Names and adjacency -------------------------------------------------------
    g_names <- sprintf("G%d", seq_len(G))
    all_v   <- c(labels, g_names)

    adj <- matrix(0L, n + G, n + G, dimnames = list(all_v, all_v))
    # actor indices: 1..n ; group indices: n + (1..G)
    idx_actor <- seq_len(n)
    idx_group <- n + as.integer(partition)

    # Add bipartite edges for observed membership only; extra groups remain empty.
    adj[cbind(idx_actor, idx_group)] <- 1L
    adj[cbind(idx_group, idx_actor)] <- 1L

    # -- Build network object ------------------------------------------------------
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

  #' Translate a single ERPM term into an {ergm} term
  #'
  #' Handles:
  #' - `groups(...)` → `b2degrange(from, to)`
  #' - `cliques(...)` passthrough with normalized args
  #' - generic renames via `rename_map` and optional `Proj1`/`B` wrapping
  #'
  #' @param term_call a call object representing a single RHS term
  #' @param rename_map named character vector mapping old → new term names
  #' @param wrap_proj1 character vector of term names to wrap in `Proj1(~ term)`
  #' @param wrap_B character vector of term names to wrap in `B(~ term, form="nonzero")`
  #' @return a translated call suitable for {ergm}
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
    if (is.symbol(fun_sym) && identical(fun_sym, as.name("cliques"))) {
      k <- 2L; norm <- FALSE
      if (length(args_list) == 1L && (is.null(names(args_list)) || names(args_list)[1L] == "")) {
        k <- as.integer(args_list[[1L]])
      } else {
        al <- as.pairlist(args_list)
        if (!is.null(al$k))           k    <- as.integer(al$k)
        if (!is.null(al$clique_size)) k    <- as.integer(al$clique_size)
        if (!is.null(al$normalized))  norm <- isTRUE(al$normalized)
      }
      return(as.call(list(as.name("cliques"),
                          clique_size = k,
                          normalized  = norm)))
    }

    # Generic rename, then optional wrappers.
    fname <- if (is.symbol(fun_sym)) as.character(fun_sym) else deparse(fun_sym)[1L]
    fname <- if (fname %in% names(rename_map)) rename_map[[fname]] else fname
    out <- as.call(c(as.name(fname), args_list))
    if (fname %in% wrap_B)     out <- call("B",     call("~", out), form = "nonzero")
    if (fname %in% wrap_proj1) out <- call("Proj1", call("~", out))
    out
  }

  #' ERPM main wrapper: translate and optionally fit with {ergm}
  #'
  #' @param formula `lhs ~ <ERPM terms>`, where `lhs` is either a partition
  #'        vector or an already built bipartite `network` object.
  #' @param eval_call logical. If TRUE evaluate the {ergm} call, else return it.
  #' @param verbose logical. Print translation log and options.
  #' @param estimate one of `c("MLE","CD","MPLE","MCMLE")` or NULL.
  #'        If NULL, `ergm()` receives no `estimate` argument and uses its defaults.
  #' @param eval.loglik logical or NULL. If NULL and `estimate` is non-NULL,
  #'        default is chosen from `estimate`: CD→FALSE, MLE→TRUE, MPLE→FALSE.
  #'        If `estimate` is NULL and `eval.loglik` is NULL, `ergm()` decides.
  #' @param control list or `control.ergm` object, or NULL.
  #'        If NULL, `ergm()` receives no `control` argument and uses its defaults.
  #' @param timeout numeric seconds; if set, evaluation runs under `R.utils::withTimeout`.
  #' @param nodes optional data.frame for actor attributes and labels.
  #' @param dyads optional list of n×n matrices to attach to `%n%`.
  #' @return If `eval_call=TRUE`, an {ergm} fit. Otherwise, the unevaluated call.
  #' @examples
  #' \dontrun{
  #'   # Translate only
  #'   erpm(partition ~ groups, eval_call = FALSE)
  #'
  #'   # Fit immediately with explicit estimate
  #'   erpm(partition ~ groups, estimate = "MLE")
  #' }
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

    # --- 0) Normalize options ---------------------------------------------------
    if (!is.null(estimate)) {
      estimate <- match.arg(estimate, c("MLE","CD","MPLE","MCMLE"))
      if (identical(estimate, "MCMLE")) estimate <- "MLE"  # normalize to {ergm}'s estimate flag
    }

    # Default eval.loglik depends on estimate, if provided
    if (is.null(eval.loglik) && !is.null(estimate)) {
      eval.loglik <- switch(estimate,
                            "CD"   = FALSE,
                            "MLE"  = TRUE,
                            "MPLE" = FALSE,
                            TRUE)
    }

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
    wrap_with_proj1 <- c()
    wrap_with_B     <- c()

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

      # Optional: preserve previous behaviour when estimate="MLE" and control supplied.
      if (!is.null(estimate) && identical(estimate, "MLE") && is.null(ctrl$init)) {
        ctrl$init <- 0
      }

      # Guard: drop `init` when its length does not match the number of stats.
      k <- length(summary(new_formula, constraints = ~ b1part))
      if (!is.null(ctrl$init) && length(ctrl$init) != k) ctrl$init <- NULL

      ctrl_sym  <- as.name(sprintf(".ctrl_%08x", as.integer(runif(1, 0, .Machine$integer.max))))
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
    if (!is.null(control))     call_args$control     <- ctrl_sym

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
          if (is.null(control)) "NULL" else as.character(ctrl_sym),
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
  # Export functions to the global environment
  # ============================================================================
  assign("build_bipartite_from_inputs", build_bipartite_from_inputs,  envir = .GlobalEnv)
  assign(".normalize_groups_args",      .normalize_groups_args,       envir = .GlobalEnv)
  assign(".split_sum_terms",            .split_sum_terms,             envir = .GlobalEnv)
  assign(".translate_one_term",         .translate_one_term,          envir = .GlobalEnv)
  assign("erpm",                        erpm,                         envir = .GlobalEnv)
  assign(".__erpm_wrapper_loaded",      TRUE,                         envir = .GlobalEnv)
}

# #' ERPM wrapper: translate ERPM formulas to {ergm} calls and optionally fit
# #'
# #' @file erpm_wrapper.R
# #' @description
# #' This module provides a thin wrapper to:
# #' 1) Build a bipartite network from a partition, with optional node and dyadic inputs.
# #' 2) Translate ERPM RHS terms (e.g., `groups`, `cov_match`, `cliques`) into {ergm} terms,
# #'    with optional encapsulations (`Proj1`, `B`).
# #' 3) Compose a standard `ergm(nw ~ <translated RHS>, constraints=~b1part, ...)` call.
# #' 4) Either return the call (dry-run) or evaluate it and return the fitted model.
# #'
# #' The file also includes argument normalizers and internal validation helpers.
# #' All user-facing console messages are in English for consistency with {ergm}.
# #'
# #' @keywords ERPM ERGM wrapper bipartite translation
# #' @md

# if (!exists(".__erpm_wrapper_loaded", envir = .GlobalEnv)) {

#   # ============================================================================
#   # Internal helpers
#   # ============================================================================

#   # Inline coalescing operator: return `a` if not NULL, else `b`.
#   # Used to avoid nested if/else for option defaults.
#   `%||%` <- function(a, b) if (!is.null(a)) a else b

#   # One-line deparser for stable logging and error messages.
#   .oneline <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")

#   # Compacts whitespace for logs: "a   b\nc" -> "a b c".
#   .tight   <- function(s) gsub("\\s+", "", s)

#   #' Minimal bipartite test
#   #'
#   #' @param x any R object
#   #' @return TRUE if `x` is a network with a finite `bipartite` attribute
#   #' @keywords internal
#   is_bipartite_network <- function(x) {
#     inherits(x, "network") && !is.null(x %n% "bipartite") && is.finite(x %n% "bipartite")
#   }

#   # Try to locate a column to be used as node labels.
#   # Preference order: `prefer`, then common aliases.
#   .get_label_col <- function(nodes, prefer = "label") {
#     stopifnot(is.data.frame(nodes))
#     nms <- trimws(names(nodes))
#     aliases <- c(prefer, "label", "nom", "name", "id")
#     hit <- intersect(aliases, nms) # Avoid duplicate "label" entries when prefer = "label"
#     if (length(hit)) hit[1L] else NULL
#   }

#   # Validate a nodes data.frame:
#   # - must be a data.frame
#   # - when a label-like column exists, it cannot contain duplicates
#   .check_nodes_df <- function(nodes) {
#     stopifnot(is.data.frame(nodes))
#     lab <- .get_label_col(nodes)           # may be NULL if no obvious label column
#     if (!is.null(lab) && anyDuplicated(nodes[[lab]]))
#       stop(sprintf("nodes$%s contains duplicates.", lab))
#     invisible(TRUE)
#   }

#   #' Validate dyadic n×n matrices
#   #'
#   #' Ensures matrices are square, of size n, and optionally that row/colnames
#   #' match the actor label order. Stored on the network as `%n%` attributes.
#   #'
#   #' @param dyads named list of matrices
#   #' @param n integer number of actors (mode 1)
#   #' @param labels character vector of actor labels, used for name checks
#   #' @keywords internal
#   .check_dyads <- function(dyads, n, labels) {
#     if (length(dyads) == 0L) return(invisible(TRUE))
#     stopifnot(is.list(dyads))
#     for (nm in names(dyads)) {
#       M <- dyads[[nm]]
#       stopifnot(is.matrix(M), nrow(M) == n, ncol(M) == n)
#       # If names are present, enforce exact matching to actor label order.
#       if (!is.null(rownames(M)) && !is.null(colnames(M))) {
#         if (!identical(rownames(M), labels) || !identical(colnames(M), labels))
#           stop(sprintf("dyads['%s']: row/colnames must match the actor order.", nm))
#       }
#     }
#     invisible(TRUE)
#   }

#   #' Normalize `groups(...)` into `b2degrange(from, to)`
#   #'
#   #' Accepted ERPM signatures:
#   #'   groups
#   #'   groups(k) / groups(size = k)
#   #'   groups(from = ..., to = ...)
#   #' Rules:
#   #'   - from: integer >= 0
#   #'   - to: integer > from  OR Inf
#   #' Returns list(from = <int>, to = <int|quote(Inf)>)
#   .normalize_groups_args <- function(args_list) {
#     nm <- names(args_list)

#     # -- helpers ---------------------------------------------------------------
#     .deparse1 <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
#     .eval_num_scalar <- function(x, env = parent.frame()) {
#       if (is.numeric(x) && length(x) == 1L) return(x)
#       vx <- try(eval(x, envir = env), silent = TRUE)
#       if (inherits(vx, "try-error")) return(NULL)
#       if (is.numeric(vx) && length(vx) == 1L) return(vx)
#       NULL
#     }
#     .as_from_int <- function(x, env = parent.frame()) {
#       if (is.symbol(x) && identical(x, as.name("Inf")))
#         stop("groups(from,to): 'from' cannot be Inf.")
#       v <- .eval_num_scalar(x, env)
#       if (is.null(v) || !is.finite(v))
#         stop(sprintf("groups(from): must be a finite integer >= 0. Got: %s", .deparse1(x)))
#       iv <- as.integer(round(v))
#       if (!isTRUE(all.equal(v, iv))) stop(sprintf("groups(from): integer required. Got: %s", format(v)))
#       if (iv < 0L) stop("groups(from): must be >= 0.")
#       iv
#     }
#     .as_to_val <- function(x, env = parent.frame()) {
#       # Accept both symbol Inf and numeric Inf.
#       if ((is.symbol(x) && identical(x, as.name("Inf"))) ||
#           (is.numeric(x) && length(x) == 1L && is.infinite(x))) {
#         return(quote(Inf))
#       }
#       v <- .eval_num_scalar(x, env)
#       if (is.null(v) || !is.finite(v))
#         stop(sprintf("groups(to): must be a finite integer or Inf. Got: %s", .deparse1(x)))
#       iv <- as.integer(round(v))
#       if (!isTRUE(all.equal(v, iv))) stop(sprintf("groups(to): integer or Inf required. Got: %s", format(v)))
#       iv
#     }
#     # -------------------------------------------------------------------------

#     # No args: groups ≡ [1, Inf)
#     if (length(args_list) == 0L) return(list(from = 1L, to = quote(Inf)))

#     # Alias support: groups(size = k)
#     if (!is.null(nm) && "size" %in% nm && !("from" %in% nm) && !("to" %in% nm)) {
#       args_list <- list(args_list[["size"]]); nm <- NULL
#     }

#     # Single positional: groups(k) ≡ [k, k+1)
#     if (length(args_list) == 1L && (is.null(nm) || isTRUE(nm[1L] == ""))) {
#       k <- .as_from_int(args_list[[1L]], parent.frame())
#       return(list(from = k, to = k + 1L))
#     }

#     # Named pair: groups(from=..., to=...)
#     if (!is.null(nm) && all(c("from","to") %in% nm)) {
#       from <- .as_from_int(args_list[["from"]], parent.frame())
#       to   <- .as_to_val  (args_list[["to"]],   parent.frame())

#       # Resolve numeric to for check; keep quote(Inf) if Inf.
#       to_num <- if (is.language(to)) Inf else to
#       if (!(is.infinite(to_num) || (is.finite(to_num) && to_num > from))) {
#         to_str <- if (is.language(to)) "Inf" else as.character(to_num)
#         stop(sprintf("groups(from,to): requires 'from' < 'to'. Got: from=%d, to=%s", from, to_str))
#       }
#       return(list(from = from, to = to))
#     }

#     stop("groups(): use `groups`, `groups(k)`/`groups(size=k)`, or `groups(from=..,to=..)`. ")
#   }

#   #' Split a sum-of-terms RHS recursively
#   #'
#   #' Turns `edges + triangles + nodematch("grp")` into a list of calls.
#   #'
#   #' @param expr a language object representing the RHS expression
#   #' @return list of call objects
#   #' @keywords internal
#   .split_sum_terms <- function(expr) {
#     if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
#       return(c(.split_sum_terms(expr[[2]]),
#                .split_sum_terms(expr[[3]])))
#     }
#     list(expr)
#   }

#   # Helper interne : stop message + usage bloc
#   .stop_build_bipartite <- function(msg) {
#     usage <- paste(
#       "[ERPM] build_bipartite_from_inputs usage:",
#       "",
#       "  build_bipartite_from_inputs(partition, nodes = NULL, dyads = list())",
#       "",
#       "  # Exemple minimal :",
#       "  partition <- c(1, 1, 2, 2, 2, 3)",
#       "  nodes <- data.frame(",
#       "    label  = c('A','B','C','D','E','F'),",
#       "    gender = c(1,1,2,1,2,2),",
#       "    age    = c(20,22,25,30,30,31)",
#       "  )",
#       "  friendship <- matrix(c(",
#       "    0,1,1,1,0,0,",
#       "    1,0,0,0,1,0,",
#       "    1,0,0,0,1,0,",
#       "    1,0,0,0,0,0,",
#       "    0,1,1,0,0,1,",
#       "    0,0,0,0,1,0",
#       "  ), 6, 6, byrow = TRUE)",
#       "  dyads <- list(friendship = friendship)",
#       "",
#       "  built <- build_bipartite_from_inputs(partition, nodes = nodes, dyads = dyads)",
#       "  nw    <- built$network",
#       sep = "\n"
#     )
#     stop(paste0(msg, "\n\n", usage), call. = FALSE)
#   }

#   #' Build a padded bipartite: G = n groups to maximize ERGM variability.
#   #' Requires helper functions present in the wrapper:
#   #'   .check_nodes_df(), .get_label_col(), .check_dyads()  
#   #' 
#   #' Actors are mode 1 (size n), groups are mode 2 (size G).
#   #' Adjacency matrix is built symmetrically (undirected bipartite).
#   #'
#   #' @param partition integer vector of group ids, length n, values in 1..G
#   #' @param nodes optional data.frame with a label column (auto-detected)
#   #' @param dyads optional named list of n×n matrices to attach as `%n%` attributes
#   #' @return list(network=..., partition=..., actor_labels=..., group_labels=...)
#   #' @keywords internal 
#   build_bipartite_from_inputs <- function(partition = NULL,
#                                           nodes     = NULL,
#                                           dyads     = list()) {
#     # -- Checks on partition -------------------------------------------------------
#     if (is.null(partition) || !is.atomic(partition) || length(partition) < 1L) {
#       .stop_build_bipartite("partition must be a non-empty atomic vector.")
#     }
#     if (any(!is.finite(partition)) || any(partition < 1)) {
#       .stop_build_bipartite("partition must contain finite positive integers.")
#     }
#     n <- length(partition)

#     # -- Nodes / labels ------------------------------------------------------------
#     if (is.null(nodes)) {
#       labels <- sprintf("A%d", seq_len(n))
#       nodes  <- data.frame(label = labels, stringsAsFactors = FALSE)
#     } else {
#       # erreurs de .check_nodes_df enveloppées
#       ok_nodes <- try(.check_nodes_df(nodes), silent = TRUE)
#       if (inherits(ok_nodes, "try-error")) {
#         .stop_build_bipartite(
#           paste0("Invalid `nodes` data.frame: ", conditionMessage(attr(ok_nodes, "condition")))
#         )
#       }
#       if (nrow(nodes) != n) {
#         .stop_build_bipartite("nrow(nodes) must equal length(partition).")
#       }
#       labcol <- .get_label_col(nodes)
#       if (is.null(labcol)) {
#         labels <- sprintf("A%d", seq_len(n))
#         nodes$label <- labels
#       } else {
#         labels <- as.character(nodes[[labcol]])
#         if (anyNA(labels) || any(!nzchar(labels)))
#           .stop_build_bipartite("empty/NA labels are not allowed in `nodes`.")
#         if (anyDuplicated(labels))
#           .stop_build_bipartite("duplicate labels are not allowed in `nodes`.")
#         if (!("label" %in% names(nodes))) nodes$label <- labels
#       }
#     }

#     # -- Padded groups: G = n ------------------------------------------------------
#     G_obs <- max(partition, na.rm = TRUE)
#     if (G_obs > n) {
#       .stop_build_bipartite("max(partition) cannot exceed n when padding with G = n.")
#     }
#     G <- n  # force as many groups as actors

#     # contrôle des matrices dyadiques
#     ok_dyads <- try(.check_dyads(dyads, n, labels), silent = TRUE)
#     if (inherits(ok_dyads, "try-error")) {
#       .stop_build_bipartite(
#         paste0("Invalid `dyads` list: ", conditionMessage(attr(ok_dyads, "condition")))
#       )
#     }

#     # -- Names and adjacency -------------------------------------------------------
#     g_names <- sprintf("G%d", seq_len(G))
#     all_v   <- c(labels, g_names)

#     adj <- matrix(0L, n + G, n + G, dimnames = list(all_v, all_v))
#     # actor indices: 1..n ; group indices: n + (1..G)
#     idx_actor <- seq_len(n)
#     idx_group <- n + as.integer(partition)

#     # Add bipartite edges for observed membership only; extra groups remain empty.
#     adj[cbind(idx_actor, idx_group)] <- 1L
#     adj[cbind(idx_group, idx_actor)] <- 1L

#     # -- Build network object ------------------------------------------------------
#     nw <- network::network(adj, directed = FALSE, matrix.type = "adjacency")
#     network::set.network.attribute(nw, "bipartite", n)
#     network::set.vertex.attribute(nw, "vertex.names", all_v)

#     # Push actor attributes; pad groups with NA
#     if (ncol(nodes) > 1L) {
#       for (a in setdiff(names(nodes), c("label"))) {
#         vals <- nodes[[a]]
#         network::set.vertex.attribute(nw, a, c(vals, rep(NA, G)))
#       }
#     }

#     # Attach dyadic n×n matrices as %n% attributes, with enforced dimnames order
#     if (length(dyads)) {
#       for (nm in names(dyads)) {
#         M <- dyads[[nm]]
#         dimnames(M) <- list(labels, labels)
#         nw %n% nm <- M
#       }
#     }

#     list(
#       network      = nw,
#       partition    = partition,
#       actor_labels = labels,
#       group_labels = g_names
#     )
#   }

#   #   nw <- network::network(adj, directed = FALSE, matrix.type = "adjacency")
#   #   network::set.network.attribute(nw, "bipartite", n)
#   #   network::set.vertex.attribute(nw, "vertex.names", all_v)

#   #   # Push actor attributes; pad groups with NA
#   #   if (ncol(nodes) > 1L) {
#   #     for (a in setdiff(names(nodes), c("label"))) {
#   #       vals <- nodes[[a]]
#   #       network::set.vertex.attribute(nw, a, c(vals, rep(NA, G)))
#   #     }
#   #   }

#   #   # Attach dyadic n×n matrices as %n% attributes, with enforced dimnames order
#   #   if (length(dyads)) {
#   #     for (nm in names(dyads)) {
#   #       M <- dyads[[nm]]
#   #       dimnames(M) <- list(labels, labels)
#   #       nw %n% nm <- M
#   #     }
#   #   }

#   #   list(network = nw,
#   #       partition = partition,
#   #       actor_labels = labels,
#   #       group_labels = g_names)
#   # }
  
#   #' Translate a single ERPM term into an {ergm} term
#   #'
#   #' Handles:
#   #' - `groups(...)` → `b2degrange(from, to)`
#   #' - `cliques(...)` passthrough with normalized args
#   #' - generic renames via `rename_map` and optional `Proj1`/`B` wrapping
#   #'
#   #' @param term_call a call object representing a single RHS term
#   #' @param rename_map named character vector mapping old → new term names
#   #' @param wrap_proj1 character vector of term names to wrap in `Proj1(~ term)`
#   #' @param wrap_B character vector of term names to wrap in `B(~ term, form="nonzero")`
#   #' @return a translated call suitable for {ergm}
#   #' @keywords internal
#   .translate_one_term <- function(term_call,
#                                   rename_map,
#                                   wrap_proj1 = character(),
#                                   wrap_B     = character()) {
#     # Turn bare symbols into zero-arg calls to unify processing.
#     if (is.symbol(term_call)) term_call <- as.call(list(term_call))
#     if (!is.call(term_call)) return(term_call)

#     fun_sym   <- term_call[[1L]]
#     args_list <- as.list(term_call)[-1L]

#     # Special case: `groups(...)`
#     if (is.symbol(fun_sym) && identical(fun_sym, as.name("groups"))) {
#       gt <- .normalize_groups_args(args_list)
#       return(as.call(list(as.name("b2degrange"), from = gt$from, to = gt$to)))
#     }

#     # Special case: `cliques(...)` (normalize aliases, keep as ERPM term name)
#     if (is.symbol(fun_sym) && identical(fun_sym, as.name("cliques"))) {
#       k <- 2L; norm <- FALSE
#       if (length(args_list) == 1L && (is.null(names(args_list)) || names(args_list)[1L] == "")) {
#         k <- as.integer(args_list[[1L]])
#       } else {
#         al <- as.pairlist(args_list)
#         if (!is.null(al$k))           k    <- as.integer(al$k)
#         if (!is.null(al$clique_size)) k    <- as.integer(al$clique_size)
#         if (!is.null(al$normalized))  norm <- isTRUE(al$normalized)
#       }
#       return(as.call(list(as.name("cliques"),
#                           clique_size = k,
#                           normalized  = norm)))
#     }

#     # Generic rename, then optional wrappers.
#     fname <- if (is.symbol(fun_sym)) as.character(fun_sym) else deparse(fun_sym)[1L]
#     fname <- if (fname %in% names(rename_map)) rename_map[[fname]] else fname
#     out <- as.call(c(as.name(fname), args_list))
#     if (fname %in% wrap_B)     out <- call("B",     call("~", out), form = "nonzero")
#     if (fname %in% wrap_proj1) out <- call("Proj1", call("~", out))
#     out
#   }

#   #' ERPM main wrapper: translate and optionally fit with {ergm}
#   #'
#   #' @param formula `lhs ~ <ERPM terms>`, where `lhs` is either a partition
#   #'        vector or an already built bipartite `network` object.
#   #' @param eval_call logical. If TRUE evaluate the {ergm} call, else return it.
#   #' @param verbose logical. Print translation log and options.
#   #' @param estimate one of `c("MLE","CD","MPLE","MCMLE")`. "MCMLE" is normalized to "MLE".
#   #' @param eval.loglik logical or NULL. If NULL, chosen from `estimate`:
#   #'        CD→FALSE, MLE→TRUE, MPLE→FALSE.
#   #' @param control list or `control.ergm` object. Passed to {ergm}.
#   #' @param timeout numeric seconds; if set, evaluation runs under `R.utils::withTimeout`.
#   #' @param nodes optional data.frame for actor attributes and labels.
#   #' @param dyads optional list of n×n matrices to attach to `%n%`.
#   #' @return If `eval_call=TRUE`, an {ergm} fit. Otherwise, the unevaluated call.
#   #' @examples
#   #' \dontrun{
#   #'   # Translate only
#   #'   erpm(partition ~ groups, eval_call = FALSE)
#   #'
#   #'   # Fit immediately
#   #'   erpm(partition ~ groups, estimate = "MLE")
#   #' }
#   #' @export
#   erpm <- function(formula, 
#                   eval_call   = TRUE, 
#                   verbose     = TRUE,
#                   estimate    = c("MLE","CD","MPLE","MCMLE"),
#                   eval.loglik = NULL,
#                   control     = NULL,
#                   timeout     = NULL,
#                   nodes       = NULL,
#                   dyads       = list()) {

#     # --- 0) Normalize options ---------------------------------------------------
#     estimate <- match.arg(estimate)
#     if (identical(estimate, "MCMLE")) estimate <- "MLE"  # unify into {ergm}'s estimate flag

#     # Default eval.loglik depends on estimate
#     if (is.null(eval.loglik)) {
#       eval.loglik <- switch(estimate,
#                             "CD"   = FALSE,
#                             "MLE"  = TRUE,
#                             "MPLE" = FALSE,
#                             TRUE)
#     }

#     # Guard: require a formula with a valid LHS
#     if (!inherits(formula, "formula"))
#       stop("Expected a `lhs ~ ...` formula with lhs = partition OR bipartite network.")

#     # Compact string version of the original formula for logs and error reports.
#     user_formula_str <- .tight(.oneline(formula))

#     # --- 1) Inspect LHS type ----------------------------------------------------
#     env0     <- environment(formula) %||% parent.frame()
#     lhs      <- formula[[2]]
#     rhs_expr <- formula[[3]]

#     # Evaluate LHS safely in the formula env
#     lhs_val <- tryCatch(eval(lhs, envir = env0), error = function(e) e)

#     # --- 2) Case A: LHS is a partition -> build the bipartite network ----------
#     if (!(inherits(lhs_val, "error")) &&
#         is.atomic(lhs_val) && !inherits(lhs_val, "network")) {

#       built <- build_bipartite_from_inputs(partition = lhs_val, nodes = nodes, dyads = dyads)
#       nw2   <- built$network

#       # Rebind the RHS over `nw` in a new formula environment to avoid global leaks.
#       new_formula <- as.formula(bquote(nw ~ .(rhs_expr)))
#       environment(new_formula) <- list2env(list(nw = nw2), parent = env0)

#     } else {
#       # --- 3) Case B: LHS is already a network -> do NOT rebuild ----------------
#       if (!(inherits(lhs_val, "error")) && inherits(lhs_val, "network")) {
#         bip <- tryCatch(lhs_val %n% "bipartite", error = function(e) NULL)
#         if (is.null(bip) || is.na(bip))
#           stop("LHS network is not bipartite or missing `%n% 'bipartite'` attribute.")
#       }
#       new_formula <- formula
#       if (!(inherits(lhs_val, "error"))) environment(new_formula) <- env0
#     }

#     # --- 4) Translate RHS terms -------------------------------------------------
#     orig_formula <- new_formula
#     rhs_expr     <- new_formula[[3]]

#     # Rename ERPM terms to {ergm} terms
#     effect_rename_map <- c()
#     # Terms that must run in the 1-mode projection or as a biased constraint
#     wrap_with_proj1 <- c()
#     wrap_with_B     <- c()

#     # Translate and validate shapes immediately to fail fast on malformed inputs.
#     safe_translate <- function() {
#       rhs_terms <- .split_sum_terms(rhs_expr)
#       translated_terms <- lapply(
#         rhs_terms,
#         .translate_one_term,
#         rename_map = effect_rename_map,
#         wrap_proj1 = wrap_with_proj1,
#         wrap_B     = wrap_with_B
#       )

#       # Validate `b2degrange` intervals after translation.
#       validate_groups_shape <- function(term_call) {
#         if (is.symbol(term_call)) term_call <- as.call(list(term_call))
#         if (!is.call(term_call)) return(term_call)
#         f <- term_call[[1L]]
#         if (is.symbol(f) && identical(f, as.name("b2degrange"))) {
#           al <- as.pairlist(as.list(term_call)[-1L])
#           env_eval <- environment(new_formula)
#           from <- eval(al$from, envir = env_eval)
#           to   <- if (is.null(al$to)) Inf else eval(al$to, envir = env_eval)
#           if (is.finite(from) && from < 0L)
#             stop(sprintf("groups(from|k): 'from'/'k' must be >= 0. Got: %d", from))
#           if (!(is.infinite(to) || (is.finite(to) && to > from))) {
#             to_str <- if (is.infinite(to)) "Inf" else as.character(to)
#             stop(sprintf("groups(from,to): requires 'from' < 'to'. Got: from=%d, to=%s",
#                         from, to_str))
#           }
#         }
#         term_call
#       }
      
#       lapply(translated_terms, validate_groups_shape)
#     }

#     translated_terms <- tryCatch(
#       safe_translate(),
#       error = function(e) {
#         if (isTRUE(verbose)) {
#           cat(sprintf("[ERPM] call initial : erpm(%s)\n", user_formula_str))
#           cat("\t error during translation: ", conditionMessage(e), "\n", sep = "")
#         }
#         stop(e)
#       }
#     )

#     # Note: no artificial upper bound related to the number of groups.
#     # If the requested range selects no groups, the statistic is 0 without error.
#     # Keep only consistency checks (from>=0, from<to) as above.
#     validate_groups_shape <- function(term_call) {
#       if (is.symbol(term_call)) term_call <- as.call(list(term_call))
#       if (!is.call(term_call)) return(term_call)
#       f <- term_call[[1L]]
#       if (is.symbol(f) && identical(f, as.name("b2degrange"))) {
#         al <- as.pairlist(as.list(term_call)[-1L])
#         env_eval <- environment(new_formula)
#         from <- eval(al$from, envir = env_eval)
#         to   <- if (is.null(al$to)) Inf else eval(al$to, envir = env_eval)
#         if (is.finite(from) && from < 0L) stop(sprintf("groups(from|k): 'from'/'k' must be >= 0. Got: %d", from))
#         if (!(is.infinite(to) || (is.finite(to) && to > from)))
#           stop("groups(from,to): requires 'from' < 'to'.")
#       }
#       term_call
#     }
#     translated_terms <- lapply(translated_terms, validate_groups_shape)

#     # Recompose the RHS from translated pieces.
#     new_rhs <- if (length(translated_terms) == 1L) translated_terms[[1L]]
#                else Reduce(function(x, y) call("+", x, y), translated_terms)

#     # Inject translated RHS back into the formula.
#     env <- environment(new_formula)
#     new_formula[[3]] <- new_rhs
#     if (!is.null(env)) environment(new_formula) <- env

#     # --- 5) Constraints and control --------------------------------------------
#     constraint_expression <- as.formula(~ b1part)

#     # Accept either a `control.ergm` object or a raw list of control args.
#     ctrl <- if (inherits(control, "control.ergm")) control
#             else if (is.null(control)) ergm::control.ergm()
#             else do.call(ergm::control.ergm, as.list(control))
            
#     if (identical(estimate, "MLE") && is.null(ctrl$init)) {
#       ctrl$init <- 0  # évite l’appel au MPLE
#     }

#     # Guard: drop `init` when its length does not match the number of stats.
#     k <- length(summary(new_formula, constraints = ~ b1part))
#     if (!is.null(ctrl$init) && length(ctrl$init) != k) ctrl$init <- NULL

#     # Bind control into the evaluation environment and pass it as a symbol.
#     ctrl_sym  <- as.name(sprintf(".ctrl_%08x", as.integer(runif(1, 0, .Machine$integer.max))))
#     eval_env  <- if (!is.null(env)) env else parent.frame()
#     assign(as.character(ctrl_sym), ctrl, envir = eval_env)

#     # --- 6) Build the ergm(...) call -------------------------------------------
#     ergm_call <- as.call(list(
#       as.name("ergm"),
#       new_formula,
#       constraints = constraint_expression,
#       estimate    = estimate,
#       eval.loglik = eval.loglik,
#       control     = ctrl_sym
#     ))

#     # --- 7) Compact logging -----------------------------------------------------
#     if (isTRUE(verbose)) {
#       final_fun <- as.call(list(as.name("ergm"), new_formula))
#       final_str <- .tight(.oneline(final_fun))
#       cat(sprintf("[ERPM] call initial : erpm(%s) -> call final : %s\n",
#                   user_formula_str, final_str))
#       cat("\t Constraints:  ~ b1part\n")
#       cat("\t Options: estimate=", estimate,
#           ", eval.loglik=", eval.loglik,
#           ", control=", as.character(ctrl_sym), "\n", sep = "")
#     }

#     # --- 8) Evaluate or return the call ----------------------------------------
#     if (!isTRUE(eval_call)) {
#       if (isTRUE(verbose)) {
#         cat("\t dry-run ergm call : ",
#             paste(deparse(ergm_call, width.cutoff = 500L), collapse = " "),
#             "\n", sep = "")
#       }
#       return(ergm_call)
#     }

#     # Evaluate with optional timeout, trap errors to reformat as ERPM errors.
#     res <- try({
#       if (is.null(timeout)) {
#         eval(ergm_call, envir = eval_env)
#       } else {
#         R.utils::withTimeout(
#           eval(ergm_call, envir = eval_env),
#           timeout   = as.numeric(timeout),
#           onTimeout = "silent"
#         )
#       }
#     }, silent = TRUE)

#     if (inherits(res, "try-error")) {
#       msg <- paste0(
#         "[ERPM ERROR]\n",
#         "  user call : erpm(", user_formula_str, ")\n",
#         "  ergm call : ", paste(deparse(ergm_call, width.cutoff = 500L), collapse = " "), "\n",
#         "  message   : ", conditionMessage(attr(res, "condition"))
#       )
#       stop(msg, call. = FALSE)
#     }
#     return(res)
#   }

#   # ============================================================================
#   # Export functions to the global environment
#   # ============================================================================
#   assign("build_bipartite_from_inputs", build_bipartite_from_inputs,  envir = .GlobalEnv)
#   # assign(".normalize_cliques_args",     .normalize_cliques_args,      envir = .GlobalEnv)
#   assign(".normalize_groups_args",      .normalize_groups_args,       envir = .GlobalEnv)
#   assign(".split_sum_terms",            .split_sum_terms,             envir = .GlobalEnv)
#   assign(".translate_one_term",         .translate_one_term,          envir = .GlobalEnv)
#   assign("erpm",                        erpm,                         envir = .GlobalEnv)
#   assign(".__erpm_wrapper_loaded",      TRUE,                         envir = .GlobalEnv)
# }