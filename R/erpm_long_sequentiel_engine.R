# ################################################################################
# # FILE: R/erpm_long_sequentiel_engine.R
# ################################################################################
# #' ERPM longitudinal engine (PLS / sequential) helpers (internal)
# #'
# #' @name erpm_long_sequentiel_engine
# #' @note erpm_long_sequentiel_engine.R
# #'
# #' @description
# #' This file contains the PLS (sequential) internal engine used by \code{erpm_long()}.
# #'
# #' At each time step, it:
# #' \enumerate{
# #'   \item extracts the partition and optional time-varying inputs (nodes, dyads, group labels);
# #'   \item builds a padded bipartite \pkg{network} from those inputs;
# #'   \item attaches past partitions (raw) as network attributes for inertial InitErgmTerm.* terms;
# #'   \item constructs a time-specific formula \code{nw ~ RHS} and calls \code{erpm()}.
# #' }
# #'
# #' IMPORTANT:
# #' Inertial computations are delegated to InitErgmTerm.* terms.
# #' This engine only detects the requested past depth (past_influence) from the RHS to:
# #'   - validate feasibility vs T;
# #'   - skip early times t <= d (no enough past partitions) in PLS.
# #'
# #' @keywords ERPM ERGM longitudinal engine
# NULL

# #' Run the PLS longitudinal loop over t = 1..T (internal)
# #'
# #' @noRd
# .erpm_long_sequentiel_run <- function(parts,
#                                      T,
#                                      env0,
#                                      rhs_expr,
#                                      eval.call,
#                                      verbose,
#                                      debug,
#                                      estimate,
#                                      eval.loglik,
#                                      control,
#                                      timeout,
#                                      seed,
#                                      nodes,
#                                      dyads,
#                                      group_labels) {

#   # ---------------------------------------------------------------------------
#   # Developer debug helper
#   # ---------------------------------------------------------------------------
#   debug_deep <- is.character(debug) && identical(debug, "deep")
#   debug_any  <- isTRUE(debug) || debug_deep

#   .dbg <- function(...) {
#     if (!isTRUE(debug_any)) return(invisible(NULL))
#     cat("[ERPM_LONG|PLS|DEBUG] ", sprintf(...), "\n", sep = "")
#     invisible(NULL)
#   }

#   # ---------------------------------------------------------------------------
#   # RHS analysis: detect inertial terms and extract requested past depth d
#   # ---------------------------------------------------------------------------
#   .detect_inertial_depth <- function(rhs, env_eval) {
#     # Heuristic/conservative detector:
#     # - collect calls whose function name starts with "inertia_"
#     # - extract past_influence (aliases: pi, d) if present
#     # - return max requested depth across inertial terms (0 if none)
#     if (is.null(rhs)) return(list(has_inertia = FALSE, depth = 0L, terms = character(0)))

#     is_inertia_call <- function(cc) {
#       if (!is.call(cc)) return(FALSE)
#       if (!is.symbol(cc[[1L]])) return(FALSE)
#       fn <- as.character(cc[[1L]])
#       startsWith(fn, "inertia_")
#     }

#     calls <- .erpm_long_collect_calls(rhs, is_inertia_call)
#     if (!length(calls)) return(list(has_inertia = FALSE, depth = 0L, terms = character(0)))

#     # Extract numeric scalar safely (allow literals and simple symbols).
#     .eval_scalar_int <- function(x) {
#       if (is.null(x)) return(NULL)
#       if (is.numeric(x) && length(x) == 1L && is.finite(x)) {
#         iv <- as.integer(round(x))
#         if (isTRUE(all.equal(x, iv))) return(iv)
#         return(NULL)
#       }
#       v <- try(eval(x, envir = env_eval), silent = TRUE)
#       if (inherits(v, "try-error")) return(NULL)
#       if (is.numeric(v) && length(v) == 1L && is.finite(v)) {
#         iv <- as.integer(round(v))
#         if (isTRUE(all.equal(v, iv))) return(iv)
#       }
#       NULL
#     }

#     depths <- integer(0)
#     terms  <- character(0)

#     for (cc in calls) {
#       fn <- as.character(cc[[1L]])
#       al <- try(as.pairlist(as.list(cc)[-1L]), silent = TRUE)
#       if (inherits(al, "try-error")) al <- NULL

#       pi_expr <- NULL
#       if (!is.null(al) && length(al)) {
#         if (!is.null(al$past_influence)) pi_expr <- al$past_influence
#         else if (!is.null(al$pi))        pi_expr <- al$pi
#         else if (!is.null(al$d))         pi_expr <- al$d
#       }

#       if (is.null(pi_expr)) {
#         stop(sprintf(
#           "[ERPM_LONG|PLS] inertial term %s() detected but missing `past_influence` (or alias pi/d).",
#           fn
#         ), call. = FALSE)
#       }

#       d <- .eval_scalar_int(pi_expr)
#       if (is.null(d) || d < 1L) {
#         stop(sprintf(
#           "[ERPM_LONG|PLS] invalid past_influence for %s(): must be an integer >= 1.",
#           fn
#         ), call. = FALSE)
#       }

#       depths <- c(depths, d)
#       terms  <- c(terms, fn)
#     }

#     list(
#       has_inertia = TRUE,
#       depth       = max(depths),
#       terms       = unique(terms)
#     )
#   }

#   rhs_info <- .detect_inertial_depth(rhs_expr, env_eval = env0)
#   has_inertia <- isTRUE(rhs_info$has_inertia)
#   d <- as.integer(rhs_info$depth)

#   .dbg("RHS inertia detected: %s", if (has_inertia) "YES" else "NO")
#   if (has_inertia) {
#     .dbg("Inertial terms: %s", paste(rhs_info$terms, collapse = ", "))
#     .dbg("Requested past depth d = %d", d)
#   }

#   # Feasibility checks vs T.
#   if (has_inertia) {
#     if (d >= T) {
#       stop(sprintf(
#         "[ERPM_LONG|PLS] past_influence d=%d is not compatible with T=%d partitions (no estimable time points).",
#         d, T
#       ), call. = FALSE)
#     }
#   }

#   # In PLS, we skip early times t <= d when inertia is required.
#   t0 <- if (has_inertia) (d + 1L) else 1L

#   # ---------------------------------------------------------------------------
#   # Storage
#   # ---------------------------------------------------------------------------
#   calls            <- vector("list", T)
#   fits             <- vector("list", T)
#   nets             <- vector("list", T)
#   history_timeline <- vector("list", T)

#   # ---------------------------------------------------------------------------
#   # Optional user-facing header complement
#   # ---------------------------------------------------------------------------
#   if (isTRUE(verbose) && has_inertia) {
#     cat(sprintf("[ERPM_LONG|PLS] inertial terms detected: %s\n", paste(rhs_info$terms, collapse = ", ")))
#     cat(sprintf("[ERPM_LONG|PLS] past_influence depth d = %d -> will fit times t=%d..%d (skip t<=%d)\n",
#                 d, t0, T, d))
#     cat("------------------------------------------------------------\n")
#   }

#   # ---------------------------------------------------------------------------
#   # Small helper: standardized longitudinal attachments for InitErgmTerm.*
#   # ---------------------------------------------------------------------------
#   .attach_history_attrs <- function(nw, mode, t, T, p_current, p_past) {

#     # Always attach a stable set of attributes (even when no inertia is requested).
#     # InitErgmTerm.* code can rely on these labels without special-casing.
#     network::set.network.attribute(nw, "erpm_mode", mode)
#     network::set.network.attribute(nw, "erpm_time_index", as.integer(t))
#     network::set.network.attribute(nw, "erpm_T",          as.integer(T))

#     # Current partition (raw, integer vector).
#     network::set.network.attribute(nw, "erpm_partition_current", as.integer(p_current))

#     # Past partitions (raw, list of integer vectors; lag 1..d, most recent first).
#     # Keep the list label stable even when d=0.
#     network::set.network.attribute(nw, "erpm_past_depth", as.integer(length(p_past)))
#     network::set.network.attribute(nw, "erpm_past_partitions", p_past)

#     # Convenience per-lag labels (optional but stable when present).
#     if (length(p_past)) {
#       for (lag in seq_along(p_past)) {
#         network::set.network.attribute(nw, paste0("erpm_partition_lag", lag), p_past[[lag]])
#       }
#     }

#     invisible(TRUE)
#   }

#   # ---------------------------------------------------------------------------
#   # Main loop
#   # ---------------------------------------------------------------------------
#   for (t in seq_len(T)) {

#     # Skip early time points if inertia requires past partitions.
#     if (has_inertia && t <= d) {
#       if (isTRUE(verbose)) {
#         cat(sprintf("\n[ERPM_LONG|PLS] --- t=%d ---\n", t))
#         cat(sprintf("[ERPM_LONG|PLS] skip: t=%d <= d=%d (not enough past partitions for inertial terms)\n", t, d))
#       }
#       .dbg("t=%d skipped (requires past partitions up to d=%d)", t, d)

#       nets[[t]] <- NULL
#       calls[[t]] <- NULL
#       fits[[t]] <- NULL
#       history_timeline[[t]] <- list(
#         t            = t,
#         skipped      = TRUE,
#         skip_reason  = sprintf("t<=d (%d<=%d)", t, d),
#         rhs          = .erpm_long_rhs_oneline(rhs_expr),
#         inertia      = TRUE,
#         past_depth   = d
#       )
#       next
#     }

#     # Partition at time t.
#     p_t <- as.integer(round(parts[[t]]))

#     # Per-time inputs (or shared object).
#     nodes_t <- .erpm_long_get_t(nodes, t, T, what = "nodes")
#     dyads_t <- .erpm_long_get_t(dyads, t, T, what = "dyads")
#     glab_t  <- .erpm_long_get_t(group_labels, t, T, what = "group_labels")

#     .dbg("t=%d: extracting inputs", t)
#     .dbg("t=%d: nodes_t=%s | dyads_t=%s | group_labels_t=%s",
#          t,
#          if (is.null(nodes_t)) "NULL" else "data.frame",
#          if (is.null(dyads_t)) "NULL" else if (is.list(dyads_t)) sprintf("list(%d)", length(dyads_t)) else class(dyads_t)[1L],
#          if (is.null(glab_t)) "NULL" else if (is.list(glab_t)) sprintf("list(%d)", length(glab_t)) else class(glab_t)[1L])

#     # Integrity checks (engine-level).
#     if (!is.null(nodes_t)) {
#       if (!is.data.frame(nodes_t)) {
#         stop(sprintf("[ERPM_LONG|PLS] nodes[[%d]] must be a data.frame.", t), call. = FALSE)
#       }
#       if (nrow(nodes_t) != length(p_t)) {
#         stop(sprintf("[ERPM_LONG|PLS] nrow(nodes[[%d]]) must equal length(partition[[%d]]).", t, t), call. = FALSE)
#       }
#     }

#     if (is.null(dyads_t)) dyads_t <- list()
#     if (!is.list(dyads_t)) {
#       stop(sprintf("[ERPM_LONG|PLS] internal error: dyads[[%d]] must be a list after normalization.", t), call. = FALSE)
#     }

#     if (isTRUE(verbose)) {
#       cat(sprintf("\n[ERPM_LONG|PLS] --- t=%d ---\n", t))
#       cat("[ERPM_LONG|PLS] partition: ", .erpm_long_partition_info(p_t), "\n", sep = "")
#       if (is.null(nodes_t)) {
#         cat("[ERPM_LONG|PLS] nodes: NULL (auto labels will be used)\n")
#       } else {
#         cat("[ERPM_LONG|PLS] nodes: cols={", paste(names(nodes_t), collapse = ","), "}\n", sep = "")
#       }

#       if (length(dyads_t) == 0L) {
#         cat("[ERPM_LONG|PLS] dyads: empty\n")
#       } else {
#         dims <- vapply(dyads_t, function(M) paste(dim(M), collapse = "x"), "")
#         cat("[ERPM_LONG|PLS] dyads: names={", paste(names(dyads_t), collapse = ","), "}\n", sep = "")
#         cat("[ERPM_LONG|PLS] dyads: dims ={", paste(sprintf("%s:%s", names(dims), dims), collapse = " "), "}\n", sep = "")
#       }
#     }

#     # Build current bipartite network.
#     .dbg("t=%d: building bipartite network", t)
#     built <- build_bipartite_from_inputs(
#       partition    = p_t,
#       nodes        = nodes_t,
#       dyads        = dyads_t,
#       group_labels = glab_t
#     )
#     nw_t <- built$network

#     # Always attach standardized current + past partitions for InitErgmTerm.* terms.
#     past_parts <- list()
#     if (has_inertia) {
#       # Attach list of past partitions with lags 1..d, most recent first.
#       past_parts <- vector("list", d)
#       for (lag in seq_len(d)) {
#         past_parts[[lag]] <- as.integer(round(parts[[t - lag]]))
#       }
#       .dbg("t=%d: prepared past partitions: lags=1..%d", t, d)
#     }

#     .attach_history_attrs(
#       nw        = nw_t,
#       mode      = "sequentiel",
#       t         = t,
#       T         = T,
#       p_current = p_t,
#       p_past    = past_parts
#     )

#     if (isTRUE(debug_any) && !debug_deep) {
#       .dbg("t=%d: network summary (shallow debug)", t)
#       print(summary(nw_t))
#     }

#     # RHS is shared in PLS; inertial terms read past partitions from nw attributes.
#     rhs_t <- rhs_expr

#     if (isTRUE(verbose)) {
#       cat("[ERPM_LONG|PLS] RHS(t): ", .erpm_long_rhs_oneline(rhs_t), "\n", sep = "")
#       cat("[ERPM_LONG|PLS] calling erpm()...\n")
#     }
#     .dbg("t=%d: building formula nw ~ RHS", t)

#     eval_env_t <- list2env(list(nw = nw_t), parent = env0)
#     f_t <- as.formula(bquote(nw ~ .(rhs_t)))
#     environment(f_t) <- eval_env_t

#     .dbg("t=%d: calling erpm(eval.call=%s)", t, if (isTRUE(eval.call)) "TRUE" else "FALSE")
#     res_t <- erpm(
#       formula      = f_t,
#       eval.call    = eval.call,
#       verbose      = verbose,
#       estimate     = estimate,
#       eval.loglik  = eval.loglik,
#       control      = control,
#       timeout      = timeout,
#       seed         = seed,
#       nodes        = NULL,
#       dyads        = list(),
#       group_labels = NULL
#     )

#     nets[[t]] <- nw_t

#     if (isTRUE(has_inertia)) {
#       history_timeline[[t]] <- list(
#         t            = t,
#         skipped      = FALSE,
#         rhs          = .erpm_long_rhs_oneline(rhs_t),
#         has_nodes    = !is.null(nodes_t),
#         dyads_names  = if (length(dyads_t)) names(dyads_t) else character(0),
#         inertia      = TRUE,
#         past_depth   = d,
#         active_terms = rhs_info$terms
#       )
#     } else {
#       history_timeline[[t]] <- list(NULL)
#     }

#     if (isTRUE(eval.call)) {
#       fits[[t]] <- res_t
#       calls[[t]] <- try(getCall(res_t), silent = TRUE)
#       if (inherits(calls[[t]], "try-error") || is.null(calls[[t]])) calls[[t]] <- NULL
#       .dbg("t=%d: fit returned (class=%s)", t, paste(class(res_t), collapse = "/"))
#     } else {
#       calls[[t]] <- res_t
#       fits[[t]]  <- NULL
#       if (isTRUE(verbose)) {
#         cat("[ERPM_LONG|PLS] dry-run call:\n")
#         cat("  ", paste(deparse(res_t, width.cutoff = 500L), collapse = " "), "\n", sep = "")
#       }
#       .dbg("t=%d: dry-run call returned", t)
#     }
#   }

#   if (isTRUE(verbose)) {
#     cat("\n------------------------------------------------------------\n")
#     cat("[ERPM_LONG|PLS] Done\n")
#     if (has_inertia) {
#       cat(sprintf("[ERPM_LONG|PLS] partitions: %d | fitted times: %d..%d | skipped: %d\n",
#                   T, t0, T, d))
#     } else {
#       cat(sprintf("[ERPM_LONG|PLS] partitions: %d | fitted times: 1..%d\n", T, T))
#     }
#     cat(sprintf("[ERPM_LONG|PLS] dry-run: %s\n", if (isTRUE(eval.call)) "NO" else "YES"))
#     cat("============================================================\n")
#   }

#   structure(
#     list(
#       mode             = "sequentiel",
#       calls            = calls,
#       fits             = fits,
#       networks         = nets,
#       history_timeline = history_timeline,
#       rhs_inertia      = list(
#         has_inertia = has_inertia,
#         past_depth  = d,
#         terms       = if (has_inertia) rhs_info$terms else character(0)
#       )
#     ),
#     class = "erpm_long"
#   )
# }