# =====================================================================
# FILE: R/erpm_long_engine.R
# =====================================================================
#' ERPM longitudinal engine helpers (internal)
#'
#' @name erpm_long_engine
#' @note erpm_long_engine.R
#' 
#' @description
#' This file contains the core internal engine used by \code{erpm_long()} to fit
#' (or dry-run) a sequence of ERPM models across time \eqn{t = 1,\dots,T}.
#'
#' At each time step, it:
#' \enumerate{
#'   \item extracts the partition and optional time-varying inputs (nodes, dyads, group labels);
#'   \item builds a padded bipartite \pkg{network} from those inputs;
#'   \item activates inertial terms only when enough past is available (e.g., no lagged effects at \eqn{t=1});
#'   \item attaches per-lag “inertia” attributes onto the current network;
#'   \item builds a time-specific formula \code{nw ~ RHS(t)} and calls \code{erpm()}.
#' }
#'
#' The resulting object stores the per-time call(s), fit(s), built network(s),
#' and a timeline describing which inertial terms were active and what attributes
#' were attached at each \eqn{t}.
#'
#' @keywords ERPM ERGM longitudinal engine

# ============================================================================
# Helpers
# ============================================================================

#' Select inertial calls that are active at time t.
#'
#' For each inertial term call, this checks its registered specification
#' (in \code{.erpm_long_inertia_registry}) and its required history depth.
#'
#' A term becomes active only once \code{t} is strictly greater than the number
#' of past steps it needs, so the engine never requests unavailable lags.
#'
#' @param inertial_calls List of inertial term calls detected on the shared RHS.
#' @param t Integer time index (1-based).
#' @param env0 Evaluation environment of the original formula.
#' @param debug Logical; if TRUE, emit optional debug traces.
#' @param .dbg Optional debug function called with compact messages.
#' @return A list of inertial term calls that are active at time \code{t}.
#' @noRd
.erpm_long_active_inertial_calls_for_t <- function(inertial_calls, t, env0, debug = FALSE, .dbg = NULL) {
  # No inertial terms were provided : nothing to do .
  if (!length(inertial_calls)) return(list())

  out <- list()
  for (cc in inertial_calls) {
    # Map the call to a standard term name and retrieve its inertia spec.
    nm <- .erpm_long_term_name(cc)
    spec <- .erpm_long_inertia_registry[[nm]]

    # Unknown term name => silently skip (allows partial registries).
    if (is.null(spec)) next

    # "past_influence" returns how many past steps are needed (depth d).
    d <- spec$past_influence(cc, env0)

    # Activate only when all required past networks exist: nets[[t-lag]] for lag=1..d.
    if (t > d) out[[length(out) + 1L]] <- cc
  }

  # Optional debug trace
  if (isTRUE(debug) && is.function(.dbg)) {
    .dbg("active_inertial_calls_for_t t=", as.integer(t),
         " -> ", if (length(out)) paste(vapply(out, .erpm_long_term_name, ""), collapse = ", ") else "(none)")
  }

  out
}

#' Attach inertia attributes for active inertial terms at time t.
#'
#' For each active inertial term, this reads its required history depth and
#' attaches one attribute per lag \eqn{\ell = 1,\dots,d} onto the current network.
#'
#' Attribute names follow a stable convention:
#' \code{erpm_inertia__<term>__lag<k>}.
#'
#' These attributes are later consumed by translated ERGM terms (e.g., via nodal
#' or dyadic covariates stored as network attributes).
#'
#' @param nw_t The current bipartite membership network at time \code{t}.
#' @param nets List of previously built networks (length \code{t-1} or more).
#' @param parts List of partitions (one per time point), used when available.
#' @param t Integer time index (1-based).
#' @param active_calls List of inertial term calls active at time \code{t}.
#' @param env0 Evaluation environment of the original formula.
#' @param debug Logical; if TRUE, enable additional diagnostics within builders.
#' @return A list with components:
#'   \itemize{
#'     \item \code{network}: the updated network with attached inertia attributes;
#'     \item \code{attached}: an index describing what was attached (term, lag, name).
#'   }
#' @noRd
.erpm_long_attach_inertia_for_t <- function(nw_t, nets, parts, t, active_calls, env0, debug = FALSE) {

  # No active inertial terms => return the network unchanged and no index.
  if (!length(active_calls)) return(list(network = nw_t, attached = NULL))

  attached <- list()

  for (cc in active_calls) {

    # Resolve inertia specification for this term.
    nm <- .erpm_long_term_name(cc)
    spec <- .erpm_long_inertia_registry[[nm]]
    if (is.null(spec)) next

    # Required history depth for this term.
    d <- spec$past_influence(cc, env0)

    # Optional static ERGM terms used by fallback summary path.
    st <- NULL
    if (is.function(spec$static_terms)) {
      st <- spec$static_terms(cc, env0)

      # Normalize to list so downstream code is uniform.
      if (!is.list(st)) st <- list(st)
    }

    # Attach one attribute per lag: lag=1 refers to t-1, lag=2 refers to t-2, etc.
    for (lag in seq_len(d)) {
      j <- t - lag

      # Do not access time <= 0 (no past networks available).
      if (j < 1L) break

      # Apply standard naming scheme so translated terms can predict names (should match initergmterm parser).
      attr_name <- paste0("erpm_inertia__", nm, "__lag", lag)

      if (is.function(spec$build_attr)) { # build_attr is a function defined into inertial registry

        # Preferred path: spec knows how to construct the attribute value.
        val <- spec$build_attr(nw_t, nets, parts, t, lag, cc, env0, debug = debug)

        # Store it on the current network so ERGM terms can read it.
        nw_t <- network::set.network.attribute(nw_t, attr_name, val)

        # If the attribute is a named vector/list, keep the names for reporting.
        stat_names <- names(val)
        if (is.null(stat_names)) stat_names <- character(0)

        # Track exactly what was attached for later inspection and verbose logs.
        attached[[length(attached) + 1L]] <- list(
          term       = nm,
          lag        = as.integer(lag),
          attr_name  = attr_name,
          stat_names = stat_names
        )
      } else { # NOT TESTED : no function is defined into inertial registry but a changestat

        # Fallback path: compute a static summary on the previous network and attach it.
        nw_prev <- nets[[j]]
        vec <- .erpm_long_compute_static_summary_one(nw_prev, st, env0)
        nw_t <- network::set.network.attribute(nw_t, attr_name, vec)

        attached[[length(attached) + 1L]] <- list(
          term       = nm,
          lag        = as.integer(lag),
          attr_name  = attr_name,
          stat_names = names(vec)
        )
      }
    }
  }

  # Attach a single index attribute so users can discover what was attached without guessing names.
  nw_t <- network::set.network.attribute(nw_t, "erpm_inertia_attached_index", attached)
  list(network = nw_t, attached = attached)
}

#' Run the longitudinal ERPM fitting loop over \eqn{t = 1,\dots,T}.
#'
#' Internal engine used by \code{erpm_long()} to iterate across time points and
#' fit (or dry-run) a sequence of ERPM models on time-indexed partitions.
#'
#' For each time point \code{t}, the routine:
#' \enumerate{
#'   \item extracts the current partition and time-specific inputs (\code{nodes},
#'         \code{dyads}, \code{group_labels});
#'   \item builds the padded bipartite membership network via
#'         \code{build_bipartite_from_inputs()};
#'   \item determines which inertial terms are active at \code{t} (based on their
#'         \code{past_influence} and the current time index);
#'   \item attaches per-lag inertia attributes onto the current network when needed;
#'   \item constructs the time-specific RHS as static terms plus active inertial terms;
#'   \item calls \code{erpm()} in fit mode (\code{eval.call=TRUE}) or dry-run mode
#'         (\code{eval.call=FALSE}).
#' }
#'
#' The function stores, for each \code{t}, the network actually used for fitting,
#' the resulting fit (if any), the produced call, and a timeline describing inertial
#' activation and attribute attachments.
#'
#' @param parts List of partitions (one per time point), length \code{T}.
#'   Each partition is expected to be integer-like; values are rounded and coerced
#'   to integer defensively.
#' @param T Integer; number of time points.
#' @param env0 Evaluation environment of the original formula; used as parent for
#'   per-time evaluation environments.
#' @param rhs_expr Shared RHS expression (unevaluated), used for deep-debug translation.
#' @param static_terms List of static RHS term calls (always active).
#' @param inertial_calls List of inertial RHS term calls (activation depends on \code{t}).
#' @param det Inertia detection descriptor produced upstream (e.g. enabled flag, names).
#' @param eval.call Logical; if TRUE, fit models, else return unevaluated calls.
#' @param verbose Logical; if TRUE, print compact per-time traces.
#' @param debug Logical or character; if TRUE, prints additional diagnostics.
#'   If \code{debug="deep"}, performs deep checks including RHS translation and
#'   nodecov/nodefactor readiness checks on the built network.
#' @param .dbg Optional debug function used by subroutines when \code{debug=TRUE}.
#' @param estimate Character or NULL, forwarded to \code{erpm()}.
#' @param eval.loglik Logical or NULL, forwarded to \code{erpm()}.
#' @param control Control object or list, forwarded to \code{erpm()}.
#' @param timeout Numeric seconds or NULL, forwarded to \code{erpm()}.
#' @param seed Integer or NULL, forwarded to \code{erpm()}.
#' @param nodes Optional node tables, either time-invariant or time-indexed.
#'   If provided at time \code{t}, must be a data.frame with \code{nrow(nodes[[t]])}
#'   equal to \code{length(parts[[t]])}.
#' @param dyads Optional dyadic inputs, either time-invariant or time-indexed.
#'   Normalized upstream to a list (possibly empty) for each time point.
#' @param group_labels Optional group label vectors, either time-invariant or time-indexed.
#'
#' @return An object of class \code{"erpm_long"} with components:
#'   \itemize{
#'     \item \code{calls}: list of per-time \code{ergm()} calls (or extracted calls from fits);
#'     \item \code{fits}: list of per-time fitted models (or NULLs in dry-run mode);
#'     \item \code{networks}: list of per-time networks used for fitting;
#'     \item \code{history_timeline}: list describing inertial activation and inertia
#'           attribute attachments per time.
#'   }
#'
#' @noRd
.erpm_long_run <- function(parts,
                          T,
                          env0,
                          rhs_expr,
                          static_terms,
                          inertial_calls,
                          det,
                          eval.call,
                          verbose,
                          debug,
                          .dbg,
                          estimate,
                          eval.loglik,
                          control,
                          timeout,
                          seed,
                          nodes,
                          dyads,
                          group_labels) {

  # Pre-allocate result containers for performance and predictable structure.
  calls             <- vector("list", T)
  fits              <- vector("list", T)
  nets              <- vector("list", T)
  history_timeline  <- vector("list", T)

  debug_deep <- is.character(debug) && identical(debug, "deep")
  debug_any  <- isTRUE(debug) || debug_deep

  for (t in seq_len(T)) {
    # Partition is expected to be integer-like; round defensively to avoid floating artifacts.
    p_t <- as.integer(round(parts[[t]]))

    # Allow nodes/dyads/group_labels to be either constant or time-indexed.
    nodes_t <- .erpm_long_get_t(nodes, t, T, what = "nodes")
    dyads_t <- .erpm_long_get_t(dyads, t, T, what = "dyads")
    glab_t  <- .erpm_long_get_t(group_labels, t, T, what = "group_labels")

    # Validate nodes_t only when provided (NULL means "auto labels").
    if (!is.null(nodes_t)) {
      if (!is.data.frame(nodes_t)) {
        stop(sprintf("[ERPM_LONG] nodes[[%d]] must be a data.frame.", t), call. = FALSE)
      }
      # Node table must align with the partition length (one row per node).
      if (nrow(nodes_t) != length(p_t)) {
        stop(sprintf("[ERPM_LONG] nrow(nodes[[%d]]) must equal length(partition[[%d]]).", t, t), call. = FALSE)
      }
    }

    # dyads_t is normalized upstream to be a list (possibly empty).
    if (is.null(dyads_t)) dyads_t <- list()
    if (!is.list(dyads_t)) {
      stop(sprintf("[ERPM_LONG] internal error: dyads[[%d]] must be a list after normalization.", t), call. = FALSE)
    }

    # Verbose logging for user-facing tracing of inputs.
    if (isTRUE(verbose)) {
      cat(sprintf("\n[ERPM_LONG] --- t=%d ---\n", t))
      cat("[ERPM_LONG] partition: ", .erpm_long_partition_info(p_t), "\n", sep = "")
      if (is.null(nodes_t)) {
        cat("[ERPM_LONG] nodes: NULL (auto labels will be used)\n")
      } else {
        cat("[ERPM_LONG] nodes: cols={", paste(names(nodes_t), collapse = ","), "}\n", sep = "")
      }

      if (length(dyads_t) == 0L) {
        cat("[ERPM_LONG] dyads: empty\n")
      } else {
        # Print dyad matrix dimensions to catch mismatches early.
        dims <- vapply(dyads_t, function(M) paste(dim(M), collapse = "x"), "")
        cat("[ERPM_LONG] dyads: names={", paste(names(dyads_t), collapse = ","), "}\n", sep = "")
        cat("[ERPM_LONG] dyads: dims ={", paste(sprintf("%s:%s", names(dims), dims), collapse = " "), "}\n", sep = "")
      }
    }

    # Build the bipartite network representation expected by the wrapper.
    built <- build_bipartite_from_inputs(
      partition    = p_t,
      nodes        = nodes_t,
      dyads        = dyads_t,
      group_labels = glab_t
    )
    nw_t <- built$network

    # Debug mode prints a network summary (useful when diagnosing attribute attachment).
    if (isTRUE(debug_any) && !debug_deep) print(summary(nw_t))

    # Determine which inertial terms are active at this time step.
    active_calls <- .erpm_long_active_inertial_calls_for_t(inertial_calls, t, env0, debug = debug_any, .dbg = .dbg)
    attached_index <- NULL

    if (length(active_calls)) {
      # Attach per-lag inertia attributes onto nw_t using past networks and partitions.
      ret <- .erpm_long_attach_inertia_for_t(nw_t, nets, parts, t, active_calls, env0, debug = debug_any)
      nw_t <- ret$network
      attached_index <- ret$attached

      # Record what happened at time t for later inspection and reproducibility.
      history_timeline[[t]] <- list(
        t              = t,
        active_terms   = vapply(active_calls, .erpm_long_term_name, character(1)),
        attached_index = attached_index
      )

      # Human-readable description of attached attributes.
      if (isTRUE(verbose)) {
        cat("[ERPM_LONG] inertia: attached per-lag attributes:\n")
        for (it in attached_index) {
          cat("  - ", it$attr_name, " (", it$term, ", lag=", it$lag, "): ",
              paste(it$stat_names, collapse = ", "), "\n", sep = "")
        }
      }
    } else {
      # Keep timeline structure aligned even when nothing is attached.
      history_timeline[t] <- list(NULL)
      # Optional informational message when inertial terms exist but cannot apply at t=1.
      if (isTRUE(verbose) && isTRUE(det$enabled) && t == 1L) {
        cat("[ERPM_LONG] inertial terms detected but t=1 => no inertial effects applied.\n")
      }
    }

    # Store the final network used at time t (after inertia attachment if any).
    nets[[t]] <- nw_t

    # RHS is base terms plus inertial terms that are currently active.
    rhs_terms_t <- c(static_terms, active_calls)
    rhs_t <- .erpm_long_combine_terms(rhs_terms_t)

    if (isTRUE(verbose)) {
      cat("[ERPM_LONG] RHS(t): ", .erpm_long_rhs_oneline(rhs_t), "\n", sep = "")
      if (length(inertial_calls)) {
        # Report which inertial terms are present but still inactive at this t.
        inactive <- setdiff(
          vapply(inertial_calls, .erpm_long_term_name, character(1)),
          vapply(active_calls,   .erpm_long_term_name, character(1))
        )
        if (length(inactive)) {
          cat("[ERPM_LONG] inertial terms not yet active at this t: ",
              paste(unique(inactive), collapse = ", "), "\n", sep = "")
        }
      }
      cat("[ERPM_LONG] calling erpm()...\n")
    }

    # Evaluate the formula in an environment that exposes 'nw' as the built network.
    eval_env_t <- list2env(list(nw = nw_t), parent = env0)

    # -----------------------------------------------------------------------
    # Deep debug: show nodecov()/nodefactor() readiness by mode.
    # We inspect the TRANSLATED RHS so it matches what ergm() will see.
    # -----------------------------------------------------------------------
    if (isTRUE(debug_deep)) {
      cat(sprintf("[ERPM_LONG][deep] --- t=%d deep data checks ---\n", t))
      print(summary(nw_t))

      # Translate RHS(t) using the same translation pipeline as erpm(),
      # so deep debug follows the effective ergm terms.
      rhs_tr <- try(
        .erpm_translate_rhs_expr(
          rhs_expr  = rhs_t,
          env_eval  = eval_env_t,
          effect_rename_map = c(),
          wrap_with_proj1   = c(),
          wrap_with_B       = c()
        ),
        silent = TRUE
      )

      if (inherits(rhs_tr, "try-error")) {
        cat("[ERPM_LONG][deep] RHS translation failed while preparing deep debug.\n")
        cat("[ERPM_LONG][deep] message: ", conditionMessage(attr(rhs_tr, "condition")), "\n", sep = "")
      } else {
        cat("[ERPM_LONG][deep] RHS(t) translated: ", .erpm_long_rhs_oneline(rhs_tr), "\n", sep = "")
        .erpm_long_deep_debug_nodecov_inputs(nw_t, rhs_tr)
      }

      cat(sprintf("[ERPM_LONG][deep] --- t=%d end deep data checks ---\n", t))
    }

    f_t <- as.formula(bquote(nw ~ .(rhs_t)))
    environment(f_t) <- eval_env_t

    # Delegate to the main wrapper. nodes/dyads/group_labels are already baked into nw_t.
    res_t <- erpm(
      formula      = f_t,
      eval.call    = eval.call,
      verbose      = verbose,
      estimate     = estimate,
      eval.loglik  = eval.loglik,
      control      = control,
      timeout      = timeout,
      seed         = seed,
      nodes        = NULL,
      dyads        = list(),
      group_labels = NULL
    )

    if (isTRUE(eval.call)) {
      # Fit object path: store fit and attempt to extract the underlying ergm call.
      fits[[t]] <- res_t
      calls[[t]] <- try(getCall(res_t), silent = TRUE)
      if (inherits(calls[[t]], "try-error") || is.null(calls[[t]])) calls[[t]] <- NULL

      # Verbose reporting of coefficients when available.
      if (isTRUE(verbose)) {
        cf <- try(stats::coef(res_t), silent = TRUE)
        if (!inherits(cf, "try-error") && length(cf)) {
          cat("[ERPM_LONG] fitted coefficients:\n")
          for (nm in names(cf)) {
            cat(sprintf("  - %s = %s\n", nm, format(cf[[nm]], digits = 6)))
          }
        } else {
          cat("[ERPM_LONG] fitted coefficients: unavailable\n")
        }
      }
    } else {
      # Dry-run path: res_t is the unevaluated call produced by erpm().
      calls[[t]] <- res_t
      fits[t] <- list(NULL)

      if (isTRUE(verbose)) {
        cat("[ERPM_LONG] dry-run call:\n")
        # Deparse to keep the call readable in logs even if it spans multiple lines.
        cat("  ", paste(deparse(res_t, width.cutoff = 500L), collapse = " "), "\n", sep = "")
      }
    }
  }

  if (isTRUE(verbose)) {
    cat("\n------------------------------------------------------------\n")
    cat("[ERPM_LONG] Done\n")
    cat(sprintf("[ERPM_LONG] models: %d | dry-run: %s\n", T, if (isTRUE(eval.call)) "NO" else "YES"))
    cat("============================================================\n")
  }

  # Return a structured object so downstream code can inspect networks, calls, and timeline.
  structure(
    list(
      calls    = calls,
      fits     = fits,
      networks = nets,
      history_timeline  = history_timeline
    ),
    class = "erpm_long"
  )
}