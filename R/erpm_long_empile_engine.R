################################################################################
# FILE: R/erpm_long_empile_engine.R
################################################################################
#' ERPM PLE engine: build stacked meta-network and attach PLE attributes
#'
#' @name erpm_long_empile_engine
#' @note erpm_long_empile_engine.R
#'
#' @description
#' This file implements the PLE ("empile") engine used by \code{erpm_long()} to
#' turn a timeline of partitions into a single stacked bipartite meta-network.
#' The engine is deliberately structured as three modules:
#' \enumerate{
#'   \item \strong{Standard meta-network build}: concatenate selected partitions into a
#'         single partition, optionally bind node tables, and block-diagonalize dyadic
#'         matrices, then call \code{build_bipartite_from_inputs()}.
#'   \item \strong{Timeblock construction}: attach a \code{timeblock} vertex attribute so
#'         each actor/group vertex can be mapped back to its originating time index.
#'   \item \strong{Inertial plumbing}: when inertial terms are present, attach the
#'         additional network attributes expected by inertial InitErgmTerms
#'         (currently tailored to \code{inertia_groups} in PLE mode).
#' }
#'
#' The implementation is intentionally explicit about sizes, offsets, and invariants
#' (actor ranges per block, disjoint group-id ranges, dyads dimensional checks) to make
#' debugging PLE construction failures straightforward.
#'
#' @keywords ERPM ERGM longitudinal PLE empile engine
################################################################################

# ==============================================================================
# Low-level helpers (existing behavior)
# ==============================================================================

#' Enforce constant actor block size across estimation blocks (internal helper)
#' @param partitions Timeline of partitions (list).
#' @param selected_partition_indices Integer vector of time indices included in the meta-network.
#' @return Integer scalar: the common actor count per selected block.
#' @noRd
.erpm_ple_assert_constant_nbr_actors_per_selected_partition <- function(partitions, selected_partition_indices) {
  nbr_actors_by_t <- vapply(partitions, length, integer(1))
  nbr_actors_per_selected_partition <- nbr_actors_by_t[selected_partition_indices[1L]]
  if (any(nbr_actors_by_t[selected_partition_indices] != nbr_actors_per_selected_partition)) {
    stop(sprintf(
      "[ERPM_PLE] PLE(inertia_groups) requires constant nbr_actors_per_selected_partition across selected_partition_indices. Got nbr_actors_by_t[selected_partition_indices]=%s.",
      paste(nbr_actors_by_t[selected_partition_indices], collapse = ",")
    ))
  }
  nbr_actors_per_selected_partition
}

#' Build the per-block past-partitions container expected by inertia_groups (internal helper)
#' @param partitions Timeline of partitions (list).
#' @param selected_partition_indices Integer vector of current-time indices used for estimation.
#' @param d Past influence depth (number of lags).
#' @param nbr_actors_per_selected_partition Actor count per selected block.
#' @return List of length B; each element is a list of length d with past partitions.
#' @noRd
.erpm_ple_make_erpm_block_past_partitions <- function(partitions, selected_partition_indices, d, nbr_actors_per_selected_partition) {
  B <- length(selected_partition_indices)
  out <- vector("list", B)

  for (b in seq_len(B)) {
    t_cur <- selected_partition_indices[b]
    pb <- vector("list", d)

    for (lag in seq_len(d)) {
      t_past <- t_cur - lag
      if (t_past < 1L) {
        stop(sprintf("[ERPM_PLE] internal: t_past=%d < 1 (b=%d, lag=%d). selected_partition_indicesmust start at d+1.",
                     t_past, b, lag))
      }

      p <- partitions[[t_past]]

      if (is.null(p) || !is.atomic(p) || length(p) != nbr_actors_per_selected_partition) {
        stop(sprintf("[ERPM_PLE] past partition invalid at time=%d (expected atomic length nbr_actors_per_selected_partition=%d).",
                     t_past, nbr_actors_per_selected_partition))
      }
      if (anyNA(p)) {
        stop(sprintf("[ERPM_PLE] past partition at time=%d contains NA.", t_past))
      }

      pb[[lag]] <- p
    }

    out[[b]] <- pb
  }

  out
}

#' Count observed groups in a partition (bookkeeping helper)
#' @param partition Atomic partition vector.
#' @return Integer scalar: number of observed groups (by first occurrence order).
#' @noRd
.erpm_ple_G_obs <- function(partition) {
  f <- factor(partition, levels = unique(partition))
  nlevels(f)
}

#' Build concatenated + shifted meta partition (internal helper)
#' @param partitions Timeline of partitions (list).
#' @param selected_partition_indices Integer vector of time indices included in the meta-network.
#' @return List with meta_partition and offset bookkeeping vectors.
#' @noRd
.erpm_ple_make_meta_partition <- function(partitions, selected_partition_indices) {

  nbr_actors_by_t <- vapply(partitions, length, integer(1)) # vapply enforces the output type
  nbr_actors_by_selected_partitions <- nbr_actors_by_t[selected_partition_indices]
  actor_start_index_by_block <- c(0L, cumsum(nbr_actors_by_selected_partitions))[seq_along(nbr_actors_by_selected_partitions)]  # length B
  actor_start_index_by_t <- rep(NA_integer_, length(nbr_actors_by_t)) # defined only for selected times
  actor_start_index_by_t[selected_partition_indices] <- actor_start_index_by_block

  meta <- integer(0)

  for (b in seq_along(selected_partition_indices)) {
    t <- selected_partition_indices[b]
    p <- partitions[[t]]

    # Stable relabeling within block: groups become 1..G_obs
    f <- factor(p, levels = unique(p))
    g <- as.integer(f)

    # Shift group ids so each block uses disjoint group-id ranges
    meta <- c(meta, g + actor_start_index_by_block[b])
  }


  list(
    meta_partition = meta,
    nbr_actors_by_t        = nbr_actors_by_t,
    selected_partition_indices       = selected_partition_indices,
    nbr_actors_by_selected_partitions         = nbr_actors_by_selected_partitions,
    actor_start_index_by_block    = actor_start_index_by_block,
    actor_start_index_by_t   = actor_start_index_by_t
  )
}

#' Build meta nodes data.frame for build_bipartite_from_inputs (internal helper)
#' @param nodes NULL or list of per-time node data.frames.
#' @param partitions Timeline of partitions (list).
#' @param selected_partition_indices Integer vector of time indices included in the meta-network.
#' @return NULL or a single bound data.frame (one row per meta-actor).
#' @noRd
.erpm_ple_make_meta_nodes_df <- function(nodes, partitions, selected_partition_indices) {
  if (is.null(nodes)) return(NULL)

  T <- length(partitions)
  if (!is.list(nodes) || length(nodes) != T) {
    stop("[ERPM_PLE] nodes must be NULL or a list of length T (one data.frame per partition).")
  }

  cols_ref <- NULL
  out <- NULL

  for (t in selected_partition_indices) {
    df <- nodes[[t]]
    if (!is.data.frame(df)) stop(sprintf("[ERPM_PLE] nodes[[%d]] must be a data.frame.", t))
    if (nrow(df) != length(partitions[[t]])) {
      stop(sprintf("[ERPM_PLE] nodes[[%d]] has %d rows but partition has %d actors.",
                   t, nrow(df), length(partitions[[t]])))
    }
    if (!("label" %in% colnames(df))) {
      stop(sprintf("[ERPM_PLE] nodes[[%d]] must contain a 'label' column.", t))
    }

    # Enforce schema consistency (same columns, same order) across selected times
    if (is.null(cols_ref)) {
      cols_ref <- colnames(df)
    } else if (!identical(colnames(df), cols_ref)) {
      stop(sprintf("[ERPM_PLE] nodes[[%d]] columns differ from nodes[[%d]].", t, selected_partition_indices[1L]))
    }

    # Preserve original labels as a monadic covariate
    df2 <- df
    names(df2)[names(df2) == "label"] <- "label_raw"

    out <- if (is.null(out)) df2 else rbind(out, df2)
  }

  nbr_actors_meta <- nrow(out)

  # Builder requires unique, non-empty labels
  out$label <- paste0("A", seq_len(nbr_actors_meta))

  # Put label first
  out <- out[, c("label", setdiff(colnames(out), "label")), drop = FALSE]
  out
}

#' Block-diagonal bind of a timeline of square matrices (internal helper)
#' @param mats List of matrices (timeline).
#' @param idx Integer vector of indices to extract and bind.
#' @param label Label used in error messages.
#' @return Numeric block-diagonal matrix.
#' @noRd
.erpm_ple_blockdiag_mats <- function(mats, idx, label = "dyads") {
  if (!is.list(mats)) stop(sprintf("[ERPM_PLE] %s: mats must be a list.", label))
  if (!length(idx)) stop(sprintf("[ERPM_PLE] %s: idx is empty.", label))

  blocks <- lapply(idx, function(t) mats[[t]])
  if (any(vapply(blocks, is.null, logical(1)))) {
    bad <- idx[vapply(blocks, is.null, logical(1))]
    stop(sprintf("[ERPM_PLE] %s: missing matrix for t=%s.", label, paste(bad, collapse = ",")))
  }

  # Basic structural checks (square + numeric)
  dims <- lapply(blocks, dim)
  if (any(vapply(dims, function(d) length(d) != 2L, logical(1)))) {
    stop(sprintf("[ERPM_PLE] %s: some blocks have no dim().", label))
  }
  if (any(vapply(dims, function(d) d[1L] != d[2L], logical(1)))) {
    stop(sprintf("[ERPM_PLE] %s: some blocks are not square.", label))
  }
  if (any(!vapply(blocks, is.numeric, logical(1)))) {
    stop(sprintf("[ERPM_PLE] %s: all blocks must be numeric.", label))
  }

  ns <- vapply(dims, `[`, integer(1), 1L)
  N  <- sum(ns)
  out <- matrix(0, nrow = N, ncol = N)

  r0 <- 0L
  c0 <- 0L
  for (k in seq_along(blocks)) {
    n <- ns[k]
    out[(r0 + 1L):(r0 + n), (c0 + 1L):(c0 + n)] <- blocks[[k]]
    r0 <- r0 + n
    c0 <- c0 + n
  }

  out
}

#' Build meta dyads (block-diagonal over selected_partition_indices) (internal helper)
#' @param dyads_mode Either "timeline" or "meta".
#' @param dyads_data Normalized dyads container (depends on dyads_mode).
#' @param partitions Timeline of partitions (list).
#' @param selected_partition_indices Integer vector of selected times.
#' @return Named list of meta dyadic matrices, or NULL.
#' @noRd
.erpm_ple_make_meta_dyads <- function(dyads_mode, dyads_data, partitions, selected_partition_indices) {
  if (is.null(dyads_data)) return(NULL)

  nbr_actors_by_t <- vapply(partitions, length, integer(1))
  nbr_actors_meta <- sum(nbr_actors_by_t[selected_partition_indices])

  if (identical(dyads_mode, "timeline")) {
    out <- list()
    for (nm in names(dyads_data)) {
      mats <- dyads_data[[nm]]

      # Dimension checks per selected time
      for (t in selected_partition_indices) {
        M <- mats[[t]]
        nbr_actors_in_block <- nbr_actors_by_t[t]
        if (!is.matrix(M) || nrow(M) != nbr_actors_in_block || ncol(M) != nbr_actors_in_block) {
          stop(sprintf("[ERPM_PLE] dyads '%s' at t=%d must be %dx%d.", nm, t, nbr_actors_in_block, nbr_actors_in_block))
        }
      }

      out[[nm]] <- .erpm_ple_blockdiag_mats(mats, idx = selected_partition_indices, label = nm)
    }
    return(out)
  }

  if (identical(dyads_mode, "meta")) {
    out <- list()
    for (nm in names(dyads_data)) {
      M <- dyads_data[[nm]]
      if (!is.matrix(M) || nrow(M) != nbr_actors_meta || ncol(M) != nbr_actors_meta) {
        stop(sprintf("[ERPM_PLE] meta dyads '%s' must be %dx%d (nbr_actors_meta=%d).",
                     nm, nbr_actors_meta, nbr_actors_meta, nbr_actors_meta))
      }
      out[[nm]] <- M
    }
    return(out)
  }

  stop("[ERPM_PLE] Internal error: unknown dyads_mode.")
}

# ==============================================================================
# Module 1: build standard bipartite meta-network
# ==============================================================================

#' Build a standard bipartite meta-network (PLE module 1)
#'
#' @description
#' Aggregates partitions/nodes/dyads across selected blocks and builds a bipartite
#' meta-network using \code{build_bipartite_from_inputs()}.
#'
#' @param partitions List of partitions (timeline).
#' @param selected_partition_indices Integer vector of selected time indices included in meta-network.
#' @param nodes NULL or list length T of per-time node data.frames.
#' @param dyads NULL or dyadic timeline input (current behavior preserved).
#' @param group_labels Group labels (currently not used in meta build; preserved for API).
#' @param verbose Verbosity flag.
#'
#' @return List containing \code{meta_nw} plus derived objects used downstream.
#' @noRd
.erpm_long_empile_build_standard_meta_network <- function(partitions,
                                                         selected_partition_indices,
                                                         nodes,
                                                         dyads,
                                                         group_labels,
                                                         verbose) {
  # NOTE: dyads normalization behavior is preserved; only the comment structure is updated.
  nodes_n <- nodes

  dyads_mode <- if (is.null(dyads)) NULL else "timeline"
  dyads_data <- NULL
  if (!is.null(dyads)) {
    nms <- unique(unlist(lapply(dyads, names)))
    dyads_data <- setNames(vector("list", length(nms)), nms)
    for (nm in nms) dyads_data[[nm]] <- lapply(dyads, `[[`, nm)
  }

  if (isTRUE(verbose)) {
    T <- length(partitions)
    message(sprintf("[ERPM_PLE] T=%d | meta blocks=%s",
                    T, paste(selected_partition_indices, collapse = ",")))
  }

  # Meta partition (concatenate + shift)
  mp <- .erpm_ple_make_meta_partition(partitions, selected_partition_indices= selected_partition_indices)
  meta_partition <- mp$meta_partition
  nbr_actors_by_t        <- mp$nbr_actors_by_t
  actor_start_index_by_t   <- mp$actor_start_index_by_t
  nbr_actors_by_selected_partitions         <- mp$nbr_actors_by_selected_partitions
  actor_start_index_by_block    <- mp$actor_start_index_by_block
  nbr_actors_meta        <- sum(nbr_actors_by_selected_partitions)

  # Meta nodes / dyads
  meta_nodes_df <- .erpm_ple_make_meta_nodes_df(nodes_n, partitions, selected_partition_indices= selected_partition_indices)
  meta_dyads    <- .erpm_ple_make_meta_dyads(dyads_mode, dyads_data, partitions, selected_partition_indices= selected_partition_indices)

  # Canonical bipartite build
  built_meta <- build_bipartite_from_inputs(
    partition    = meta_partition,
    nodes        = meta_nodes_df,
    dyads        = if (is.null(meta_dyads)) list() else meta_dyads,
    group_labels = NULL
  )
  meta_nw <- built_meta$network

  list(
    meta_nw        = meta_nw,
    meta_partition = meta_partition,
    nbr_actors_meta        = nbr_actors_meta,
    nbr_actors_by_t        = nbr_actors_by_t,
    actor_start_index_by_t   = actor_start_index_by_t,
    nbr_actors_by_selected_partitions         = nbr_actors_by_selected_partitions,
    actor_start_index_by_block    = actor_start_index_by_block,
    meta_nodes_df  = meta_nodes_df,
    meta_dyads     = meta_dyads,
    dyads_mode     = dyads_mode,
    nodes_n        = nodes_n,
    dyads_data     = dyads_data
  )
}

#' Sanity checks for PLE meta-network actor/group construction (internal helper)
#'
#' This is a defensive validator for the meta-network plumbing:
#' sizes, labels, membership edges, per-block disjointness, node attributes, and dyads.
#'
#' @param built Result of .erpm_long_empile_build_standard_meta_network().
#' @param debug Logical. If TRUE, print detailed state to console.
#' @return TRUE (invisibly) if all checks pass; otherwise stops with an error.
#' @noRd
.erpm_ple_check_meta_network_actors_groups <- function(built, debug = FALSE) {
  .dbg <- function(...) if (isTRUE(debug)) message(...)

  # Unpack
  meta_nw  <- built$meta_nw
  meta_partition <- built$meta_partition
  nbr_actors_by_selected_partitions <- built$nbr_actors_by_selected_partitions
  actor_start_index_by_block <- built$actor_start_index_by_block

  meta_nodes_df <- built$meta_nodes_df
  meta_dyads    <- built$meta_dyads
  dyads_mode    <- built$dyads_mode
  nodes_n       <- built$nodes_n
  dyads_data    <- built$dyads_data

  .dbg("[ERPM_PLE][CHECK][DEBUG] ---- built keys ----")
  .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] names(built)=%s", paste(names(built), collapse = ", ")))

  # Helpers
  .stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

  .as_char_head <- function(x, n = 20L) {
    if (is.factor(x)) x <- as.character(x)
    if (is.list(x)) {
      x <- vapply(x, function(z) paste0(z, collapse = "|"), character(1))
    }
    paste(utils::head(x, n), collapse = ",")
  }

  # Compare dyads robustly: ignore attributes, coerce to double, allow tiny tolerance
  .dyads_equal <- function(A, B, tol = 0) {
    if (is.null(A) || is.null(B)) return(FALSE)
    if (!is.matrix(A) || !is.matrix(B)) return(FALSE)
    if (!identical(dim(A), dim(B))) return(FALSE)

    storage.mode(A) <- "double"
    storage.mode(B) <- "double"
    dimnames(A) <- NULL
    dimnames(B) <- NULL

    if (identical(A, B)) return(TRUE) # fast path

    ok <- isTRUE(all.equal(A, B, check.attributes = FALSE, tolerance = tol))
    if (ok) return(TRUE)

    fin <- is.finite(A) & is.finite(B)
    if (!any(fin)) return(TRUE)
    maxdiff <- max(abs(A[fin] - B[fin]))
    maxdiff <= tol
  }

  .dyads_diff_max <- function(A, B) {
    if (is.null(A) || is.null(B) || !is.matrix(A) || !is.matrix(B) || !identical(dim(A), dim(B))) return(NA_real_)
    storage.mode(A) <- "double"
    storage.mode(B) <- "double"
    dimnames(A) <- NULL
    dimnames(B) <- NULL
    fin <- is.finite(A) & is.finite(B)
    if (!any(fin)) return(0)
    max(abs(A[fin] - B[fin]))
  }

  # Global sizes
  nA <- meta_nw %n% "bipartite"
  if (is.na(nA)) .stopf("[ERPM_PLE][CHECK] meta_nw has no 'bipartite' attribute.")

  expected_nA <- sum(nbr_actors_by_selected_partitions)

  .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] nA(bipartite)=%d | expected_nA=%d", nA, expected_nA))
  .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] network.size=%d | 2*nA=%d",
               network::network.size(meta_nw), 2L * nA))
  .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] edgecount=%d | expected=%d",
               network::network.edgecount(meta_nw), nA))

  if (nA != expected_nA) {
    .stopf("[ERPM_PLE][CHECK] bipartite=%d but expected %d actors.", nA, expected_nA)
  }
  if (network::network.size(meta_nw) != 2L * nA) {
    .stopf("[ERPM_PLE][CHECK] meta_nw size != 2 * bipartite.")
  }
  if (network::network.edgecount(meta_nw) != nA) {
    .stopf("[ERPM_PLE][CHECK] meta_nw must have exactly one edge per actor.")
  }

  # Actor / group indices
  actors_idx <- seq_len(nA)
  groups_idx <- (nA + 1L):(2L * nA)

  vnames <- network::get.vertex.attribute(meta_nw, "vertex.names")
  if (isTRUE(debug)) {
    .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] head(vnames,20)=%s", paste(utils::head(vnames, 20), collapse = ", ")))
    if (length(vnames) >= 2L * nA) {
      .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] actor vnames[1:%d]=%s",
                   min(nA, 20L), paste(vnames[utils::head(actors_idx, 20L)], collapse = ", ")))
      .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] group vnames[%d:%d]=%s",
                   nA + 1L, nA + min(nA, 20L),
                   paste(vnames[groups_idx[seq_len(min(nA, 20L))]], collapse = ", ")))
    }
  }

  if (!all(vnames[actors_idx] == paste0("A", actors_idx))) {
    .stopf("[ERPM_PLE][CHECK] actor labels are not A1..An.")
  }
  if (!all(vnames[groups_idx] == paste0("G", seq_len(nA)))) {
    .stopf("[ERPM_PLE][CHECK] group labels are not G1..Gn.")
  }

  # Meta partition consistency
  .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] length(meta_partition)=%d | nA=%d", length(meta_partition), nA))
  if (isTRUE(debug)) {
    .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] meta_partition head=%s", paste(utils::head(meta_partition, 20L), collapse = ",")))
    .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] meta_partition range=[%d,%d]",
                 suppressWarnings(min(meta_partition)), suppressWarnings(max(meta_partition))))
  }

  if (length(meta_partition) != nA) {
    .stopf("[ERPM_PLE][CHECK] meta_partition length != number of actors.")
  }
  if (any(meta_partition < 1L | meta_partition > nA)) {
    .stopf("[ERPM_PLE][CHECK] meta_partition contains invalid group ids.")
  }

  # Edge list consistency: actor i -> group nA + meta_partition[i]
  el <- network::as.edgelist(meta_nw)
  if (isTRUE(debug)) {
    .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] edgelist dim=%dx%d", nrow(el), ncol(el)))
    .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] edgelist head=%s",
                 paste(apply(utils::head(el, 10L), 1L, function(r) paste0(r[1L], "->", r[2L])), collapse = " | ")))
    .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] all(el[,1]<=nA)=%s", as.character(all(el[,1L] <= nA))))
    .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] all(el[,2]>nA)=%s", as.character(all(el[,2L] > nA))))
  }

  if (!all(el[, 1L] <= nA)) .stopf("[ERPM_PLE][CHECK] some edges do not start from actor nodes.")
  if (!all(el[, 2L] > nA))  .stopf("[ERPM_PLE][CHECK] some edges do not point to group nodes.")

  if (!all(el[, 2L] == nA + meta_partition[el[, 1L]])) {
    if (isTRUE(debug)) {
      bad <- which(el[, 2L] != nA + meta_partition[el[, 1L]])
      bad <- utils::head(bad, 20L)
      .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] bad edges idx=%s", paste(bad, collapse = ",")))
      .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] bad edges=%s",
                   paste(apply(el[bad, , drop = FALSE], 1L, function(r) paste0(r[1L], "->", r[2L])), collapse = " | ")))
    }
    .stopf("[ERPM_PLE][CHECK] edge list inconsistent with meta_partition.")
  }

  # Block structure: group-id ranges must be disjoint across blocks
  if (isTRUE(debug)) {
    .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] B=%d blocks | nbr_actors_by_selected_partitions=%s",
                 length(nbr_actors_by_selected_partitions),
                 paste(nbr_actors_by_selected_partitions, collapse = ",")))
    .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] actor_start_index_by_block=%s",
                 paste(actor_start_index_by_block, collapse = ",")))
  }

  for (b in seq_along(nbr_actors_by_selected_partitions)) {
    a0 <- actor_start_index_by_block[b]
    n  <- nbr_actors_by_selected_partitions[b]
    idx <- (a0 + 1L):(a0 + n)
    g   <- meta_partition[idx]

    if (isTRUE(debug)) {
      .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] block %d: actors %d..%d | group_ids=%s",
                   b, a0 + 1L, a0 + n, paste(g, collapse = ",")))
    }

    if (length(intersect(g, meta_partition[-idx])) > 0L) {
      .stopf("[ERPM_PLE][CHECK] group id overlap detected for block %d.", b)
    }
  }

  # Nodes: conformity checks (only if present)
  if (isTRUE(debug)) {
    .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] nodes_n class=%s | meta_nodes_df class=%s",
                 paste(class(nodes_n), collapse = "/"),
                 if (is.null(meta_nodes_df)) "NULL" else paste(class(meta_nodes_df), collapse = "/")))
  }

  if (!is.null(meta_nodes_df)) {
    if (!is.data.frame(meta_nodes_df)) .stopf("[ERPM_PLE][CHECK] meta_nodes_df must be a data.frame when non-NULL.")
    if (nrow(meta_nodes_df) != nA) {
      .stopf("[ERPM_PLE][CHECK] nrow(meta_nodes_df)=%d but expected %d (one row per actor).", nrow(meta_nodes_df), nA)
    }
    if (!("label" %in% names(meta_nodes_df))) .stopf("[ERPM_PLE][CHECK] meta_nodes_df must contain a 'label' column.")

    lab <- as.character(meta_nodes_df$label)
    if (anyNA(lab) || any(!nzchar(lab))) .stopf("[ERPM_PLE][CHECK] meta_nodes_df$label contains NA/empty labels.")

    if (isTRUE(debug)) {
      .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] meta_nodes_df nrow=%d ncol=%d cols=%s",
                   nrow(meta_nodes_df), ncol(meta_nodes_df), paste(names(meta_nodes_df), collapse = ",")))
      for (col in names(meta_nodes_df)) {
        .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] meta_nodes_df$%s head=%s", col, .as_char_head(meta_nodes_df[[col]], 20L)))
      }
    }

    node_cols <- setdiff(names(meta_nodes_df), "label")
    if (length(node_cols)) {
      for (a in node_cols) {
        va <- network::get.vertex.attribute(meta_nw, a)
        if (is.null(va)) .stopf("[ERPM_PLE][CHECK] actor attribute '%s' missing on meta_nw.", a)

        if (!isTRUE(all.equal(va[actors_idx], meta_nodes_df[[a]], check.attributes = FALSE))) {
          if (isTRUE(debug)) {
            .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] mismatch attr='%s' actor-side", a))
            .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] meta_nodes_df head=%s", .as_char_head(meta_nodes_df[[a]], 20L)))
            .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] meta_nw attr head=%s", .as_char_head(va[actors_idx], 20L)))
          }
          .stopf("[ERPM_PLE][CHECK] actor attribute '%s' differs from meta_nodes_df.", a)
        }

        if (any(!is.na(va[groups_idx]))) {
          if (isTRUE(debug)) {
            bad <- utils::head(which(!is.na(va[groups_idx])), 20L)
            .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] non-NA padding detected for attr='%s' at groups idx=%s",
                         a, paste(bad, collapse = ",")))
          }
          .stopf("[ERPM_PLE][CHECK] group-side padding for node attribute '%s' is not NA.", a)
        }
      }
    }
  }

  # Dyads: conformity checks (only if present)
  dyads_att <- tryCatch(meta_nw %n% "dyads", error = function(e) NULL)

  .dbg_block <- function(lines, prefix = "[ERPM_PLE][CHECK][DEBUG] ") {
    if (!isTRUE(debug)) return(invisible(NULL))
    if (!length(lines)) return(invisible(NULL))
    for (ln in lines) .dbg(paste0(prefix, ln))
    invisible(NULL)
  }

  .dbg_matrix <- function(M, title, max_print = Inf) {
    if (!isTRUE(debug)) return(invisible(NULL))
    if (is.null(M)) {
      .dbg(paste0("[ERPM_PLE][CHECK][DEBUG] ", title, " = NULL"))
      return(invisible(NULL))
    }
    if (!is.matrix(M)) {
      .dbg(paste0("[ERPM_PLE][CHECK][DEBUG] ", title, " (not a matrix) class=", paste(class(M), collapse = "/")))
      return(invisible(NULL))
    }

    oldw <- getOption("width")
    oldm <- getOption("max.print")
    options(width = max(200L, oldw), max.print = if (is.finite(max_print)) max_print else 1e9)
    on.exit(options(width = oldw, max.print = oldm), add = TRUE)

    out <- capture.output(print(M))
    .dbg(paste0("[ERPM_PLE][CHECK][DEBUG] ", title, " (", nrow(M), "x", ncol(M), ")"))
    .dbg_block(out)
    invisible(NULL)
  }

  if (isTRUE(debug)) {
    .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] dyads_mode=%s", if (is.null(dyads_mode)) "NULL" else as.character(dyads_mode)))
    .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] dyads_data is %s", if (is.null(dyads_data)) "NULL" else "non-NULL"))
    .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] meta_dyads is %s | nw%%n%%'dyads' is %s",
                if (is.null(meta_dyads)) "NULL" else "non-NULL",
                if (is.null(dyads_att)) "NULL" else "non-NULL"))
    if (!is.null(meta_dyads) && length(meta_dyads)) .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] meta_dyads names=%s", paste(names(meta_dyads), collapse = ",")))
    if (!is.null(dyads_att) && length(dyads_att))   .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] attached dyads names=%s", paste(names(dyads_att), collapse = ",")))
  }

  if (!is.null(meta_dyads)) {
    if (!is.list(meta_dyads)) .stopf("[ERPM_PLE][CHECK] meta_dyads must be a list when non-NULL.")
    if (length(meta_dyads) && (is.null(names(meta_dyads)) || any(!nzchar(names(meta_dyads))))) {
      .stopf("[ERPM_PLE][CHECK] meta_dyads must be a NAMED list (e.g., fm=..., Z1=...).")
    }

    if (is.null(dyads_att)) .stopf("[ERPM_PLE][CHECK] meta_dyads is non-NULL but meta_nw %%n%% 'dyads' is NULL.")
    if (!is.list(dyads_att)) .stopf("[ERPM_PLE][CHECK] meta_nw %%n%% 'dyads' must be a list.")
    if (!setequal(names(dyads_att), names(meta_dyads))) {
      .stopf("[ERPM_PLE][CHECK] meta_nw %%n%% 'dyads' names differ from meta_dyads names.")
    }

    for (nm in names(meta_dyads)) {
      M <- meta_dyads[[nm]]
      A <- dyads_att[[nm]]

      if (!is.matrix(M))  .stopf("[ERPM_PLE][CHECK] meta_dyads[['%s']] must be a matrix.", nm)
      if (!is.numeric(M)) .stopf("[ERPM_PLE][CHECK] meta_dyads[['%s']] must be numeric.", nm)
      if (nrow(M) != nA || ncol(M) != nA) {
        .stopf("[ERPM_PLE][CHECK] meta_dyads[['%s']] has dim %dx%d but expected %dx%d.", nm, nrow(M), ncol(M), nA, nA)
      }
      if (any(!is.finite(M))) .stopf("[ERPM_PLE][CHECK] meta_dyads[['%s']] contains non-finite values.", nm)

      if (is.null(A) || !is.matrix(A)) {
        .stopf("[ERPM_PLE][CHECK] meta_nw %%n%% 'dyads'[['%s']] is missing or not a matrix.", nm)
      }
      if (!identical(dim(A), c(nA, nA))) {
        .stopf("[ERPM_PLE][CHECK] meta_nw %%n%% 'dyads'[['%s']] dim %dx%d but expected %dx%d.",
              nm, nrow(A), ncol(A), nA, nA)
      }

      if (isTRUE(debug)) {
        .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] dyads '%s': dim=%dx%d | finite=%s",
                    nm, nrow(M), ncol(M), as.character(all(is.finite(M)))))
        .dbg_matrix(A, paste0("meta_nw %n% 'dyads'[['", nm, "']]"))
      }

      if (!.dyads_equal(A, M, tol = 0)) {
        if (isTRUE(debug)) {
          .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] dyads mismatch '%s': attached vs meta_dyads", nm))
          .dbg(sprintf("[ERPM_PLE][CHECK][DEBUG] max|diff|=%s", as.character(.dyads_diff_max(A, M))))

          AA <- A; MM <- M
          storage.mode(AA) <- "double"; storage.mode(MM) <- "double"
          dimnames(AA) <- NULL; dimnames(MM) <- NULL
          D <- abs(AA - MM)
          D[!is.finite(D)] <- NA_real_
          if (any(D > 0, na.rm = TRUE)) {
            ord <- order(D, decreasing = TRUE, na.last = NA)
            ord <- utils::head(ord, 20L)
            ij  <- arrayInd(ord, dim(D))
            lines <- vapply(seq_len(nrow(ij)), function(k) {
              i <- ij[k, 1L]; j <- ij[k, 2L]
              sprintf("diff[%d,%d]=%g (attached=%g, meta=%g)", i, j, D[i,j], AA[i,j], MM[i,j])
            }, character(1))
            .dbg_block(lines)
          }
        }
        .stopf("[ERPM_PLE][CHECK] meta_nw %%n%% 'dyads'[['%s']] differs from meta_dyads.", nm)
      }
    }
  } else {
    if (!is.null(dyads_att) && length(dyads_att)) {
      .stopf("[ERPM_PLE][CHECK] meta_dyads is NULL but meta_nw %%n%% 'dyads' is non-empty.")
    }
  }

  # Input-mode invariants (light checks to catch plumbing mistakes)
  if (is.null(nodes_n) && !is.null(meta_nodes_df)) {
    .stopf("[ERPM_PLE][CHECK] meta_nodes_df is non-NULL but built$nodes_n is NULL (unexpected).")
  }

  if (is.null(dyads_mode) && !is.null(meta_dyads)) {
    .stopf("[ERPM_PLE][CHECK] meta_dyads is non-NULL but dyads_mode is NULL (unexpected).")
  }

  if (!is.null(dyads_mode) && !identical(dyads_mode, "timeline")) {
    .stopf("[ERPM_PLE][CHECK] dyads_mode='%s' (expected NULL or 'timeline').", as.character(dyads_mode))
  }

  if (!is.null(dyads_data)) {
    if (!is.list(dyads_data) || (length(dyads_data) && (is.null(names(dyads_data)) || any(!nzchar(names(dyads_data)))))) {
      .stopf("[ERPM_PLE][CHECK] dyads_data must be a named list when non-NULL.")
    }
  }

  .dbg("[ERPM_PLE][CHECK][DEBUG] OK")
  invisible(TRUE)
}

# ==============================================================================
# Module 2: build and attach timeblock
# ==============================================================================

#' Build timeblock for PLE meta-networks (internal helper)
#'
#' @description
#' The meta-network has \code{2 * nbr_actors_meta} vertices: actors
#' (1..nbr_actors_meta) and padded groups (nbr_actors_meta+1 .. 2*nbr_actors_meta).
#' This helper assigns the originating time index to each actor and each padded group.
#'
#' @param selected_partition_indices Selected time indices included in meta-network.
#' @param nbr_actors_by_t Actor counts per time.
#' @return Integer vector of length \code{2 * nbr_actors_meta}.
#' @noRd
.erpm_ple_make_timeblock <- function(selected_partition_indices, nbr_actors_by_t) {
  nbr_actors_by_selected_partitions  <- nbr_actors_by_t[selected_partition_indices]
  nbr_actors_meta <- sum(nbr_actors_by_selected_partitions)

  tb_actor <- integer(nbr_actors_meta)
  tb_group <- integer(nbr_actors_meta)

  a0 <- 0L
  g0 <- 0L

  for (b in seq_along(selected_partition_indices)) {
    t  <- selected_partition_indices[b]
    nbr_actors_in_block <- nbr_actors_by_selected_partitions[b]

    tb_actor[(a0 + 1L):(a0 + nbr_actors_in_block)] <- t
    tb_group[(g0 + 1L):(g0 + nbr_actors_in_block)] <- t

    a0 <- a0 + nbr_actors_in_block
    g0 <- g0 + nbr_actors_in_block
  }

  c(tb_actor, tb_group)
}

#' Attach timeblock to meta-network (PLE module 2)
#' @param meta_nw Meta-network (bipartite).
#' @param selected_partition_indices Selected time indices included in meta-network.
#' @param nbr_actors_by_t Actor counts per time.
#' @return Updated meta-network with \code{timeblock} vertex attribute.
#' @noRd
.erpm_long_empile_attach_timeblock <- function(meta_nw, selected_partition_indices, nbr_actors_by_t) {
  timeblock <- .erpm_ple_make_timeblock(selected_partition_indices, nbr_actors_by_t = nbr_actors_by_t)
  network::set.vertex.attribute(meta_nw, "timeblock", timeblock)
  meta_nw
}

# ==============================================================================
# Module 3: attach inertial attributes
# ==============================================================================

#' Attach inertial attributes required by inertial InitErgmTerms (PLE module 3)
#'
#' @description
#' This module currently targets \code{inertia_groups} in PLE mode.
#' It attaches the exact network attributes expected by your InitErgmTerm:
#' \code{erpm_mode}, \code{erpm_B}, \code{erpm_n}, \code{erpm_G},
#' and \code{erpm_block_past_partitions}.
#'
#' @param meta_nw Meta-network to decorate.
#' @param partitions Timeline of partitions (list).
#' @param selected_partition_indices Selected time indices (estimation blocks).
#' @param d Past influence depth.
#' @param verbose Verbosity flag.
#' @return Updated meta-network with inertial attributes attached.
#' @noRd
.erpm_long_empile_attach_inertial_attributes <- function(meta_nw,
                                                        partitions,
                                                        selected_partition_indices,
                                                        d,
                                                        verbose) {
  if (d < 1L) stop("[ERPM_PLE] inertia_groups requires past_influence >= 1.")

  nbr_actors_per_selected_partition <- .erpm_ple_assert_constant_nbr_actors_per_selected_partition(partitions, selected_partition_indices= selected_partition_indices)

  B <- length(selected_partition_indices)
  G_block <- nbr_actors_per_selected_partition

  # Sanity check: bipartite size must be nbr_actors_per_selected_partition * B
  n1_total <- as.integer(network::get.network.attribute(meta_nw, "bipartite"))
  if (is.na(n1_total) || n1_total != nbr_actors_per_selected_partition * B) {
    stop(sprintf(
      "[ERPM_PLE] inconsistent bipartite size for inertia_groups: bipartite=%s but nbr_actors_per_selected_partition*B=%d*%d=%d.",
      as.character(n1_total), nbr_actors_per_selected_partition, B, nbr_actors_per_selected_partition * B
    ))
  }

  erpm_block_past_partitions <- .erpm_ple_make_erpm_block_past_partitions(
    partitions = partitions,
    selected_partition_indices   = selected_partition_indices,
    d          = d,
    nbr_actors_per_selected_partition    = nbr_actors_per_selected_partition
  )

  # Attributes expected by InitErgmTerm.inertia_groups.R (PLE mode)
  network::set.network.attribute(meta_nw, "erpm_mode", "empile")
  network::set.network.attribute(meta_nw, "erpm_B", B)
  network::set.network.attribute(meta_nw, "erpm_n", nbr_actors_per_selected_partition)
  network::set.network.attribute(meta_nw, "erpm_G", G_block)
  network::set.network.attribute(meta_nw, "erpm_block_past_partitions", erpm_block_past_partitions)

  if (isTRUE(verbose)) {
    message(sprintf("[ERPM_PLE] inertia_groups attrs attached: erpm_mode=empile | erpm_B=%d | erpm_n=%d | erpm_G=%d",
                    B, nbr_actors_per_selected_partition, G_block))
  }

  meta_nw
}

# ==============================================================================
# Orchestrator: called by erpm_long()
# ==============================================================================

#' Build PLE meta-network and attach required attributes (engine entry point)
#'
#' @description
#' Orchestrates the PLE engine pipeline:
#' \enumerate{
#'   \item build standard bipartite meta-network (partitions/nodes/dyads);
#'   \item attach \code{timeblock} vertex attribute;
#'   \item attach inertial attributes if inertial terms were detected.
#' }
#'
#' @param partitions List of partitions.
#' @param rhs RHS expression (passed through in the return bundle).
#' @param inertial_present Logical: inertial terms detected upstream.
#' @param past_influence Integer past influence (d).
#' @param nodes Per-time nodes input.
#' @param dyads Dyadic inputs.
#' @param group_labels Group labels (kept for API compatibility).
#' @param directed Kept for compatibility.
#' @param verbose Verbosity flag.
#' @param debug Debug flag.
#'
#' @return List with at least: \code{meta_nw}, \code{selected_partition_indices}, \code{rhs}.
#' @noRd
.erpm_long_empile_build_meta_nw <- function(partitions,
                                           rhs,
                                           inertial_present = FALSE,
                                           past_influence = 0L,
                                           nodes = NULL,
                                           dyads = NULL,
                                           group_labels = NULL,
                                           directed = FALSE,
                                           verbose = FALSE,
                                           debug   = FALSE) {

  if (!is.list(partitions) || !length(partitions)) {
    stop("[ERPM_PLE] partitions must be a non-empty list.")
  }

  T <- length(partitions)
  d <- as.integer(past_influence)

  if (is.na(d) || d < 0L) stop("[ERPM_PLE] past_influence must be a non-negative integer.")

  if (inertial_present && d >= T) {
    stop(sprintf("[ERPM_PLE] past_influence=%d but T=%d: cannot build estimable meta-network.", d, T))
  }

  # Blocks included in the meta-network
  selected_partition_indices<- if (inertial_present) seq.int(d + 1L, T) else seq_len(T)

  .erpm_long_vcat(verbose, sprintf("[ERPM_PLE] selected_partition_indices=%s",
                                  paste(selected_partition_indices, collapse = ",")))
  .erpm_long_dcat(debug, sprintf("[ERPM_PLE][DEBUG] T=%d | d=%d | inertial_present=%s",
                                T, d, as.character(inertial_present)))
  .erpm_long_dcat(debug, sprintf("[ERPM_PLE][DEBUG] nodes=%s | dyads=%s",
                                if (is.null(nodes)) "NULL" else "non-NULL",
                                if (is.null(dyads)) "NULL" else "non-NULL"))

  # ---------------------------------------------------------------------------
  # 1) Standard meta-network build
  # ---------------------------------------------------------------------------
  b <- .erpm_long_empile_build_standard_meta_network(
    partitions   = partitions,
    selected_partition_indices     = selected_partition_indices,
    nodes        = nodes,
    dyads        = dyads,
    group_labels = group_labels,
    verbose      = verbose
  )

  .erpm_long_vcat(verbose, sprintf("[ERPM_PLE] meta built: bipartite(nA)=%d | size=%d | edges=%d",
                                  as.integer(b$meta_nw %n% "bipartite"),
                                  network::network.size(b$meta_nw),
                                  network::network.edgecount(b$meta_nw)))
  .erpm_long_dcat(debug, sprintf("[ERPM_PLE][DEBUG] dyads_mode=%s",
                                if (is.null(b$dyads_mode)) "NULL" else as.character(b$dyads_mode)))
  if (!is.null(b$meta_dyads)) {
    .erpm_long_dcat(debug, sprintf("[ERPM_PLE][DEBUG] meta_dyads names=%s", paste(names(b$meta_dyads), collapse=",")))
  }

  .erpm_ple_check_meta_network_actors_groups(b, debug=debug)

  meta_nw <- b$meta_nw

  # Keep a consistent marker even without inertial terms
  if (!inertial_present) {
    network::set.network.attribute(meta_nw, "erpm_mode", "empile")
  }

  # ---------------------------------------------------------------------------
  # 2) Attach timeblock
  # ---------------------------------------------------------------------------
  meta_nw <- .erpm_long_empile_attach_timeblock(
    meta_nw = meta_nw,
    selected_partition_indices= selected_partition_indices,
    nbr_actors_by_t = b$nbr_actors_by_t
  )

  # ---------------------------------------------------------------------------
  # 3) Attach inertial attributes (if requested)
  # ---------------------------------------------------------------------------
  if (inertial_present) {
    meta_nw <- .erpm_long_empile_attach_inertial_attributes(
      meta_nw    = meta_nw,
      partitions = partitions,
      selected_partition_indices   = selected_partition_indices,
      d          = d,
      verbose    = verbose
    )
  }

  # ---------------------------------------------------------------------------
  # Bookkeeping attributes (existing behavior preserved)
  # ---------------------------------------------------------------------------
  network::set.network.attribute(meta_nw, "erpm_long.mode", "PLE")
  network::set.network.attribute(meta_nw, "erpm_long.T", T)
  network::set.network.attribute(meta_nw, "erpm_long.d", if (inertial_present) d else 0L)
  network::set.network.attribute(meta_nw, "erpm_long.selected_partition_indices", selected_partition_indices)

  nbr_actors_by_selected_partitions <- b$nbr_actors_by_t[selected_partition_indices]
  actor_offsets <- c(0L, cumsum(nbr_actors_by_selected_partitions))[seq_along(nbr_actors_by_selected_partitions)]

  network::set.network.attribute(meta_nw, "erpm_long.nbr_actors_meta", b$nbr_actors_meta)
  network::set.network.attribute(meta_nw, "erpm_long.nbr_actors_by_t", b$nbr_actors_by_t)
  network::set.network.attribute(meta_nw, "erpm_long.actor_start_index_by_t", b$actor_start_index_by_t)
  network::set.network.attribute(meta_nw, "erpm_long.actor_offsets", actor_offsets)
  network::set.network.attribute(meta_nw, "erpm_long.group_id_actor_start_index_by_t", b$actor_start_index_by_t)

  if (!is.null(b$meta_nodes_df)) {
    nm_nodes <- setdiff(colnames(b$meta_nodes_df), "label")
    network::set.network.attribute(meta_nw, "erpm_long.meta_nodes_names", nm_nodes)
  }

  if (!is.null(b$meta_dyads)) {
    network::set.network.attribute(meta_nw, "erpm_long.meta_dyads_names", names(b$meta_dyads))
    network::set.network.attribute(meta_nw, "erpm_long.dyads_mode", b$dyads_mode)
  }

  list(
    meta_nw      = meta_nw,
    selected_partition_indices     = selected_partition_indices,
    timeline_nws = if (inertial_present) network::get.network.attribute(meta_nw, "erpm_long.timeline_nws") else NULL,
    rhs          = rhs
  )
}