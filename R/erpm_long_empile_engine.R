################################################################################
# FILE: R/erpm_long_empile_engine.R
################################################################################

#' ERPM PLE engine: build stacked (empile) meta-network and attach attributes
#'
#' @name erpm_long_empile_engine
#' @note erpm_long_empile_engine.R
#'
#' @description
#' This module provides the PLE-only engine used by \code{erpm_long()}.
#' The engine is intentionally split into three explicit modules:
#' \enumerate{
#'   \item build a standard bipartite meta-network by aggregating partitions/nodes/dyads;
#'   \item build and attach a \code{timeblock} vertex attribute;
#'   \item if inertial terms are present, attach additional network attributes
#'         required by inertial InitErgmTerms (e.g., \code{inertia_groups}).
#' }
#'
#' The current code preserves your existing behavior; only structure and
#' nomenclature are targeted here.
#'
#' @keywords ERPM ERGM PLE empile engine

# ==============================================================================
# Low-level helpers (existing behavior)
# ==============================================================================

#' Enforce constant actor block size across estimation blocks (internal helper)
#' @noRd
.erpm_ple_assert_constant_n_block <- function(partitions, idx_est) {
  nA_by_t <- vapply(partitions, length, integer(1))
  n_block <- nA_by_t[idx_est[1L]]
  if (any(nA_by_t[idx_est] != n_block)) {
    stop(sprintf(
      "[ERPM_PLE] PLE(inertia_groups) requires constant n_block across idx_est. Got nA_by_t[idx_est]=%s.",
      paste(nA_by_t[idx_est], collapse = ",")
    ))
  }
  n_block
}

#' Build per-block past partitions container expected by inertia_groups (internal helper)
#' @noRd
.erpm_ple_make_erpm_block_past_partitions <- function(partitions, idx_est, d, n_block) {
  B <- length(idx_est)
  out <- vector("list", B)

  for (b in seq_len(B)) {
    t_cur <- idx_est[b]
    pb <- vector("list", d)

    for (lag in seq_len(d)) {
      t_past <- t_cur - lag
      if (t_past < 1L) {
        stop(sprintf("[ERPM_PLE] internal: t_past=%d < 1 (b=%d, lag=%d). idx_est must start at d+1.",
                     t_past, b, lag))
      }

      p <- partitions[[t_past]]

      if (is.null(p) || !is.atomic(p) || length(p) != n_block) {
        stop(sprintf("[ERPM_PLE] past partition invalid at time=%d (expected atomic length n_block=%d).",
                     t_past, n_block))
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
#' @noRd
.erpm_ple_G_obs <- function(partition) {
  f <- factor(partition, levels = unique(partition))
  nlevels(f)
}

#' Build concatenated+shifted meta partition (internal helper)
#' @noRd
.erpm_ple_make_meta_partition <- function(partitions, idx_est) {
  nA_by_t <- vapply(partitions, length, integer(1))

  nA_est <- nA_by_t[idx_est]
  offsets_est <- c(0L, cumsum(nA_est))[seq_along(nA_est)]  # length B

  # offsets_by_t[t] is defined for t in idx_est, NA otherwise
  offsets_by_t <- rep(NA_integer_, length(nA_by_t))
  offsets_by_t[idx_est] <- offsets_est

  meta <- integer(0)

  for (b in seq_along(idx_est)) {
    t <- idx_est[b]
    p <- partitions[[t]]

    # Stable relabeling within block: groups become 1..G_obs
    f <- factor(p, levels = unique(p))
    g <- as.integer(f)

    # Shift group ids so each block uses disjoint group-id ranges
    meta <- c(meta, g + offsets_est[b])
  }

  list(
    meta_partition = meta,
    nA_by_t        = nA_by_t,
    idx_est        = idx_est,
    nA_est         = nA_est,
    offsets_est    = offsets_est,
    offsets_by_t   = offsets_by_t
  )
}

#' Build meta nodes data.frame for build_bipartite_from_inputs (internal helper)
#' @noRd
.erpm_ple_make_meta_nodes_df <- function(nodes, partitions, idx_est) {
  if (is.null(nodes)) return(NULL)

  T <- length(partitions)
  if (!is.list(nodes) || length(nodes) != T) {
    stop("[ERPM_PLE] nodes must be NULL or a list of length T (one data.frame per partition).")
  }

  cols_ref <- NULL
  out <- NULL

  for (t in idx_est) {
    df <- nodes[[t]]
    if (!is.data.frame(df)) stop(sprintf("[ERPM_PLE] nodes[[%d]] must be a data.frame.", t))
    if (nrow(df) != length(partitions[[t]])) {
      stop(sprintf("[ERPM_PLE] nodes[[%d]] has %d rows but partition has %d actors.",
                   t, nrow(df), length(partitions[[t]])))
    }
    if (!("label" %in% colnames(df))) {
      stop(sprintf("[ERPM_PLE] nodes[[%d]] must contain a 'label' column.", t))
    }

    # Enforce schema consistency (same columns in same order)
    if (is.null(cols_ref)) {
      cols_ref <- colnames(df)
    } else if (!identical(colnames(df), cols_ref)) {
      stop(sprintf("[ERPM_PLE] nodes[[%d]] columns differ from nodes[[%d]].", t, idx_est[1L]))
    }

    # Preserve original labels as a monadic covariate
    df2 <- df
    names(df2)[names(df2) == "label"] <- "label_raw"

    out <- if (is.null(out)) df2 else rbind(out, df2)
  }

  meta_nA <- nrow(out)

  # Builder requires unique, non-empty labels
  out$label <- paste0("A", seq_len(meta_nA))

  # Put label first
  out <- out[, c("label", setdiff(colnames(out), "label")), drop = FALSE]
  out
}

#' Build meta dyads (block-diagonal over idx_est) (internal helper)
#' @noRd
.erpm_ple_make_meta_dyads <- function(dyads_mode, dyads_data, partitions, idx_est) {
  if (is.null(dyads_data)) return(NULL)

  nA_by_t <- vapply(partitions, length, integer(1))
  meta_nA <- sum(nA_by_t[idx_est])

  if (identical(dyads_mode, "timeline")) {
    out <- list()
    for (nm in names(dyads_data)) {
      mats <- dyads_data[[nm]]

      # Dimension checks per selected time
      for (t in idx_est) {
        M <- mats[[t]]
        nA <- nA_by_t[t]
        if (!is.matrix(M) || nrow(M) != nA || ncol(M) != nA) {
          stop(sprintf("[ERPM_PLE] dyads '%s' at t=%d must be %dx%d.", nm, t, nA, nA))
        }
      }

      out[[nm]] <- .erpm_ple_blockdiag_mats(mats, idx = idx_est, label = nm)
    }
    return(out)
  }

  if (identical(dyads_mode, "meta")) {
    out <- list()
    for (nm in names(dyads_data)) {
      M <- dyads_data[[nm]]
      if (!is.matrix(M) || nrow(M) != meta_nA || ncol(M) != meta_nA) {
        stop(sprintf("[ERPM_PLE] meta dyads '%s' must be %dx%d (meta_nA=%d).",
                     nm, meta_nA, meta_nA, meta_nA))
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
#' Aggregates partitions/nodes/dyads across selected blocks (idx_est) and builds a
#' bipartite meta-network using \code{build_bipartite_from_inputs()}.
#'
#' @param partitions List of partitions (timeline).
#' @param idx_est Integer vector of selected time indices included in meta-network.
#' @param nodes NULL or list length T of per-time node data.frames.
#' @param dyads NULL or dyadic timeline input (current behavior preserved).
#' @param group_labels Group labels (currently not used in meta build; preserved for API).
#' @param verbose Verbosity.
#'
#' @return A list containing \code{meta_nw} plus derived objects used downstream.
#'
#' @noRd
.erpm_long_empile_build_standard_meta_network <- function(partitions,
                                                         idx_est,
                                                         nodes,
                                                         dyads,
                                                         group_labels,
                                                         verbose) {
  # ---------------------------------------------------------------------------
  # NOTE: current dyads normalization behavior is kept as-is (structure-only).
  # ---------------------------------------------------------------------------
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
                    T, paste(idx_est, collapse = ",")))
  }

  # --- Meta partition (concatenate + shift) ----------------------------------
  mp <- .erpm_ple_make_meta_partition(partitions, idx_est = idx_est)
  meta_partition <- mp$meta_partition
  nA_by_t        <- mp$nA_by_t
  offsets_by_t   <- mp$offsets_by_t
  nA_est         <- mp$nA_est
  offsets_est    <- mp$offsets_est
  meta_nA        <- sum(nA_est)

  # --- Meta nodes / dyads -----------------------------------------------------
  meta_nodes_df <- .erpm_ple_make_meta_nodes_df(nodes_n, partitions, idx_est = idx_est)
  meta_dyads    <- .erpm_ple_make_meta_dyads(dyads_mode, dyads_data, partitions, idx_est = idx_est)

  # --- Canonical bipartite build ---------------------------------------------
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
    meta_nA        = meta_nA,
    nA_by_t        = nA_by_t,
    offsets_by_t   = offsets_by_t,
    nA_est         = nA_est,
    offsets_est    = offsets_est,
    meta_nodes_df  = meta_nodes_df,
    meta_dyads     = meta_dyads,
    dyads_mode     = dyads_mode,
    nodes_n        = nodes_n,
    dyads_data     = dyads_data
  )
}

# ==============================================================================
# Module 2: build and attach timeblock
# ==============================================================================

#' Build timeblock for PLE meta-networks (internal helper)
#'
#' @description
#' With \code{build_bipartite_from_inputs()}, the meta-network has \code{2*meta_nA}
#' vertices: actors (1..meta_nA) and padded groups (meta_nA+1 .. 2*meta_nA).
#' This helper assigns a time index to each actor and each padded group.
#'
#' @param idx_est Selected time indices included in meta-network.
#' @param nA_by_t Actor counts per time.
#'
#' @return Integer vector of length \code{2*meta_nA}.
#'
#' @noRd
.erpm_ple_make_timeblock <- function(idx_est, nA_by_t) {
  nA_est  <- nA_by_t[idx_est]
  meta_nA <- sum(nA_est)

  tb_actor <- integer(meta_nA)
  tb_group <- integer(meta_nA)

  a0 <- 0L
  g0 <- 0L

  for (b in seq_along(idx_est)) {
    t  <- idx_est[b]
    nA <- nA_est[b]

    tb_actor[(a0 + 1L):(a0 + nA)] <- t
    tb_group[(g0 + 1L):(g0 + nA)] <- t

    a0 <- a0 + nA
    g0 <- g0 + nA
  }

  c(tb_actor, tb_group)
}

#' Attach timeblock to meta-network (PLE module 2)
#' @noRd
.erpm_long_empile_attach_timeblock <- function(meta_nw, idx_est, nA_by_t) {
  timeblock <- .erpm_ple_make_timeblock(idx_est, nA_by_t = nA_by_t)
  network::set.vertex.attribute(meta_nw, "timeblock", timeblock)
  meta_nw
}

# ==============================================================================
# Module 3: attach inertial attributes
# ==============================================================================

#' Attach inertial attributes required by inertial InitErgmTerms (PLE module 3)
#'
#' @description
#' Current behavior is tailored to \code{inertia_groups} and attaches the exact
#' network attributes expected by your InitErgmTerm implementation.
#'
#' @noRd
.erpm_long_empile_attach_inertial_attributes <- function(meta_nw,
                                                        partitions,
                                                        idx_est,
                                                        d,
                                                        verbose) {
  if (d < 1L) stop("[ERPM_PLE] inertia_groups requires past_influence >= 1.")

  n_block <- .erpm_ple_assert_constant_n_block(partitions, idx_est = idx_est)

  B <- length(idx_est)
  G_block <- n_block

  # Sanity check: bipartite size must be n_block * B
  n1_total <- as.integer(network::get.network.attribute(meta_nw, "bipartite"))
  if (is.na(n1_total) || n1_total != n_block * B) {
    stop(sprintf(
      "[ERPM_PLE] inconsistent bipartite size for inertia_groups: bipartite=%s but n_block*B=%d*%d=%d.",
      as.character(n1_total), n_block, B, n_block * B
    ))
  }

  erpm_block_past_partitions <- .erpm_ple_make_erpm_block_past_partitions(
    partitions = partitions,
    idx_est    = idx_est,
    d          = d,
    n_block    = n_block
  )

  # Attributes expected by InitErgmTerm.inertia_groups.R (PLE mode)
  network::set.network.attribute(meta_nw, "erpm_mode", "empile")
  network::set.network.attribute(meta_nw, "erpm_B", B)
  network::set.network.attribute(meta_nw, "erpm_n", n_block)
  network::set.network.attribute(meta_nw, "erpm_G", G_block)
  network::set.network.attribute(meta_nw, "erpm_block_past_partitions", erpm_block_past_partitions)

  if (isTRUE(verbose)) {
    message(sprintf("[ERPM_PLE] inertia_groups attrs attached: erpm_mode=empile | erpm_B=%d | erpm_n=%d | erpm_G=%d",
                    B, n_block, G_block))
  }

  meta_nw
}

# ==============================================================================
# Orchestrator: called by erpm_long()
# ==============================================================================

#' Build PLE meta-network and attach required attributes (engine entry point)
#'
#' @description
#' Orchestrates the three PLE modules:
#' \enumerate{
#'   \item build standard bipartite meta-network;
#'   \item attach timeblock;
#'   \item attach inertial attributes if needed.
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
#' @param verbose Verbosity.
#'
#' @return A list with at least:
#' \itemize{
#'   \item \code{meta_nw}: the built meta-network
#'   \item \code{idx_est}: indices of blocks included in meta-network
#'   \item \code{rhs}: passthrough RHS
#' }
#'
#' @noRd
.erpm_long_empile_build_meta_nw <- function(partitions,
                                           rhs,
                                           inertial_present = FALSE,
                                           past_influence = 0L,
                                           nodes = NULL,
                                           dyads = NULL,
                                           group_labels = NULL,
                                           directed = FALSE,
                                           verbose = FALSE) {
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
  idx_est <- if (inertial_present) seq.int(d + 1L, T) else seq_len(T)

  # ---------------------------------------------------------------------------
  # 1) Standard meta-network build
  # ---------------------------------------------------------------------------
  b <- .erpm_long_empile_build_standard_meta_network(
    partitions   = partitions,
    idx_est      = idx_est,
    nodes        = nodes,
    dyads        = dyads,
    group_labels = group_labels,
    verbose      = verbose
  )
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
    idx_est = idx_est,
    nA_by_t = b$nA_by_t
  )

  # ---------------------------------------------------------------------------
  # 3) Attach inertial attributes (if requested)
  # ---------------------------------------------------------------------------
  if (inertial_present) {
    meta_nw <- .erpm_long_empile_attach_inertial_attributes(
      meta_nw    = meta_nw,
      partitions = partitions,
      idx_est    = idx_est,
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
  network::set.network.attribute(meta_nw, "erpm_long.idx_est", idx_est)

  nA_est <- b$nA_by_t[idx_est]
  actor_offsets <- c(0L, cumsum(nA_est))[seq_along(nA_est)]

  network::set.network.attribute(meta_nw, "erpm_long.meta_nA", b$meta_nA)
  network::set.network.attribute(meta_nw, "erpm_long.nA_by_t", b$nA_by_t)
  network::set.network.attribute(meta_nw, "erpm_long.offsets_by_t", b$offsets_by_t)
  network::set.network.attribute(meta_nw, "erpm_long.actor_offsets", actor_offsets)
  network::set.network.attribute(meta_nw, "erpm_long.group_id_offsets_by_t", b$offsets_by_t)

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
    idx_est      = idx_est,
    timeline_nws = if (inertial_present) network::get.network.attribute(meta_nw, "erpm_long.timeline_nws") else NULL,
    rhs          = rhs
  )
}