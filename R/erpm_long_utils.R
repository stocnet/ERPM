################################################################################
# FILE: R/erpm_long_utils.R
################################################################################
#' ERPM long utilities: shared helpers for PLE meta-networks and inertial terms
#'
#' @name erpm_long_utils
#' @note erpm_long_utils.R
#'
#' @description
#' This file groups small utilities used by \code{erpm_long()} in PLE mode and by
#' inertial InitErgmTerms that consume PLE meta-networks.
#'
#' The helpers are intentionally lightweight and focus on:
#' \itemize{
#'   \item safe access to meta-network attributes written by the PLE engine;
#'   \item standardized verbosity/debug messaging for long workflows;
#'   \item convenience getters for inertia-specific PLE attributes
#'         (currently \code{inertia_groups}).
#' }
#'
#' There is no PLS logic here by design: this module only supports the stacked
#' meta-network pathway.
#'
#' @keywords ERPM ERGM longitudinal PLE utils
################################################################################

# ------------------------------------------------------------------------------
# Get timeline networks from a meta-network (PLE)
# ------------------------------------------------------------------------------

#' Fetch timeline networks attached to a PLE meta-network (internal helper)
#' @param nw PLE meta-network.
#' @return List of timeline networks.
#' @noRd
.erpm_long_get_timeline_nws <- function(nw) {
  tl <- network::get.network.attribute(nw, "erpm_long.timeline_nws")
  if (is.null(tl)) {
    stop("[ERPM_LONG] No timeline networks attached to meta-network.")
  }
  if (!is.list(tl) || !length(tl)) {
    stop("[ERPM_LONG] Invalid 'erpm_long.timeline_nws' attribute.")
  }
  tl
}

# ------------------------------------------------------------------------------
# Get estimation block indices
# ------------------------------------------------------------------------------

#' Fetch selected partition indices from a PLE meta-network (internal helper)
#' @param nw PLE meta-network.
#' @return Integer vector of selected time indices.
#' @noRd
.erpm_long_get_idx_est <- function(nw) {
  idx <- network::get.network.attribute(nw, "erpm_long.selected_partition_indices")
  if (is.null(idx)) {
    stop("[ERPM_LONG] Missing 'erpm_long.selected_partition_indices' attribute on meta-network.")
  }
  as.integer(idx)
}

# ------------------------------------------------------------------------------
# Get past_influence (d)
# ------------------------------------------------------------------------------

#' Fetch past influence depth (d) from a PLE meta-network (internal helper)
#' @param nw PLE meta-network.
#' @return Integer scalar past influence.
#' @noRd
.erpm_long_get_past_influence <- function(nw) {
  d <- network::get.network.attribute(nw, "erpm_long.d")
  if (is.null(d)) {
    stop("[ERPM_LONG] Missing 'erpm_long.d' attribute on meta-network.")
  }
  as.integer(d)
}

# ------------------------------------------------------------------------------
# Get total number of partitions T
# ------------------------------------------------------------------------------

#' Fetch total timeline length (T) from a PLE meta-network (internal helper)
#' @param nw PLE meta-network.
#' @return Integer scalar T.
#' @noRd
.erpm_long_get_T <- function(nw) {
  T <- network::get.network.attribute(nw, "erpm_long.T")
  if (is.null(T)) {
    stop("[ERPM_LONG] Missing 'erpm_long.T' attribute on meta-network.")
  }
  as.integer(T)
}

# ------------------------------------------------------------------------------
# Get block structure (sizes + offsets)
# ------------------------------------------------------------------------------

#' Fetch block sizes/offsets metadata (internal helper)
#' @param nw PLE meta-network.
#' @return List with integer vectors: sizes and offsets.
#' @noRd
.erpm_long_get_blocks <- function(nw) {
  sizes   <- network::get.network.attribute(nw, "erpm_block_sizes")
  offsets <- network::get.network.attribute(nw, "erpm_block_offsets")

  if (is.null(sizes) || is.null(offsets)) {
    stop("[ERPM_LONG] Meta-network is missing block metadata.")
  }

  list(
    sizes   = as.integer(sizes),
    offsets = as.integer(offsets)
  )
}

# ------------------------------------------------------------------------------
# Logging helpers (verbose=user, debug=dev)
# ------------------------------------------------------------------------------

#' Verbose console printer (internal helper)
#' @param verbose Logical.
#' @param ... Passed to message().
#' @noRd
.erpm_long_vcat <- function(verbose, ...) {
  if (isTRUE(verbose)) message(...)
  invisible(NULL)
}

#' Debug console printer (internal helper)
#' @param debug Logical.
#' @param ... Passed to message().
#' @noRd
.erpm_long_dcat <- function(debug, ...) {
  if (isTRUE(debug)) message(...)
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# Check whether meta-network was built with inertial support
# ------------------------------------------------------------------------------

#' Quick predicate: does the meta-network have inertial plumbing? (internal helper)
#' @param nw PLE meta-network.
#' @return Logical.
#' @noRd
.erpm_long_has_inertia <- function(nw) {
  !is.null(network::get.network.attribute(nw, "erpm_long.timeline_nws"))
}

# ------------------------------------------------------------------------------
# Convenience: fetch meta-node attribute (monadic covariate)
# ------------------------------------------------------------------------------

#' Fetch a meta-network attribute used as a monadic covariate (internal helper)
#' @param nw PLE meta-network.
#' @param name Attribute name.
#' @return The attribute value.
#' @noRd
.erpm_long_get_meta_node <- function(nw, name) {
  x <- network::get.network.attribute(nw, name)
  if (is.null(x)) {
    stop(sprintf("[ERPM_LONG] Meta-node attribute '%s' not found.", name))
  }
  x
}

# ------------------------------------------------------------------------------
# Convenience: fetch meta-dyad attribute (dyadic covariate)
# ------------------------------------------------------------------------------

#' Fetch a meta-network attribute used as a dyadic covariate (internal helper)
#' @param nw PLE meta-network.
#' @param name Attribute name.
#' @return The attribute value.
#' @noRd
.erpm_long_get_meta_dyad <- function(nw, name) {
  x <- network::get.network.attribute(nw, name)
  if (is.null(x)) {
    stop(sprintf("[ERPM_LONG] Meta-dyad attribute '%s' not found.", name))
  }
  x
}

# ------------------------------------------------------------------------------
# Global options (namespaced)
# ------------------------------------------------------------------------------

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.erpm_long <- list(
    erpm.long.verbose = FALSE
  )
  toset <- !(names(op.erpm_long) %in% names(op))
  if (any(toset)) options(op.erpm_long[toset])
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# Inertia-specific getters (PLE)
# ------------------------------------------------------------------------------

#' Fetch PLE inertia_groups block metadata from a meta-network (internal helper)
#' @param nw PLE meta-network.
#' @return List with B, n_block, G_block, past_by_block.
#' @noRd
.erpm_long_get_erpm_blocks_for_inertia <- function(nw) {
  mode <- network::get.network.attribute(nw, "erpm_mode")
  if (is.null(mode) || is.na(mode) || !identical(as.character(mode), "empile")) {
    stop("[ERPM_LONG] inertia_groups expects PLE meta-network with %n% 'erpm_mode' == 'empile'.")
  }

  B <- as.integer(network::get.network.attribute(nw, "erpm_B"))
  n <- as.integer(network::get.network.attribute(nw, "erpm_n"))
  G <- as.integer(network::get.network.attribute(nw, "erpm_G"))

  if (any(is.na(c(B, n, G))) || any(c(B, n, G) <= 0L)) {
    stop("[ERPM_LONG] missing/invalid inertia_groups attrs: erpm_B, erpm_n, erpm_G.")
  }

  past <- network::get.network.attribute(nw, "erpm_block_past_partitions")
  if (is.null(past) || !is.list(past) || length(past) != B) {
    stop(sprintf("[ERPM_LONG] missing/invalid erpm_block_past_partitions (expected list length B=%d).", B))
  }

  list(B = B, n_block = n, G_block = G, past_by_block = past)
}