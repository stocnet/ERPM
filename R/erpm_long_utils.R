################################################################################
# FILE: R/erpm_long_utils.R
# OBJECT: Utilities shared by PLE-only erpm_long() and inertial InitErgmTerms
# NOTES :
#   - No PLS logic.
#   - Pure helpers: access to meta-network attributes built by the engine.
#   - Used by InitErgmTerm.inertia_groups (and future inertial terms).
################################################################################

# ------------------------------------------------------------------------------
# Get timeline networks from a meta-network (PLE)
# ------------------------------------------------------------------------------
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
.erpm_long_vcat <- function(verbose, ...) {
  if (isTRUE(verbose)) message(...)
  invisible(NULL)
}

.erpm_long_dcat <- function(debug, ...) {
  if (isTRUE(debug)) message(...)
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# Check whether meta-network was built with inertial support
# ------------------------------------------------------------------------------
.erpm_long_has_inertia <- function(nw) {
  !is.null(network::get.network.attribute(nw, "erpm_long.timeline_nws"))
}

# ------------------------------------------------------------------------------
# Convenience: fetch meta-node attribute (monadic covariate)
# ------------------------------------------------------------------------------
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
# Utility: read verbose option
# ------------------------------------------------------------------------------
# .erpm_long_opt_verbose <- function() {
#   isTRUE(getOption("erpm.long.verbose", FALSE))
# }

# ==============================================================================
# Add getters for inertia_groups-specific PLE attributes.
# ==============================================================================
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