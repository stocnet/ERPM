################################################################################
# FILE: R/erpm_long_inertia_registry.R
################################################################################
#' ERPM longitudinal inertia registry: declarative term specs for erpm_long()
#' @name erpm_long_inertia_registry
#' @note erpm_long_inertia_registry.R
#'
#' @description
#' This file defines a registry of *inertial* terms used by \code{erpm_long()}.
#' An inertial term summarizes information from *past* partitions/networks and
#' attaches it as per-lag attributes that later change-statistics can use.
#'
#' Concretely, each registry entry describes:
#' \enumerate{
#'   \item how many past time steps to look back (\code{past_influence});
#'   \item what to compute or attach for each lag (\code{build_attr});
#'   \item optional metadata (\code{label}, \code{description}, \code{args},
#'         \code{attr_prefix}) used for readability and debugging.
#' }
#'
#' The registry is purely declarative: it does not fit models. It only tells
#' \code{erpm_long()} how to build the lagged attributes needed by inertia-aware
#' ERGM terms.
#'
#' @keywords ERPM ERGM longitudinal inertia registry
NULL

#' @description Internal registry used by erpm_long(). See file-level docs.
#' @return A named list. Each element is a term specification with callbacks.
#' @keywords internal
#' @noRd
.erpm_long_inertia_registry <- list(

  # --------------------------------------------------------------------------
  # inertia_groups(size = NULL, past_influence = 1)
  #
  # Intended interpretation (user):
  # - A group at time t is persistent if the exact same set of actors appears
  #   as a group in at least one past partition within the past_influence window.
  #
  # Attachment (per lag):
  # - A compact representation of groups at t-lag: signatures + sizes.
  # - size filter is stored to mirror the call intent (used later by changestat).
  # --------------------------------------------------------------------------
  inertia_groups = list(
    label       = "Inertia of exact groups (set persistence)",
    description = "Persistence of groups identified by exact actor sets within a past window.",
    args        = c("size", "past_influence"),
    attr_prefix = "erpm_inertia__inertia_groups",

    # Read past_influence from the user call. Default = 1 means “look back one step”.
    past_influence = function(call, env0) .erpm_long_get_past_influence(call, env0, default = 1L),

    # Build the per-lag attribute object that will be attached to nw_t.
    build_attr = function(nw_t, nets, parts, t, lag, call, env0, debug = FALSE) {
      # Prefer partitions directly (erpm_long input). Fallback to net extraction.
      # This makes the registry robust to two common inputs:
      # - parts provided explicitly
      # - only networks provided, so we reconstruct partitions from them
      p_prev <- NULL
      idx <- t - lag
      if (is.list(parts) && length(parts) >= idx && idx >= 1L) {
        # Round + integer coercion standardizes partitions to integer labels.
        # This avoids silent issues if upstream produced numeric-but-integer-ish values.
        p_prev <- as.integer(round(parts[[idx]]))
        .erpm_long_dbg(debug, "[ERPM_LONG|DEBUG] inertia_groups: using parts[[", idx, "]]")
      } else {
        # If partitions are not available, extract them from the past network.
        nw_prev <- nets[[idx]]
        p_prev  <- .erpm_long_extract_partition_from_network(nw_prev)
        .erpm_long_dbg(debug, "[ERPM_LONG|DEBUG] inertia_groups: using network extraction at t-lag=", idx)
      }

      # Convert partition labels -> list of groups (each group is a vector of actor indices).
      groups  <- .erpm_long_groups_from_partition(p_prev)

      # Compute stable identifiers for groups (e.g., hashes or canonical signatures).
      # These signatures are the compact key used to test equality across time.
      sig <- .erpm_long_group_signatures(groups)

      # Store group sizes because some inertia terms filter by size downstream.
      sz  <- vapply(groups, length, integer(1))

      # Return a structured attribute object. Later code attaches it to the network.
      list(
        type        = "group_signature_set",
        lag         = as.integer(lag),
        signatures  = as.character(sig),
        sizes       = as.integer(sz),
        # Keep the original size filter intent to ensure changestat uses the same rule.
        size_filter = .erpm_long_get_size_filter(call, env0)
      )
    }
  )

  # Future inertial terms can be added here.
)