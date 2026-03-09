################################################################################
# FILE: R/InitErgmTerm.inertia_groups.R
################################################################################
#' ERGM inertial term: inertia_groups (PLS/PLE-ready)
#'
#' @name InitErgmTerm.inertia_groups
#' @aliases inertia_groups
#' @note InitErgmTerm.inertia_groups.R
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' \code{inertia_groups} is a longitudinal (inertial) ERGM term intended to be used
#' with \code{erpm_long()} under either:
#' \itemize{
#'   \item PLS (sequential): one bipartite network per time, with past partitions attached
#'         to the current network as network attributes;
#'   \item PLE (stacked): one stacked block-diagonal bipartite meta-network, with per-block
#'         past observed partitions attached to the network.
#' }
#'
#' This InitErgmTerm is responsible for:
#' \enumerate{
#'   \item validating user arguments (past_influence, type, size, ...);
#'   \item detecting whether the current network carries PLS or PLE longitudinal attributes;
#'   \item extracting the relevant past partitions from network attributes;
#'   \item packing \code{inputs} for the generic C changestat \code{c_inertia_groups}.
#' }
#'
#' IMPORTANT:
#' The changestat code is paradigm-agnostic for \code{type="exogenous"}.
#' All paradigm-specific work must be done here by preparing \code{inputs}.
#'
#' Debugging:
#'   options(ERPM.inertia_groups.debug = TRUE) to enable debug logs
#'
#' @keywords ERPM ERGM inertial longitudinal
#' @md
NULL

# ------------------------------------------------------------------------------
# Small internal helpers (local to this InitErgmTerm)
# ------------------------------------------------------------------------------

#' Null-coalescing helper (local)
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Debug printer (local, inspired by InitErgmTerm.cliques style)
#'
#' The initializer emits debug logs through ergm's init warning channel:
#' - it avoids polluting the console during MCMC/initialization;
#' - users can inspect messages via warnings().
#'
#' Debug is controlled by an R option (preferred):
#'   options(ERPM.inertia_groups.debug = TRUE) or "deep"
#'
#' A term argument `debug` is still accepted for backward compatibility, but
#' it is treated as an override only when explicitly provided.
#'
#' @noRd
.inertia_groups_dbg <- function(termname, debug, ...) {
  if (!isTRUE(debug)) return(invisible(NULL))
  msg <- sprintf(...)
  if (exists("ergm_Init_warn", mode = "function")) {
    ergm_Init_warn(sQuote(termname), ": ", msg)
  } else {
    cat("[", termname, "|DEBUG] ", msg, "\n", sep = "")
  }
  invisible(NULL)
}

#' Stop helper using ergm_Init_stop when available
#' @noRd
.inertia_groups_stop <- function(termname, ...) {
  msg <- paste0(...)
  if (exists("ergm_Init_stop", mode = "function")) {
    ergm_Init_stop(sQuote(termname), ": ", msg)
  }
  stop(paste0("[", termname, "] ", msg), call. = FALSE)
}

#' Coerce and validate an integer scalar >= 1
#' @noRd
.inertia_groups_as_int1 <- function(x, name = "past_influence", termname = "inertia_groups") {
  if (is.numeric(x) && length(x) == 1L && is.finite(x)) {
    iv <- as.integer(round(x))
    if (isTRUE(all.equal(x, iv)) && iv >= 1L) return(iv)
  }
  .inertia_groups_stop(termname, "`", name, "` must be an integer >= 1.")
}

#' Normalize size filter:
#' - NULL => integer(0) meaning "all sizes"
#' - scalar/vector/interval => set of positive integers, unique, sorted
#' @noRd
.inertia_groups_norm_size <- function(size, termname = "inertia_groups") {
  if (is.null(size)) return(integer(0))

  if (!(is.atomic(size) && length(size) >= 1L)) {
    .inertia_groups_stop(termname, "`size` must be NULL or an atomic vector (e.g., 4, c(1,4), 1:4).")
  }

  v_num <- suppressWarnings(as.numeric(size))
  if (anyNA(v_num) || any(!is.finite(v_num))) {
    .inertia_groups_stop(termname, "`size` contains NA/NaN/Inf.")
  }

  v_int <- as.integer(round(v_num))
  if (!isTRUE(all.equal(v_num, as.numeric(v_int)))) {
    .inertia_groups_stop(termname, "`size` must be integer-valued.")
  }
  if (any(v_int <= 0L)) {
    .inertia_groups_stop(termname, "`size` must contain integers > 0.")
  }

  sort(unique(v_int))
}

#' Partition -> list of groups (actor indices 1..n_block), as integer vectors
#' @noRd
.inertia_groups_groups_from_partition <- function(p) {
  p <- as.integer(p)
  split(seq_along(p), p)
}

#' Groups -> sorted global actor ids for block b (global actor space 1..n1_total)
#' @noRd
.inertia_groups_groups_to_global_ids <- function(groups, b, n_block) {
  off <- (b - 1L) * n_block
  lapply(groups, function(v) sort(off + as.integer(v)))
}

#' Apply size filter on a list of integer vectors
#' @noRd
.inertia_groups_filter_by_size <- function(groups, sizes_int) {
  if (!length(sizes_int)) return(groups)
  keep <- vapply(groups, function(v) length(v) %in% sizes_int, logical(1))
  groups[keep]
}

# ------------------------------------------------------------------------------
# InitErgmTerm
# ------------------------------------------------------------------------------

#' @export
InitErgmTerm.inertia_groups <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "inertia_groups"

  # ---------------------------------------------------------------------------
  # 0) Parse and validate user arguments
  # ---------------------------------------------------------------------------
  # Accept positional inertia_groups(1) as past_influence=1 (only if unambiguous).
  if (length(arglist) == 1L) {
    nm <- names(arglist)
    if (is.null(nm) || isTRUE(nm[1L] == "")) {
      arglist <- list(past_influence = arglist[[1L]])
    }
  }

  # Normalize aliases into canonical names early.
  if (!is.null(names(arglist)) && "pi" %in% names(arglist) && !"past_influence" %in% names(arglist))
    arglist[["past_influence"]] <- arglist[["pi"]]
  if (!is.null(names(arglist)) && "d" %in% names(arglist) && !"past_influence" %in% names(arglist))
    arglist[["past_influence"]] <- arglist[["d"]]
  if (!is.null(names(arglist)) && "sizes" %in% names(arglist) && !"size" %in% names(arglist))
    arglist[["size"]] <- arglist[["sizes"]]

  # Use ergm's standard checking for defaults and basic types.
  # Debug is special: accept TRUE/FALSE/"deep".
  a <- ergm::check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("past_influence", "type", "size", "debug"),
    vartypes      = c("numeric", "character", "ANY", "ANY"),
    defaultvalues = list(1, "exogenous", NULL, NULL), # <- debug default from option below
    required      = c(FALSE, FALSE, FALSE, FALSE)
  )

  # ---------------------------------------------------------------------------
  # Debug control
  # ---------------------------------------------------------------------------
  opt_dbg <- getOption("ERPM.inertia_groups.debug", TRUE)

  debug_raw <- a$debug
  if (is.null(debug_raw)) debug_raw <- opt_dbg

  debug <- isTRUE(debug_raw) ||
    (is.character(debug_raw) && length(debug_raw) == 1L && !is.na(debug_raw) &&
       tolower(trimws(debug_raw)) %in% c("deep", "true", "t", "1"))
  deep <- is.character(debug_raw) && length(debug_raw) == 1L && !is.na(debug_raw) &&
    tolower(trimws(debug_raw)) == "deep"

  d <- .inertia_groups_as_int1(a$past_influence, name = "past_influence", termname = termname)

  type <- a$type %||% "exogenous"
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    .inertia_groups_stop(termname, "`type` must be a single string: \"exogenous\" or \"endogenous\".")
  }
  type <- tolower(trimws(type))

  sizes_int <- .inertia_groups_norm_size(a$size, termname = termname)
  L <- length(sizes_int)

  .inertia_groups_dbg(termname, debug,
                      "args: type=%s | past_influence=%d | size_filter=%s | deep=%s",
                      type, d,
                      if (!L) "all" else paste(sizes_int, collapse = ","),
                      if (isTRUE(deep)) "TRUE" else "FALSE")

  # ---------------------------------------------------------------------------
  # 1) Read bipartite size and longitudinal mode
  # ---------------------------------------------------------------------------
  n1_total <- network::get.network.attribute(nw, "bipartite")
  if (is.null(n1_total) || is.na(n1_total)) {
    .inertia_groups_stop(termname, "non-bipartite network or missing %n% 'bipartite' attribute.")
  }
  n1_total <- as.integer(n1_total)
  if (n1_total <= 0L) {
    .inertia_groups_stop(termname, "invalid bipartite size (N1 <= 0).")
  }

  erpm_mode <- network::get.network.attribute(nw, "erpm_mode")
  if (is.null(erpm_mode) || is.na(erpm_mode)) erpm_mode <- "empile"
  erpm_mode <- as.character(erpm_mode)

  is_PLE <- identical(erpm_mode, "empile")
  is_PLS <- !is_PLE

  .inertia_groups_dbg(termname, debug,
                      "network: erpm_mode=%s | paradigm=%s | n1_total=%d",
                      erpm_mode, if (is_PLE) "PLE" else "PLS", n1_total)

  # type handling
  if ( identical(type, "endogenous") ) {
    if (!is_PLE) {
      .inertia_groups_stop(termname, "`type=\"endogenous\"` is only meaningful in PLE (stacked) mode.")
    }
    .inertia_groups_stop(termname, "`type=\"endogenous\"` is not implemented yet.")
  }
  if (!identical(type, "exogenous")) {
    .inertia_groups_stop(termname, "`type` must be \"exogenous\" (default) or \"endogenous\" (not implemented).")
  }

  # ---------------------------------------------------------------------------
  # 2) Determine block structure (B, n_block, G_block)
  # ---------------------------------------------------------------------------
  if (is_PLE) {
    B       <- as.integer(network::get.network.attribute(nw, "erpm_B"))
    n_block <- as.integer(network::get.network.attribute(nw, "erpm_n"))
    G_block <- as.integer(network::get.network.attribute(nw, "erpm_G"))

    if (any(is.na(c(B, n_block, G_block))) || any(c(B, n_block, G_block) <= 0L)) {
      .inertia_groups_stop(termname, "[PLE] missing/invalid stacked attributes: erpm_B, erpm_n, erpm_G.")
    }
    if (n1_total != n_block * B) {
      .inertia_groups_stop(
        termname,
        sprintf("[PLE] inconsistent actor count: bipartite=%d but erpm_n*erpm_B=%d*%d=%d.",
                n1_total, n_block, B, n_block * B)
      )
    }
  } else {
    # PLS: one block only (current network); keep the same packing layout.
    B       <- 1L
    n_block <- n1_total
    G_block <- n1_total
  }

  .inertia_groups_dbg(termname, debug,
                      "blocks: B=%d | n_block=%d | G_block=%d | d=%d | L=%d",
                      B, n_block, G_block, d, L)

  # ---------------------------------------------------------------------------
  # 3) Extract past partitions from network attributes
  # ---------------------------------------------------------------------------
  if (is_PLE) {
    # Prefer historical attribute name; accept standardized engine name as fallback.
    past_by_block <- network::get.network.attribute(nw, "erpm_block_past_partitions")
    if (is.null(past_by_block)) {
      # Accept standardized engine attribute name
      past_by_block <- network::get.network.attribute(nw, "erpm_past_partitions")
    }

    if (is.null(past_by_block) || !is.list(past_by_block) || length(past_by_block) != B) {
      .inertia_groups_stop(
        termname,
        sprintf("[PLE] expected %%n%% 'erpm_block_past_partitions' or 'erpm_past_partitions' as list(B=%d).", B)
      )
    }

    for (b in seq_len(B)) {
      pb <- past_by_block[[b]]
      if (is.null(pb) || !is.list(pb) || length(pb) < d) {
        .inertia_groups_stop(termname, sprintf("[PLE] block %d has insufficient past partitions: need past_influence=%d lags.", b, d))
      }
    }

    if (isTRUE(deep)) {
      lens <- vapply(past_by_block, length, integer(1))
      .inertia_groups_dbg(termname, debug, "past(deep): per-block lag lengths: %s", paste(lens, collapse = ","))
    }

  } else {
    past_parts <- network::get.network.attribute(nw, "erpm_past_partitions")
    past_depth <- network::get.network.attribute(nw, "erpm_past_depth")

    if (is.null(past_parts) || !is.list(past_parts)) {
      .inertia_groups_stop(termname, "[PLS] expected %n% 'erpm_past_partitions' as a list of lags.")
    }
    if (is.null(past_depth) || is.na(past_depth)) past_depth <- length(past_parts)
    past_depth <- as.integer(past_depth)

    if (past_depth < d || length(past_parts) < d) {
      .inertia_groups_stop(
        termname,
        sprintf("[PLS] insufficient past partitions: need past_influence=%d lags but have past_depth=%d.", d, past_depth)
      )
    }

    # Normalize to a PLE-like container with B=1 for the packing loop.
    past_by_block <- list(past_parts)

    if (isTRUE(deep)) {
      .inertia_groups_dbg(termname, debug,
                          "past(deep): past_depth=%d | length(erpm_past_partitions)=%d",
                          past_depth, length(past_parts))
    }
  }

  # ---------------------------------------------------------------------------
  # 4) Build past group lists (GLOBAL actor ids) per (block, lag)
  # ---------------------------------------------------------------------------
  # INPUT_PARAM layout (numeric vector):
  #   header:
  #     n1_total, n_block, G_block, B, d, L
  #   sizes (L entries):
  #     sizes_int (possibly empty)
  #   offsets table (B*d entries):
  #     offsets[block,lag] = 0-based index in INPUT_PARAM where the (block,lag)
  #     data block begins (i.e., position just before writing M for that block)
  #   data blocks (block-major, lag-major):
  #     for each (b,lag):
  #       M, then for each group m=1..M:
  #         len_m, id_1, ..., id_len_m
  offsets <- integer(B * d)

  inputs <- c(
    as.numeric(n1_total),
    as.numeric(n_block),
    as.numeric(G_block),
    as.numeric(B),
    as.numeric(d),
    as.numeric(L),
    as.numeric(sizes_int)
  )

  offsets_start <- length(inputs) + 1L
  inputs <- c(inputs, rep(0, B * d)) # placeholder offsets (filled later)

  # Deterministic order: block-major, lag-major
  for (b in seq_len(B)) {
    pb <- past_by_block[[b]]

    for (lag in seq_len(d)) {
      # 0-based start position in the final INPUT_PARAM vector
      start0 <- length(inputs)
      offsets[(b - 1L) * d + lag] <- start0

      p_lag <- pb[[lag]]
      if (is.null(p_lag) || !is.atomic(p_lag)) {
        .inertia_groups_stop(termname, sprintf("invalid past partition at block=%d lag=%d (must be an atomic vector).", b, lag))
      }
      if (length(p_lag) != n_block) {
        .inertia_groups_stop(
          termname,
          sprintf("past partition length mismatch at block=%d lag=%d: expected n_block=%d, got %d.",
                  b, lag, n_block, length(p_lag))
        )
      }

      groups <- .inertia_groups_groups_from_partition(p_lag)
      groups <- .inertia_groups_filter_by_size(groups, sizes_int)
      groups_global <- .inertia_groups_groups_to_global_ids(groups, b, n_block)

      M <- length(groups_global)

      block_vec <- numeric(0)
      block_vec <- c(block_vec, as.numeric(M))
      if (M) {
        for (g in groups_global) {
          block_vec <- c(block_vec, as.numeric(length(g)), as.numeric(g))
        }
      }

      inputs <- c(inputs, block_vec)

      .inertia_groups_dbg(termname, debug,
                          "packed block=%d lag=%d: offset0=%d | M=%d",
                          b, lag, start0, M)

      if (isTRUE(deep) && M) {
        sizes_here <- vapply(groups_global, length, integer(1))
        .inertia_groups_dbg(termname, debug,
                            "packed(deep) block=%d lag=%d: group sizes: %s",
                            b, lag, paste(sizes_here, collapse = ","))
      }
    }
  }

  # Fill offsets table (as numeric)
  offsets_end <- offsets_start + (B * d) - 1L
  inputs[offsets_start:offsets_end] <- as.numeric(offsets)

  if (isTRUE(deep)) {
    .inertia_groups_dbg(termname, debug,
                        "inputs(deep): header_len=%d | offsets=[%d..%d] | total_len=%d",
                        offsets_start - 1L, offsets_start, offsets_end, length(inputs))
  }

  # ---------------------------------------------------------------------------
  # 5) Coefficient naming
  # ---------------------------------------------------------------------------
  size_tag <- if (!L) "all" else paste0("size=", paste(sizes_int, collapse = ","))
  coef_name <- sprintf("inertia_groups[type=%s,pi=%d]_%s", type, d, size_tag)

  # ---------------------------------------------------------------------------
  # 6) Return ERGM term spec
  # ---------------------------------------------------------------------------
  list(
    name       = "inertia_groups",
    coef.names = coef_name,
    pkgname    = "ERPM",
    inputs     = inputs,
    dependence = TRUE
  )
}