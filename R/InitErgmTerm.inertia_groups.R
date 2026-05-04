################################################################################
# FILE: R/InitErgmTerm.inertia_groups.R
################################################################################
#' ERGM inertial term: inertia_groups (PLE)
#'
#' @name InitErgmTerm.inertia_groups
#' @aliases inertia_groups
#' @note InitErgmTerm.inertia_groups.R
#'
#' @description
#' \code{inertia_groups} is a longitudinal (inertial) ERGM term intended to be
#' used with \code{erpm_long()} in PLE (stacked) mode: one block-diagonal
#' bipartite meta-network, with per-block past observed partitions attached.
#'
#' This InitErgmTerm is responsible for:
#' \enumerate{
#'   \item validating user arguments (past_influence, type, size, ...);
#'   \item deriving block structure from engine-set attributes
#'         (\code{erpm_long.selected_partition_indices},
#'          \code{erpm_long.nbr_actors_by_t}, vertex attribute \code{timeblock});
#'   \item computing per-block effective comparison size \code{n_eff} and
#'         warning if partition lengths differ across time;
#'   \item packing \code{inputs} for the C changestat \code{c_inertia_groups}.
#' }
#'
#' IMPORTANT:
#' The changestat code is paradigm-agnostic for \code{type="exogenous"}.
#' All paradigm-specific work is done here by preparing \code{inputs}.
#'
#' Debugging:
#'   options(ERPM.inertia_groups.debug = TRUE) to enable debug logs
#'   options(ERPM.inertia_groups.debug = "deep") for verbose debug logs
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

#' Debug printer (local)
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

#' Partition -> list of groups (local actor indices 1..n_eff), as integer vectors
#' @noRd
.inertia_groups_groups_from_partition <- function(p) {
  p <- as.integer(p)
  split(seq_along(p), p)
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

  a <- ergm::check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("past_influence", "type", "size", "debug"),
    vartypes      = c("numeric", "character", "ANY", "ANY"),
    defaultvalues = list(1, "exogenous", NULL, NULL),
    required      = c(FALSE, FALSE, FALSE, FALSE)
  )

  # ---------------------------------------------------------------------------
  # Debug control
  # ---------------------------------------------------------------------------
  opt_dbg    <- getOption("ERPM.inertia_groups.debug", TRUE)
  debug_raw  <- a$debug
  if (is.null(debug_raw)) debug_raw <- opt_dbg

  debug <- isTRUE(debug_raw) ||
    (is.character(debug_raw) && length(debug_raw) == 1L && !is.na(debug_raw) &&
       tolower(trimws(debug_raw)) %in% c("deep", "true", "t", "1"))
  deep  <- is.character(debug_raw) && length(debug_raw) == 1L && !is.na(debug_raw) &&
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
  # 1) Read bipartite size and validate erpm_mode (PLE only)
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
  if (is.null(erpm_mode) || is.na(erpm_mode) ||
      !identical(as.character(erpm_mode), "empile")) {
    .inertia_groups_stop(termname,
      "inertia_groups requires a PLE (empile) meta-network. ",
      "Missing or invalid %%n%% 'erpm_mode' (got: ",
      if (is.null(erpm_mode)) "NULL" else as.character(erpm_mode), ").")
  }

  if (identical(type, "endogenous")) {
    .inertia_groups_stop(termname, "`type=\"endogenous\"` is not implemented yet.")
  }
  if (!identical(type, "exogenous")) {
    .inertia_groups_stop(termname,
      "`type` must be \"exogenous\" (default) or \"endogenous\" (not implemented).")
  }

  # ---------------------------------------------------------------------------
  # 2) Derive block structure from engine-set bookkeeping attributes
  # ---------------------------------------------------------------------------
  sel_idx  <- network::get.network.attribute(nw, "erpm_long.selected_partition_indices")
  nbr_by_t <- network::get.network.attribute(nw, "erpm_long.nbr_actors_by_t")

  if (is.null(sel_idx) || !length(sel_idx)) {
    .inertia_groups_stop(termname,
      "missing 'erpm_long.selected_partition_indices' on the network. ",
      "Was the meta-network built with erpm_long()?")
  }
  if (is.null(nbr_by_t) || !length(nbr_by_t)) {
    .inertia_groups_stop(termname,
      "missing 'erpm_long.nbr_actors_by_t' on the network.")
  }

  sel_idx  <- as.integer(sel_idx)
  nbr_by_t <- as.integer(nbr_by_t)
  B        <- length(sel_idx)
  n_b      <- nbr_by_t[sel_idx]  # actor count per selected block, length B

  if (sum(n_b) != n1_total) {
    .inertia_groups_stop(termname, sprintf(
      "inconsistent actor count: bipartite=%d but sum(n_b over B=%d blocks)=%d (n_b=[%s]).",
      n1_total, B, sum(n_b), paste(n_b, collapse = ",")))
  }

  # 0-based cumulative actor offsets: actor_offsets[b] = sum(n_b[1..b-1])
  actor_offsets <- as.integer(c(0L, cumsum(n_b))[seq_len(B)])

  # Build group_to_block[1..n1_total] from the timeblock vertex attribute.
  # timeblock labels use ORIGINAL time indices (from selected_partition_indices),
  # so we invert: time index t -> block index b (1-indexed).
  N  <- network::network.size(nw)   # = 2 * n1_total
  tb <- network::get.vertex.attribute(nw, "timeblock")
  if (is.null(tb) || length(tb) != N) {
    .inertia_groups_stop(termname,
      sprintf("vertex attribute 'timeblock' missing or wrong length (got %s, expected %d).",
              if (is.null(tb)) "NULL" else length(tb), N))
  }
  tb_groups <- as.integer(tb[(n1_total + 1L):N])   # timeblock of group vertices

  max_t     <- max(sel_idx)
  time_to_b <- integer(max_t)
  time_to_b[sel_idx] <- seq_len(B)
  group_to_block <- time_to_b[tb_groups]   # length n1_total, values 1..B

  if (anyNA(group_to_block) || any(group_to_block < 1L) || any(group_to_block > B)) {
    .inertia_groups_stop(termname,
      "group_to_block mapping failed: some timeblock labels in group vertices ",
      "do not correspond to any selected_partition_indices.")
  }

  .inertia_groups_dbg(termname, debug,
                      "blocks: B=%d | n_b=[%s] | actor_offsets=[%s] | d=%d | L=%d",
                      B, paste(n_b, collapse = ","),
                      paste(actor_offsets, collapse = ","), d, L)

  # ---------------------------------------------------------------------------
  # 3) Extract past partitions and compute n_eff per block
  # ---------------------------------------------------------------------------
  past_by_block <- network::get.network.attribute(nw, "erpm_block_past_partitions")
  if (is.null(past_by_block)) {
    past_by_block <- network::get.network.attribute(nw, "erpm_past_partitions")
  }

  if (is.null(past_by_block) || !is.list(past_by_block) || length(past_by_block) != B) {
    .inertia_groups_stop(termname, sprintf(
      "expected %%n%% 'erpm_block_past_partitions' as a list of length B=%d.", B))
  }

  for (b in seq_len(B)) {
    pb <- past_by_block[[b]]
    if (is.null(pb) || !is.list(pb) || length(pb) < d) {
      .inertia_groups_stop(termname,
        sprintf("block %d has insufficient past partitions: need d=%d lags.", b, d))
    }
  }

  # n_eff[b] = min of the current block size and all past partition sizes for b.
  # Actors are identified by position; if sizes differ we only compare the first
  # n_eff[b] actors (those present at all time points considered for block b).
  n_eff <- integer(B)
  for (b in seq_len(B)) {
    past_sizes <- vapply(past_by_block[[b]][seq_len(d)], length, integer(1))
    n_eff[b]   <- min(c(n_b[b], past_sizes))
  }

  all_sizes <- c(n_b, unlist(lapply(seq_len(B), function(b)
    vapply(past_by_block[[b]][seq_len(d)], length, integer(1))
  )))
  if (length(unique(all_sizes)) > 1L) {
    warning(sprintf(
      "[inertia_groups] partitions have different sizes (%s); only the first %d actor(s) (minimum partition size) will be considered for inertia matching.",
      paste(sort(unique(all_sizes)), collapse = ","),
      min(all_sizes)
    ))
  }

  if (isTRUE(deep)) {
    .inertia_groups_dbg(termname, debug,
                        "past(deep): n_eff=[%s]", paste(n_eff, collapse = ","))
  }

  # ---------------------------------------------------------------------------
  # 4) Pack inputs for the C changestat
  # ---------------------------------------------------------------------------
  #
  # INPUT_PARAM layout (0-based indices in C / doubles):
  #
  #   [0]                          n1_total
  #   [1]                          B
  #   [2]                          d
  #   [3]                          L
  #   [4 .. 3+L]                   sizes_int[0..L-1]
  #   [4+L .. 3+L+B]               actor_offsets[0..B-1]   (0-based)
  #   [4+L+B .. 3+L+2B]            n_eff[0..B-1]
  #   [4+L+2B .. 3+L+2B+n1_total]  group_to_block[0..n1_total-1]  (1..B)
  #   [4+L+2B+n1_total .. 3+L+2B+n1_total+B*d]  offsets[0..B*d-1]  (0-based)
  #   data blocks (block-major, lag-major):
  #     for each (b, lag): M, len_1, ids_1..., len_2, ids_2..., ...
  #     (ids are GLOBAL actor ids 1..n1_total, already truncated to n_eff[b])

  inputs <- c(
    as.numeric(n1_total),
    as.numeric(B),
    as.numeric(d),
    as.numeric(L),
    as.numeric(sizes_int),
    as.numeric(actor_offsets),
    as.numeric(n_eff),
    as.numeric(group_to_block)
  )

  offsets_start <- length(inputs) + 1L          # 1-based R index of first offset
  inputs        <- c(inputs, rep(0, B * d))      # placeholder offsets (filled later)

  offsets <- integer(B * d)

  for (b in seq_len(B)) {
    pb    <- past_by_block[[b]]
    off_b <- actor_offsets[b]   # 0-based global actor offset for block b

    for (lag in seq_len(d)) {
      start0 <- length(inputs)              # 0-based C index of this data block
      offsets[(b - 1L) * d + lag] <- start0

      p_lag <- pb[[lag]]
      p_eff <- p_lag[seq_len(n_eff[b])]    # truncate to n_eff[b] actors

      if (!is.atomic(p_eff)) {
        .inertia_groups_stop(termname,
          sprintf("invalid past partition at block=%d lag=%d.", b, lag))
      }

      groups        <- .inertia_groups_groups_from_partition(p_eff)
      groups        <- .inertia_groups_filter_by_size(groups, sizes_int)
      # Global actor ids: off_b (0-based offset) + local id (1-based)
      groups_global <- lapply(groups, function(v) sort(off_b + as.integer(v)))

      M         <- length(groups_global)
      block_vec <- as.numeric(M)
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

  # Fill offsets table
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
  size_tag  <- if (!L) "all" else paste0("size=", paste(sizes_int, collapse = ","))
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
