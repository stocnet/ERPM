################################################################################
# FILE: R/erpm_build_bipartite.R
################################################################################
#' ERPM bipartite builder: build a network from a partition + inputs
#' @name erpm_build_bipartite
#' @note erpm_build_bipartite.R
#'
#' @description
#' This file implements the bipartite builder used by \code{erpm()} when the
#' formula LHS is a partition vector:
#' \enumerate{
#'   \item validate partition, optional node table, and optional dyadic matrices;
#'   \item build a padded bipartite \pkg{network} object (mode-2 padded to \eqn{G=n});
#'   \item attach node attributes and dyadic matrices as network attributes.
#' }
#'
#' @keywords ERPM ERGM  bipartite builder

# ============================================================================
# Bipartite builder and helper
# ============================================================================

#' Stop with a standardized usage message
#' @noRd
.erpm_stop_build_bipartite <- function(msg) {
  usage <- paste(
    "[ERPM] build_bipartite_from_inputs usage:",
    "",
    "  build_bipartite_from_inputs(partition, nodes = NULL, dyads = list(), group_labels = NULL)",
    "",
    "  # Minimal example:",
    "  partition <- c(1, 1, 2, 2, 2, 3)",
    "  nodes <- data.frame(",
    "    label  = c('A','B','C','D','E','F'),",
    "    gender = c(1, 1, 2, 1, 2, 2),",
    "    age    = c(20, 22, 25, 30, 30, 31)",
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

#' Build a padded bipartite network from a partition
#'
#' Construct a bipartite \pkg{network} object from:
#' \itemize{
#'   \item a partition vector of length \eqn{n} (mode-1 actors),
#'   \item optional node attributes in a data.frame,
#'   \item optional dyadic n×n matrices, attached as `%n%` attributes.
#' }
#'
#' The function pads the number of groups to \eqn{G = n} to maximize ERGM
#' variability (unused groups become empty mode-2 vertices).
#'
#' @param partition integer vector of group ids, length n, values in 1..G
#' @param nodes optional data.frame with a label column (auto-detected)
#' @param dyads optional named list of n×n matrices to attach as `%n%` attributes
#' @param group_labels optional character vector of group labels of length G (= n)
#' @return list with components:
#'   \itemize{
#'     \item \code{network}: the bipartite \pkg{network} object;
#'     \item \code{partition}: the original partition vector;
#'     \item \code{actor_labels}: character vector of actor labels;
#'     \item \code{group_labels}: character vector of group labels.
#'   }
#' @examples
#' \dontrun{
#'   partition <- c(1, 1, 2, 2, 3)
#'   nodes <- data.frame(
#'     label = paste0("A", 1:5),
#'     age   = c(20, 21, 22, 23, 24)
#'   )
#'   friendship <- matrix(0, 5, 5)
#'   built <- build_bipartite_from_inputs(partition, nodes = nodes,
#'                                       dyads = list(friendship = friendship))
#'   nw <- built$network
#' }
#'
#' @note
#' This builder is used by \code{erpm()} when the LHS of the formula is
#' a partition vector rather than a pre-built network.
#' The bipartite attribute is set to \code{n} so that \code{constraints = ~ b1part}
#' is meaningful for \pkg{ergm}.
#'
#' @note
#' Self-tests construct bipartite networks from partitions and compare their
#' summaries to reference partitions.
#' @keywords bipartite
#' @export
build_bipartite_from_inputs <- function(partition    = NULL,
                                        nodes        = NULL,
                                        dyads        = list(),
                                        group_labels = NULL) {
  # -- Checks on partition ----------------------------------------------------
  if (is.null(partition) || !is.atomic(partition) || length(partition) < 1L) {
    .erpm_stop_build_bipartite("partition must be a non-empty atomic vector.")
  }
  if (any(!is.finite(partition)) || any(partition < 1)) {
    .erpm_stop_build_bipartite("partition must contain finite positive integers.")
  }

  # Enforce integer-valued group ids early.
  # This avoids silent truncation later when computing group indices.
  ip <- as.integer(round(partition))
  if (!isTRUE(all.equal(partition, ip))) {
    .erpm_stop_build_bipartite(
      "partition must contain integer-valued group ids (e.g., 1,2,3)."
    )
  }
  partition <- ip

  n <- length(partition)

  # -- Nodes / labels ---------------------------------------------------------
  if (is.null(nodes)) {
    labels <- sprintf("A%d", seq_len(n))
    nodes  <- data.frame(label = labels, stringsAsFactors = FALSE)
  } else {
    # Wrap .erpm_check_nodes_df errors in a more informative message
    # Allow nodes as named list of vectors (e.g., list(colors=..., shapes=...))
    if (is.list(nodes) && !is.data.frame(nodes)) {
      if (length(nodes) && (is.null(names(nodes)) || any(!nzchar(names(nodes))))) {
        .erpm_stop_build_bipartite("`nodes` as list must be a named list (e.g., list(colors=..., shapes=...)).")
      }
      if (length(nodes) && !all(vapply(nodes, is.atomic, logical(1)))) {
        .erpm_stop_build_bipartite("`nodes` list must contain only atomic vectors.")
      }
      # Convert to data.frame; label handled below by existing logic
      nodes <- as.data.frame(nodes, stringsAsFactors = FALSE)
    }
    ok_nodes <- try(.erpm_check_nodes_df(nodes), silent = TRUE)
    if (inherits(ok_nodes, "try-error")) {
      .erpm_stop_build_bipartite(
        paste0("Invalid `nodes` data.frame: ", conditionMessage(attr(ok_nodes, "condition")))
      )
    }
    if (nrow(nodes) != n) {
      .erpm_stop_build_bipartite("nrow(nodes) must equal length(partition).")
    }
    labcol <- .erpm_get_label_col(nodes)
    if (is.null(labcol)) {
      labels <- sprintf("A%d", seq_len(n))
      nodes$label <- labels
    } else {
      labels <- as.character(nodes[[labcol]])
      if (anyNA(labels) || any(!nzchar(labels)))
        .erpm_stop_build_bipartite("empty/NA labels are not allowed in `nodes`.")
      if (anyDuplicated(labels))
        .erpm_stop_build_bipartite("duplicate labels are not allowed in `nodes`.")
      if (!("label" %in% names(nodes))) nodes$label <- labels
    }
  }

  # -- Padded groups: G = n ---------------------------------------------------
  G_obs <- max(partition, na.rm = TRUE)
  if (G_obs > n) {
    .erpm_stop_build_bipartite("max(partition) cannot exceed n when padding with G = n.")
  }
  G <- n  # force as many groups as actors
  stopifnot(G >= n)

  # Validate dyadic matrices
  ok_dyads <- try(.erpm_check_dyads(dyads, n, labels), silent = TRUE)
  if (inherits(ok_dyads, "try-error")) {
    .erpm_stop_build_bipartite(
      paste0("Invalid `dyads` list: ", conditionMessage(attr(ok_dyads, "condition")))
    )
  }

  # -- Names and adjacency ----------------------------------------------------
  # Optional group labels allow nicer outputs/logs without changing the model.
  # When provided, we require exactly G labels (G = n here) and uniqueness.
  if (is.null(group_labels)) {
    g_names <- sprintf("G%d", seq_len(G))
  } else {
    if (!is.atomic(group_labels) || length(group_labels) != G) {
      .erpm_stop_build_bipartite("group_labels must be an atomic vector of length G (= n).")
    }
    g_names <- as.character(group_labels)
    if (anyNA(g_names) || any(!nzchar(g_names)))
      .erpm_stop_build_bipartite("empty/NA group_labels are not allowed.")
    if (anyDuplicated(g_names))
      .erpm_stop_build_bipartite("duplicate group_labels are not allowed.")
  }

  all_v <- c(labels, g_names)

  # actor indices: 1..n ; group indices: n + (1..G)
  idx_actor <- seq_len(n)
  idx_group <- n + as.integer(partition)

  # -- Build network object ---------------------------------------------------
  # PERF:
  # We avoid building a dense (n+G)^2 adjacency matrix, which becomes costly
  # when G = n (here: (2n)^2 entries). Instead we create the network from an
  # edgelist (one membership edge per actor).
  edges <- cbind(idx_actor, idx_group)

  nw <- network::network.initialize(n + G, directed = FALSE)
  network::add.edges(nw, tail = edges[, 1L], head = edges[, 2L])

  network::set.network.attribute(nw, "bipartite", n)
  network::set.vertex.attribute(nw, "vertex.names", all_v)

  # Push actor attributes; pad groups with NA
  if (ncol(nodes) > 1L) {
    for (a in setdiff(names(nodes), c("label"))) {
      vals <- nodes[[a]]
      network::set.vertex.attribute(nw, a, c(vals, rep(NA, G)))
    }
  }

  # Attach dyadic n×n matrices as a dedicated network attribute.
  #
  # IMPORTANT:
  # Previously, dyads were attached as top-level network attributes with the same
  # name as the list element (risk of collisions with other attributes).
  # We now store the whole list under a single attribute "dyads".
  # Matrices are still forced to actor label dimnames in a controlled order.
  if (length(dyads)) {
    dyads2 <- dyads
    for (nm in names(dyads2)) {
      M <- dyads2[[nm]]
      dimnames(M) <- list(labels, labels)
      dyads2[[nm]] <- M
    }
    network::set.network.attribute(nw, "dyads", dyads2)
  }

  list(
    network      = nw,
    partition    = partition,
    actor_labels = labels,
    group_labels = g_names
  )
}