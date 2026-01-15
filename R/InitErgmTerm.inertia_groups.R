################################################################################
# FILE: R/InitErgmTerm.inertia_groups.R
################################################################################
#' ERGM term: inertia_groups (longitudinal exact-group persistence)
#'
#' @name InitErgmTerm.inertia_groups
#' @aliases inertia_groups
#' @note InitErgmTerm.inertia_groups.R
#'
#' @description
#' \code{inertia_groups} is an ERGM term for bipartite membership networks used
#' by \code{erpm_long()}. It counts, at time \eqn{t}, how many current groups
#' (group-mode vertices) have an \emph{exact} actor membership set that matches
#' at least one group observed in the past, over a window of \code{past_influence}
#' lags. The required past information is expected to be attached to the network
#' as network attributes by \code{erpm_long()}.
#'
#' The term is inertial: it does not compute history by itself. It requires
#' network attributes of the form:
#' \preformatted{
#'   nw %n% "erpm_inertia__inertia_groups__lag1"
#'   nw %n% "erpm_inertia__inertia_groups__lag2"
#'   ...
#'   nw %n% "erpm_inertia__inertia_groups__lagd"
#' }
#' where \code{d = past_influence}.
#'
#' Each lag attribute must be a list with at least:
#' \itemize{
#'   \item \code{type = "group_signature_set"}
#'   \item \code{signatures}: character vector, each entry like "1,3,5"
#' }
#'
#' Optional \code{size} restricts which current group sizes are eligible to be
#' counted (same semantics as other ERPM group-level effects): NULL = all sizes.
#'
#' @param nw A bipartite \pkg{network} object for time \eqn{t}.
#' @param arglist Term arguments from \pkg{ergm}. Supported:
#'   \itemize{
#'     \item \code{past_influence}: integer >= 1 (default 1).
#'     \item \code{size}: NULL or positive integers (filter on current group size).
#'   }
#' @param ... Unused.
#' @param version ERGM API version.
#'
#' @return A term specification list for \pkg{ergm}.
#' @export
InitErgmTerm.inertia_groups <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "inertia_groups"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  # Global option:
  #   options(ERPM.inertia_groups.debug = TRUE/FALSE)
  # When TRUE, the initializer prints diagnostic messages to the console.
  dbg    <- isTRUE(getOption("ERPM.inertia_groups.debug", TRUE))
  dbgcat <- function(...) if (dbg) cat("[inertia_groups][DEBUG]", ..., "\n", sep = "")

  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("past_influence", "size"),
    vartypes      = c("numeric,integer", "numeric,integer"),
    defaultvalues = list(1L, NULL),
    required      = c(FALSE, FALSE)
  )

  # -- bipartite boundary (actor mode size) -----------------------------------
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L) stop(termname, ": strictly bipartite network required.")
  dbgcat("n1 = ", n1)

  # -- past_influence d --------------------------------------------------------
  d <- a$past_influence
  if (length(d) != 1L || is.na(d) || !is.finite(d))
    stop(termname, ": 'past_influence' must be a finite scalar.")
  d <- as.integer(round(d))
  if (d < 1L) stop(termname, ": 'past_influence' must be >= 1.")
  dbgcat("past_influence (d) = ", d)

  # -- size filter on CURRENT groups ------------------------------------------
  sizes <- a$size
  if (is.null(sizes)) {
    L <- 0L
    sizes_vec <- numeric(0)
    dbgcat("size filter: NULL (all current group sizes eligible)")
  } else {
    if (!is.numeric(sizes) || length(sizes) == 0L)
      stop(termname, ": 'size' must be NULL or a non-empty numeric/integer vector.")
    sizes <- as.integer(round(sizes))
    if (any(!is.finite(sizes)) || any(sizes <= 0L))
      stop(termname, ": 'size' must contain positive integers.")
    sizes <- sort(unique(sizes))
    L <- length(sizes)
    sizes_vec <- as.double(sizes)
    dbgcat("size filter: S = {", paste(sizes, collapse = ","), "} (L = ", L, ")")
  }

  # -- helpers ----------------------------------------------------------------
  .stop_bad_attr <- function(msg) stop(paste0(termname, ": ", msg), call. = FALSE)

  .parse_signature <- function(sig, n1) {
    # Expect "1,3,5" (no spaces). Allow empty only if truly empty group (we disallow here).
    if (!is.character(sig) || length(sig) != 1L || is.na(sig) || !nzchar(sig))
      .stop_bad_attr("invalid signature entry (must be a non-empty string).")
    parts <- strsplit(sig, ",", fixed = TRUE)[[1L]]
    if (!length(parts)) .stop_bad_attr("invalid signature: empty.")
    ids <- suppressWarnings(as.integer(parts))
    if (anyNA(ids)) .stop_bad_attr(sprintf("invalid signature (non-integer token): %s", sig))
    if (any(ids < 1L | ids > n1)) .stop_bad_attr(sprintf("signature ids out of range 1..n1: %s", sig))
    ids <- sort(unique(ids))
    ids
  }

  # -- read and validate lag attributes; pack past groups ----------------------
  # We build a single numeric INPUT_PARAM vector:
  #
  # [1]  n1
  # [2]  d
  # [3]  L
  # [4..] sizes[L]
  # then for each lag=1..d:
  #   [ ] M_lag
  #   then for each group j=1..M_lag:
  #       len_j
  #       actor_ids (len_j entries)
  #
  # All stored as doubles, cast to int in C.
  inputs <- c(as.double(n1), as.double(d), as.double(L), sizes_vec)

  coef.name <- sprintf("inertia_groups[pi=%d]%s",
                       d,
                       if (L > 0L) sprintf("_S{%s}", paste(sizes, collapse = ",")) else "_all")
  dbgcat("coef.name = ", coef.name)

  for (lag in seq_len(d)) {
    attr_name <- paste0("erpm_inertia__", termname, "__lag", lag)
    obj <- nw %n% attr_name
    dbgcat("reading attribute: ", sQuote(attr_name))

    if (is.null(obj))
      .stop_bad_attr(sprintf("missing required network attribute %s (not attached by erpm_long?).", sQuote(attr_name)))
    if (!is.list(obj))
      .stop_bad_attr(sprintf("attribute %s must be a list.", sQuote(attr_name)))

    if (is.null(obj$type) || !identical(as.character(obj$type), "group_signature_set"))
      .stop_bad_attr(sprintf("attribute %s has wrong type (expected 'group_signature_set').", sQuote(attr_name)))

    sigs <- obj$signatures
    if (is.null(sigs) || !is.character(sigs))
      .stop_bad_attr(sprintf("attribute %s must contain a character vector $signatures.", sQuote(attr_name)))

    # Past groups can include empty groups in principle, but signatures should not be empty.
    # If the attribute has length 0, it simply means "no groups in the past".
    M <- length(sigs)
    inputs <- c(inputs, as.double(M))
    dbgcat("  lag=", lag, " | past groups M = ", M)

    if (M > 0L) {
      for (s in sigs) {
        ids <- .parse_signature(s, n1)
        if (length(ids) < 1L)
          .stop_bad_attr(sprintf("empty group signature is not allowed in %s.", sQuote(attr_name)))
        inputs <- c(inputs, as.double(length(ids)), as.double(ids))
        dbgcat("    sig=", sQuote(s), " | len=", length(ids))
      }
    }
  }

  dbgcat("inputs packed: length=", length(inputs))

  list(
    name         = "inertia_groups",
    coef.names   = coef.name,
    inputs       = inputs,
    dependence   = TRUE,
    minval       = 0,
    maxval       = Inf,
    emptynwstats = 0
  )
}
