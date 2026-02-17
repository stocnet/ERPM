# ==============================================================================
# File    : InitErgmTerm.squared_sizes.R
# Purpose : Declare the ERGM term 'squared_sizes' for bipartite networks
#           (aggregation over the group mode).
# Project : ERPM / ERGM extensions
# ============================================================================

#' ERGM term: squared_sizes (group-mode degrees raised to a power)
#'
#' @name InitErgmTerm.squared_sizes
#' @aliases squared_sizes
#' @note InitErgmTerm.squared_sizes.R
#'
#' @description
#' \code{squared_sizes} is an ERGM term for bipartite networks that aggregates
#' powers of group sizes. The network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side of the bipartite graph).
#' }
#'
#' For each node in the group mode with degree \eqn{\deg(g)} (number of adjacent
#' actors), the term can accumulate contributions of the form
#' \deqn{
#'   \deg(g)^{\text{pow}},
#' }
#' restricted to a set of admissible group sizes. The argument \code{sizes}
#' specifies a (possibly multi-valued) set of group sizes over which the
#' contributions are aggregated into a single statistic.
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic.
#'
#' IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
#' - This term is intended to support multi-toggle moves (swap/split/merge
#'   represented as a list of toggles).
#' - Therefore, the compiled change-statistic must be implemented using the
#'   D_CHANGESTAT_FN API (multi-toggle).
#' - On the R side, we MUST advertise this to ergm by returning `d_func = TRUE`.
#'   Otherwise ergm will try to call the changestat as a one-toggle C_CHANGESTAT_FN,
#'   causing a signature mismatch and typically a segfault.
#'
#' Compiled symbol naming convention:
#' - Recommended: implement the C function as `d_squared_sizes` via
#'   D_CHANGESTAT_FN(d_squared_sizes).
#' - Avoid exposing a symbol named `c_squared_sizes` with a D-signature, because
#'   ergm may resolve it as the one-toggle entrypoint and crash.
#'
#' The R initializer:
#' \itemize{
#'   \item enforces that the network is bipartite (actor mode vs group mode);
#'   \item parses the user arguments \code{sizes} and \code{pow};
#'   \item checks that \code{sizes} are integer group sizes \eqn{\ge 1};
#'   \item checks that \code{pow} is a single integer and at least 1;
#'   \item builds a compact \code{inputs} vector for the C code encoding:
#'         \code{pow}, the number of admissible sizes, then the list of sizes.
#' }
#'
#' @export
InitErgmTerm.squared_sizes <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "squared_sizes"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  # Global option:
  #   options(ERPM.squared_sizes.debug = TRUE/FALSE)
  # When TRUE, the initializer prints diagnostic messages to the console.
  dbg    <- isTRUE(getOption("ERPM.squared_sizes.debug", TRUE))
  dbgcat <- function(...) if (dbg) cat("[squared_sizes][DEBUG]", ..., "\n", sep = "")

  dbgcat("InitErgmTerm.squared_sizes called with args: ",
         paste(names(arglist), collapse = ", "))

  # Guard:  typo 'size' instead of 'sizes'
  if ("size" %in% names(arglist) && !"sizes" %in% names(arglist)) {
    ergm_Init_stop(
      sQuote(termname),
      ": argument 'size' is not supported; did you mean 'sizes'?"
    )
  }

  # ---------------------------------------------------------------------------
  # Base validation and structural requirements
  # ---------------------------------------------------------------------------
  # - Enforce a bipartite network (actor mode / group mode).
  # - Declare the allowed arguments: sizes, pow.
  # - Provide default values: sizes = NULL (all sizes), pow = 2.
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      =    c("sizes",  "pow"),
    vartypes      =    c("numeric","numeric"),
    defaultvalues = list(NULL,     2),
    required      =    c(FALSE,    FALSE)
  )

  sizes <- a$sizes
  pow   <- a$pow

  # ---------------------------------------------------------------------------
  # Default sizes: all non-empty group sizes up to the actor-mode size
  # ---------------------------------------------------------------------------
  if (is.null(sizes) || length(sizes) == 0L) {
    n1 <- network::get.network.attribute(nw, "bipartite")
    if (is.null(n1) || is.na(n1)) {
      ergm_Init_stop(
        sQuote(termname),
        ": network must have a valid 'bipartite' attribute (actor-mode size)."
      )
    }
    n1 <- as.integer(n1)
    if (!is.finite(n1) || n1 < 1L) {
      ergm_Init_stop(
        sQuote(termname),
        ": 'bipartite' attribute must be a positive finite integer."
      )
    }
    dbgcat("bipartite attribute (n1) =", n1,
           " | using default sizes = 1:", n1)
    sizes <- seq_len(n1)
  } else {
    dbgcat("user-specified sizes (raw) = ",
           paste(sizes, collapse = ","))
  }

  # ---------------------------------------------------------------------------
  # Constraints and normalization on sizes
  # ---------------------------------------------------------------------------
  sizes <- as.numeric(sizes)
  if (any(is.na(sizes))) {
    ergm_Init_stop(sQuote(termname), ": 'sizes' must not contain NA.")
  }
  if (any(sizes < 1 | sizes != as.integer(sizes))) {
    ergm_Init_stop(sQuote(termname), ": 'sizes' must be integer >= 1.")
  }
  sizes <- as.integer(sizes)
  dbgcat("sizes (validated) = ", paste(sizes, collapse = ","))

  # ---------------------------------------------------------------------------
  # Constraints on pow (scalar)
  # ---------------------------------------------------------------------------
  if (length(pow) == 0L) {
    pow <- 2L
  }
  if (length(pow) > 1L) {
    ergm_Init_stop(sQuote(termname), ": 'pow' must be of length 1.")
  }
  if (any(pow < 1 | pow != as.integer(pow))) {
    ergm_Init_stop(sQuote(termname), ": 'pow' must be an integer >= 1.")
  }
  pow <- as.integer(pow[1L])
  dbgcat("pow (validated)   = ", pow)

  # ---------------------------------------------------------------------------
  # Trivial case: no sizes => no statistic
  # ---------------------------------------------------------------------------
  if (length(sizes) == 0L) {
    dbgcat("no sizes after validation -> returning NULL term")
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Coefficient name
  # ---------------------------------------------------------------------------
  if (length(sizes) == 1L) {
    size_tag <- paste0("size", sizes)
  } else {
    size_tag <- paste0("sizes", paste(sizes, collapse = "_"))
  }
  coef.name <- paste0(
    "squared_", size_tag,
    ifelse(pow != 1L, paste0("_pow", pow), "")
  )
  dbgcat("coef.names = ", coef.name)

  # ---------------------------------------------------------------------------
  # Compact INPUT_PARAM vector for the C change-statistic
  # ---------------------------------------------------------------------------
  inputs <- c(
    as.double(pow),
    as.double(length(sizes)),
    as.double(sizes)
  )
  dbgcat("inputs length = ", length(inputs),
         " (1 aggregated stat over ", length(sizes), " sizes)")

  # ---------------------------------------------------------------------------
  # Standard ERGM term initialization return value
  # ---------------------------------------------------------------------------
  # IMPORTANT:
  # - `d_func = TRUE` tells ergm to call the multi-toggle (D_) changestat entrypoint.
  # - Without it, ergm assumes a one-toggle C_ changestat and will call the function
  #   with the wrong signature if you compiled only a D_ function (=> segfault).
  list(
    name         = "squared_sizes",
    coef.names   = coef.name,
    inputs       = inputs,
    dependence   = TRUE,
    d_func       = TRUE,   # <-- REQUIRED for D_CHANGESTAT_FN (multi-toggle)
    emptynwstats = 0
  )
}