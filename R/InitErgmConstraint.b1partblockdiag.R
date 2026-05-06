#' @templateVar name b1partblockdiag
#' @title Combined b1part + block-diagonal constraint for bipartite networks
#'
#' @description For bipartite networks organised in contiguous blocks (e.g. the
#'   stacked PLE meta-network produced by \code{erpm_long()}), enforce two
#'   constraints simultaneously:
#'   \enumerate{
#'     \item Each mode-1 node (actor) has degree exactly 1, so MCMC explores the
#'           space of partitions — same behaviour as \code{b1part}.
#'     \item Only dyads \eqn{(i,j)} where \code{attr(i) == attr(j)} are allowed,
#'           preventing cross-block edges — same behaviour as
#'           \code{blockdiag(attr)}.
#'   }
#'   The vertex attribute \code{attr} must label blocks contiguously within each
#'   bipartite mode.
#'
#' @usage
#' # b1partblockdiag(attr)
#' @template ergmTerm-attr
#'
#' @template ergmConstraint-general
#'
#' @concept bipartite
#' @concept directed
#' @concept undirected
#'
#' @importFrom ergm check.ErgmTerm ERGM_VATTR_SPEC ergm_get_vattr rlebdm
#' @importFrom network is.bipartite
#' @import rle
InitErgmConstraint.b1partblockdiag <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      bipartite     = TRUE,
                      varnames      = c("attr"),
                      vartypes      = c(ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL),
                      required      = c(TRUE))

  dbg <- isTRUE(getOption("ERPM.b1partblockdiag.debug", FALSE))
  .dbg <- function(...) if (dbg) message("[b1partblockdiag] ", ...)

  contigmsg <- paste0(
    "b1partblockdiag requires contiguous blocks within each bipartite mode. ",
    "See ergmConstraint?blockdiag for details.")

  list(
    attr       = a$attr,
    dependence = TRUE,
    implies    = c("b1degrees", "edges"),
    free_dyads = {
      n <- network.size(nw)
      storage.mode(n) <- "integer"
      bip <- nw %n% "bipartite"

      .dbg(sprintf("network: n=%d | mode-1 (actors)=%d | mode-2 (groups)=%d | attr='%s'",
                   n, bip, n - bip, paste(a$attr, collapse = ",")))

      av <- c(ergm_get_vattr(a$attr, nw))  # strip attributes that confuse rle()
      ea <- av[seq_len(bip)]               # mode-1 (actor) block labels
      aa <- av[bip + seq_len(n - bip)]     # mode-2 (group) block labels

      .dbg(sprintf("mode-1 block labels: %s", paste(base::rle(ea)$values, collapse = ", ")))
      .dbg(sprintf("mode-2 block labels: %s", paste(base::rle(aa)$values, collapse = ", ")))

      if (anyDuplicated(base::rle(ea)$values) || anyDuplicated(base::rle(aa)$values))
        stop(contigmsg)

      tmp <- .erpm_double_rle(ea, aa)
      el  <- tmp$lengths1   # block lengths in mode 1
      al  <- tmp$lengths2   # block lengths in mode 2

      B <- length(tmp$values)
      .dbg(sprintf("B=%d shared blocks: %s", B, paste(tmp$values, collapse = ", ")))
      if (dbg) {
        actor_starts <- c(0L, cumsum(el)[-B]) + 1L
        group_starts <- bip + c(0L, cumsum(al)[-B]) + 1L
        for (b in seq_len(B)) {
          fmt <- "[b1partblockdiag] block %-4s: actors %d..%d (n=%d) <-> groups %d..%d (n=%d) | dyads=%d" # nolint: line_length_linter.
          message(sprintf(fmt,
                          tmp$values[b],
                          actor_starts[b], actor_starts[b] + el[b] - 1L, el[b],
                          group_starts[b], group_starts[b] + al[b] - 1L, al[b],
                          el[b] * al[b]))
        }
      }

      # Bipartite blockdiag RLEBDM — mirrors ergm::InitErgmConstraint.blockdiag.
      # Uses the rle package's c.rle / rep.rle methods (loaded via @import rle).
      o <- rlebdm(
        c(rep(rle(FALSE), bip * n, scale = "run"),
          do.call(c, rep(
            Map(function(blen, bend) {
                  rep(rle(c(FALSE, TRUE, FALSE)),
                      c(bend - blen, blen, n - bend), scale = "run")
                },
                el, cumsum(el)),
            al))),
        n)

      ot <- rlebdm(
        c(do.call(c, rep(
            Map(function(blen, bend) {
                  rep(rle(c(FALSE, TRUE, FALSE)),
                      c(bip + bend - blen, blen, n - bip - bend), scale = "run")
                },
                al, cumsum(al)),
            el)),
          rep(rle(FALSE), (n - bip) * n, scale = "run")),
        n)

      fd <- compress(o | ot)

      .dbg(sprintf("free_dyads: total=%d (expected=%d) | class=%s",
                   sum(as.matrix(fd)),
                   sum(el * al) * 2L,
                   paste(class(fd), collapse = "/")))
      fd
    }
  )
}

# Safe reimplementation of ergm's internal .double.rle using base::rle explicitly,
# avoiding conflicts with the rle package's enhanced rle class.
.erpm_double_rle <- function(a1, a2) {
  e1 <- base::rle(a1)
  e2 <- base::rle(a2)
  o <- intersect(e1$values, e2$values)
  if (!all(e1$values[e1$values %in% o] == e2$values[e2$values %in% o]))
    stop("b1partblockdiag: common blocks must have the same order in both bipartite modes.")
  l1 <- e1$lengths[match(o, e1$values)]; l1[is.na(l1)] <- 0L
  l2 <- e2$lengths[match(o, e2$values)]; l2[is.na(l2)] <- 0L
  list(values = o, lengths1 = l1, lengths2 = l2)
}

#' @templateVar name B1PartBlockdiag
#' @aliases InitErgmProposal.B1PartBlockdiag
#' @title MHp for b1partblockdiag constraints
#'
#' @description MHp for \code{constraints = ~b1partblockdiag}. Reuses the
#'   \code{B1Part} C proposal; the \code{free_dyads} field from the constraint
#'   restricts sampling to within-block dyads.
#'
#' @template ergmProposal-general
#' @importFrom ergm ergm_Init_stop
#' @importFrom network is.bipartite
NULL
InitErgmProposal.B1PartBlockdiag <- function(arguments, nw) {
  if (!is.bipartite(nw))
    ergm_Init_stop("B1PartBlockdiag requires a bipartite network.")
  list(name = "B1Part", inputs = NULL)
}
