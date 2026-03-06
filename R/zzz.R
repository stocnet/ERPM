# ==============================================================================
# File    : R/zzz.R
# Purpose : Register ERPM Metropolis–Hastings proposals in ergm's proposal table
# ==============================================================================
#' Package load hook for ERPM
#'
#' This hook is executed when the package is loaded (both with a regular
#' 'library(ERPM)' and with 'devtools::load_all()'). Its role is to make sure
#' that ERPM-specific MH proposals are visible to 'ergm' under the
#' '~b1part' constraint.
#'
#' Concretely, we:
#'   - register 'ErpmToggleStep',
#'   - register 'ErpmSwapStep',
#'   - register 'ErpmMergeStep',
#'   - register 'ErpmSplitStep',
#'   - register 'ErpmMix',
#'   - ensure that the legacy 'B1Part' row is still present.
#'
#' Duplicate rows are avoided to keep repeated dev reloads clean.
#'
#' @param libname Path to the library.
#' @param pkgname Package name.
#' @keywords internal
.onLoad <- function(libname, pkgname) {

  # Register ERPM-specific proposals (identical to C idempotence ).
  .register_erpm_proposal("ErpmToggleStep", priority = 6,  weights = "default")
  .register_erpm_proposal("ErpmSwapStep",   priority = 6,  weights = "default")
  .register_erpm_proposal("ErpmMergeStep",  priority = 7,  weights = "default")
  .register_erpm_proposal("ErpmSplitStep",  priority = 7,  weights = "default")
  .register_erpm_proposal("ErpmMix",        priority = 9, weights = "default")

  # Keep legacy B1Part registration. 
  .RegisterProposals()

  # Optional verbose output for debugging proposal wiring.
  if (isTRUE(getOption("ERPM.zzz.verbose", FALSE))) {
    tab <- ergm::ergm_proposal_table()
    tab_erpm <- tab[tab$Package == "ERPM", ]
    packageStartupMessage(
      "[ERPM] proposal table rows (Package==ERPM):\n",
      paste(capture.output(tab_erpm), collapse = "\n")
    )
  }
}

#' Ensure legacy B1Part proposal row is present
#'
#' This function checks whether the row already exists
#' before attempting to add it again.
#' 
#' Note : The naming convention of this function follows the pre-existing code.
#'
#' @return Invisibly returns 'TRUE' if the row was added, 'FALSE' otherwise.
#' @keywords internal
.RegisterProposals <- function() {

  tab <- ergm::ergm_proposal_table()

  already <- any(
    tab$Proposal    == "B1Part"    &
      tab$Class       == "c"         &
      tab$Reference   == "Bernoulli" &
      tab$Constraints == "&b1part"   &
      tab$Priority    == 10          &
      tab$Weights     == "random"
  )

  if (!already) {
    ergm::ergm_proposal_table(
      Class       = "c",
      Reference   = "Bernoulli",
      Constraints = "&b1part",
      Priority    = 10,
      Weights     = "random",
      Proposal    = "B1Part"
    )
  }

  invisible(!already)
}

#' Register an ERPM proposal in ergm's proposal table
#'
#' This helper adds a row to 'ergm::ergm_proposal_table()' for a given
#' ERPM proposal (e.g. 'ErpmToggleStep', 'ErpmSwapStep') under the
#' '~b1part' constraint and Bernoulli reference.
#'
#' The function is safe under repeated calls,
#' as it checks whether an identical row already exists before inserting
#' a new one.
#'
#' @param proposal Character scalar. Name of the MH proposal.
#' @param priority Integer priority used by ergm’s proposal selection logic.
#' @param weights Character weight category (typically '"default"').
#'
#' @return Invisibly returns 'TRUE' if the row was added, 'FALSE' otherwise.
#' @keywords internal
.register_erpm_proposal <- function(proposal, priority, weights) {

  tab <- ergm::ergm_proposal_table()

  already <- any(
    tab$Proposal     == proposal      &
      tab$Package    == "ERPM"        &
      tab$Class      == "c"           &
      tab$Reference  == "Bernoulli"   &
      tab$Constraints == "&b1part"    &
      tab$Weights    == weights
  )

  if (already) {
    return(invisible(FALSE))
  }

  ergm::ergm_proposal_table(
    Class       = "c",
    Reference   = "Bernoulli",
    Constraints = "&b1part",
    Priority    = priority,
    Weights     = weights,
    Proposal    = proposal,
    Package     = "ERPM"
  )

  invisible(TRUE)
}

#' Package unload hook for ERPM
#'
#' Ensures that the ERPM shared library is properly unloaded when the
#' package is detached.
#'
#' @param libpath Path to the library.
#' @keywords internal
.onUnload <- function(libpath) {
  library.dynam.unload("ERPM", libpath)
}