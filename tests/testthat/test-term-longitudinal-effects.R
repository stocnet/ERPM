### HELPER FUNCTIONS ###

longitudinal_summary <- function(partitions, rhs) {
  formula <- as.formula(call("~", quote(partitions), rhs))
  environment(formula) <- environment()

  as.numeric(erpm_long(formula, eval.call = TRUE, verbose = FALSE, debug = FALSE, summary_test = TRUE))
}

ref_longitudinal_groups <- function(partition) {
  split(seq_along(partition), as.integer(partition))
}

ref_inertia_groups <- function(partitions, past_influence = 1L, size = NULL) {
  d <- as.integer(past_influence)
  total <- 0

  for (t in seq.int(d + 1L, length(partitions))) {
    current_groups <- ref_longitudinal_groups(partitions[[t]])

    for (current_group in current_groups) {
      if (!is.null(size) && !(length(current_group) %in% as.integer(size))) next

      matched <- FALSE

      for (lag in seq_len(d)) {
        past_groups <- ref_longitudinal_groups(partitions[[t - lag]])
        matched <- any(vapply(past_groups, function(past_group) {
          length(past_group) == length(current_group) && all(sort(past_group) == sort(current_group))
        }, logical(1)))

        if (matched) break
      }

      if (matched) total <- total + length(current_group)
    }
  }

  total
}

### TERM: INERTIA_GROUPS ###

test_that("inertia_groups summary matches analytic reference values", {
  withr::local_options(
    ERPM.inertia_groups.debug = FALSE,
    ERPM.b1partblockdiag.debug = FALSE
  )

  fixtures <- list(
    list(
      partitions = list(
        c(1, 1, 2, 2, 3, 3),
        c(1, 1, 2, 3, 3, 2),
        c(1, 1, 2, 3, 3, 2),
        c(1, 2, 2, 3, 3, 1)
      )
    ),
    list(
      partitions = list(
        c(1, 1, 1, 2, 2),
        c(1, 1, 2, 2, 2),
        c(1, 2, 2, 2, 1),
        c(1, 2, 2, 3, 3)
      )
    )
  )

  cases <- list(
    list(rhs = quote(inertia_groups(past_influence = 1)), d = 1, size = NULL),
    list(rhs = quote(inertia_groups(1)), d = 1, size = NULL),
    list(rhs = quote(inertia_groups(past_influence = 2)), d = 2, size = NULL),
    list(rhs = quote(inertia_groups(past_influence = 1, size = 2)), d = 1, size = 2),
    list(rhs = quote(inertia_groups(past_influence = 1, sizes = c(2, 3))), d = 1, size = c(2, 3))
  )

  for (fixture in fixtures) {
    for (case in cases) {
      observed <- longitudinal_summary(fixture$partitions, case$rhs)
      expected <- ref_inertia_groups(fixture$partitions, past_influence = case$d, size = case$size)

      #cat("\n")
      #cat("partitions:\n")
      #for (t in seq_along(fixture$partitions)) cat("  t", t, ":", paste(fixture$partitions[[t]], collapse = ", "), "\n", sep = "")
      #cat("term:", paste(deparse(case$rhs), collapse = " "), "\n")
      #cat("observed:", paste(observed, collapse = ", "), "\n")
      #cat("expected:", paste(expected, collapse = ", "), "\n")

      expect_equal(observed, expected)
    }
  }
})

test_that("inertia_groups validates longitudinal and term arguments", {
  withr::local_options(
    ERPM.inertia_groups.debug = FALSE,
    ERPM.b1partblockdiag.debug = FALSE
  )

  partitions <- list(
    c(1, 1, 2, 2),
    c(1, 2, 2, 1),
    c(1, 2, 1, 2)
  )

  expect_error(
    erpm_long(partitions ~ inertia_groups(past_influence = 0), summary_test = TRUE),
    "past_influence must be >= 1",
    fixed = TRUE
  )
  expect_error(
    erpm_long(partitions ~ inertia_groups(past_influence = 3), summary_test = TRUE),
    "need at least T >= past_influence + 1",
    fixed = TRUE
  )
  #expect_error(
  #  erpm_long(partitions ~ inertia_groups(type = "endogenous"), summary_test = TRUE),
  #  "not implemented yet",
  #  fixed = TRUE
  #)
  #expect_error(
  #  erpm_long(partitions ~ inertia_groups(size = 0), summary_test = TRUE),
  #  "`size` must contain integers > 0",
  #  fixed = TRUE
  #)
  expect_error(
    erpm_long(partitions ~ inertia_groups(size = 1.5), summary_test = TRUE),
    "`size` must be integer-valued",
    fixed = TRUE
  )
})
