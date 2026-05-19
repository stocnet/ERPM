### HELPER FUNCTIONS ###

dyadic_covariate_network <- function(partition, dyads) {
  nodes <- data.frame(label = paste0("A", seq_along(partition)), stringsAsFactors = FALSE)
  build_bipartite_from_inputs(partition = partition, nodes = nodes, dyads = dyads)$network
}

dyadic_covariate_summary <- function(nw, rhs) {
  formula <- as.formula(call("~", quote(nw), rhs))
  environment(formula) <- environment()

  as.numeric(summary(formula))
}

ref_dyadic_groups <- function(partition) {
  split(seq_along(partition), as.integer(partition))
}

ref_dyadic_clique_sum <- function(indices, Z, k) {
  if (length(indices) < k) return(0)

  combinations <- utils::combn(indices, k)
  sum(apply(combinations, 2L, function(actors) {
    pairs <- utils::combn(actors, 2L)
    prod(apply(pairs, 2L, function(pair) Z[pair[1], pair[2]] + Z[pair[2], pair[1]]))
  }))
}

ref_dyadcov <- function(partition, Z, clique_size = 2L, normalize = c("none", "global", "by_group")) {
  normalize <- match.arg(normalize)
  groups <- ref_dyadic_groups(partition)
  k <- as.integer(clique_size)
  total <- 0

  for (indices in groups) {
    group_size <- length(indices)
    if (group_size < k) next

    group_total <- ref_dyadic_clique_sum(indices, Z, k)

    if (normalize == "none") {
      total <- total + group_total
    } else if (normalize == "global") {
      total <- total + group_total / group_size
    } else {
      total <- total + group_total / choose(group_size, k)
    }
  }

  total
}

ref_dyadcov_GW <- function(partition, Z, lambda = 2) {
  groups <- ref_dyadic_groups(partition)
  total <- 0

  for (indices in groups) {
    group_size <- length(indices)
    if (group_size < 2L) next

    factor <- 1
    for (k in 2:group_size) {
      total <- total + factor * ref_dyadic_clique_sum(indices, Z, k)
      factor <- factor * (-1 / lambda)
    }
  }

  total
}

ref_dyadcov_full <- function(partition, Z, size = NULL) {
  groups <- ref_dyadic_groups(partition)
  total <- 0

  for (indices in groups) {
    group_size <- length(indices)
    if (group_size < 2L) next
    if (!is.null(size) && !(group_size %in% as.integer(size))) next

    ordered_pairs <- expand.grid(i = indices, j = indices)
    ordered_pairs <- ordered_pairs[ordered_pairs$i != ordered_pairs$j, , drop = FALSE]
    total <- total + sum(Z[cbind(ordered_pairs$i, ordered_pairs$j)])
  }

  total
}

dyadic_covariate_fixtures <- function() {
  list(
    list(
      partition = c(1, 1, 1, 2, 2, 3),
      Z1 = matrix(c(
        0, 1, 2, 0, 1, 3,
        4, 0, 1, 2, 0, 1,
        2, 3, 0, 1, 4, 0,
        1, 0, 2, 0, 3, 1,
        0, 2, 1, 4, 0, 2,
        3, 1, 0, 2, 1, 0
      ), 6, 6, byrow = TRUE),
      Z2 = matrix(c(
        0, 2, 1, 1, 0, 2,
        1, 0, 3, 0, 2, 1,
        4, 1, 0, 2, 1, 0,
        0, 3, 1, 0, 2, 2,
        2, 1, 0, 1, 0, 3,
        1, 0, 2, 4, 1, 0
      ), 6, 6, byrow = TRUE)
    ),
    list(
      partition = c(1, 1, 2, 2, 2, 3, 3),
      Z1 = matrix(c(
        0, 1, 0, 2, 1, 3, 0,
        2, 0, 1, 0, 3, 1, 2,
        1, 4, 0, 2, 0, 1, 3,
        0, 1, 3, 0, 2, 0, 1,
        2, 0, 1, 4, 0, 2, 0,
        1, 3, 0, 1, 2, 0, 4,
        0, 2, 1, 3, 0, 1, 0
      ), 7, 7, byrow = TRUE),
      Z2 = matrix(c(
        0, 3, 2, 0, 1, 2, 1,
        1, 0, 0, 2, 3, 0, 2,
        2, 1, 0, 3, 0, 1, 0,
        1, 0, 2, 0, 4, 2, 1,
        0, 2, 1, 3, 0, 1, 2,
        3, 1, 0, 2, 1, 0, 4,
        1, 0, 3, 1, 2, 0, 0
      ), 7, 7, byrow = TRUE)
    )
  )
}

### TERM: DYADCOV ###

test_that("dyadcov summary matches analytic reference values", {
  withr::local_options(ERPM.dyadcov.debug = FALSE)

  cases <- list(
    list(rhs = quote(dyadcov("Z1", clique_size = 2, normalize = "none")), dyad = "Z1", k = 2L, normalize = "none"),
    list(rhs = quote(dyadcov("Z1", clique_size = 2, normalize = "global")), dyad = "Z1", k = 2L, normalize = "global"),
    list(rhs = quote(dyadcov("Z1", clique_size = 3, normalize = "none")), dyad = "Z1", k = 3L, normalize = "none"),
    list(rhs = quote(dyadcov("Z1", clique_size = 3, normalize = "by_group")), dyad = "Z1", k = 3L, normalize = "by_group"),
    list(rhs = quote(dyadcov("Z2", clique_size = 2, normalize = "none")), dyad = "Z2", k = 2L, normalize = "none")
  )

  for (fixture in dyadic_covariate_fixtures()) {
    nw <- dyadic_covariate_network(fixture$partition, dyads = list(Z1 = fixture$Z1, Z2 = fixture$Z2))

    for (case in cases) {
      observed <- dyadic_covariate_summary(nw, case$rhs)
      expected <- ref_dyadcov(fixture$partition, fixture[[case$dyad]], clique_size = case$k, normalize = case$normalize)

      #cat("\n")
      #cat("partition:", paste(fixture$partition, collapse = ", "), "\n")
      #cat("dyad:", case$dyad, "\n")
      #cat("term:", paste(deparse(case$rhs), collapse = " "), "\n")
      #cat("observed:", paste(observed, collapse = ", "), "\n")
      #cat("expected:", paste(expected, collapse = ", "), "\n")

      expect_equal(observed, expected)
    }
  }
})

test_that("dyadcov validates arguments", {
  withr::local_options(ERPM.dyadcov.debug = FALSE)

  partition <- c(1, 1, 2, 2)
  Z <- matrix(1, 4, 4)
  diag(Z) <- 0
  Z_bad <- Z
  Z_bad[1, 2] <- NA
  nw <- dyadic_covariate_network(partition, dyads = list(Z = Z, Z_bad = Z_bad))

  expect_error(summary(nw ~ dyadcov("missing")), "dyad not found", fixed = TRUE)
  expect_error(summary(nw ~ dyadcov("Z_bad")), "NA values are not allowed", fixed = TRUE)
  expect_error(summary(nw ~ dyadcov("Z", clique_size = 1)), "'clique_size' must be an integer >= 2", fixed = TRUE)
  expect_error(summary(nw ~ dyadcov("Z", normalize = "bad")), "should be one of", fixed = TRUE)
})

### TERM: DYADCOV_GW ###

test_that("dyadcov_GW summary matches analytic reference values", {
  withr::local_options(ERPM.dyadcov_GW.debug = FALSE)

  cases <- list(
    list(rhs = quote(dyadcov_GW("Z1", lambda = 2)), dyad = "Z1", lambda = 2),
    list(rhs = quote(dyadcov_GW("Z1", lambda = 3)), dyad = "Z1", lambda = 3),
    list(rhs = quote(dyadcov_GW("Z2", lambda = 2)), dyad = "Z2", lambda = 2),
    list(rhs = quote(dyadcov_GW("Z2", lambda = 1.5)), dyad = "Z2", lambda = 1.5)
  )

  for (fixture in dyadic_covariate_fixtures()) {
    nw <- dyadic_covariate_network(fixture$partition, dyads = list(Z1 = fixture$Z1, Z2 = fixture$Z2))

    for (case in cases) {
      observed <- dyadic_covariate_summary(nw, case$rhs)
      expected <- ref_dyadcov_GW(fixture$partition, fixture[[case$dyad]], lambda = case$lambda)

      #cat("\n")
      #cat("partition:", paste(fixture$partition, collapse = ", "), "\n")
      #cat("dyad:", case$dyad, "\n")
      #cat("term:", paste(deparse(case$rhs), collapse = " "), "\n")
      #cat("observed:", paste(observed, collapse = ", "), "\n")
      #cat("expected:", paste(expected, collapse = ", "), "\n")

      expect_equal(observed, expected)
    }
  }
})

test_that("dyadcov_GW validates arguments", {
  withr::local_options(ERPM.dyadcov_GW.debug = FALSE)

  partition <- c(1, 1, 2, 2)
  Z <- matrix(1, 4, 4)
  diag(Z) <- 0
  Z_bad <- Z
  Z_bad[1, 2] <- NA
  nw <- dyadic_covariate_network(partition, dyads = list(Z = Z, Z_bad = Z_bad))

  expect_error(summary(nw ~ dyadcov_GW("missing")), "not found", fixed = TRUE)
  expect_error(summary(nw ~ dyadcov_GW("Z_bad")), "NA values are not allowed", fixed = TRUE)
  #expect_error(summary(nw ~ dyadcov_GW("Z", lambda = c(2, 3))), "'lambda' must be a scalar", fixed = TRUE)
  expect_error(summary(nw ~ dyadcov_GW("Z", lambda = 0)), "strictly positive", fixed = TRUE)
})

### TERM: DYADCOV_FULL ###

test_that("dyadcov_full summary matches analytic reference values", {
  withr::local_options(ERPM.dyadcov_full.debug = FALSE)

  cases <- list(
    list(rhs = quote(dyadcov_full("Z1")), dyad = "Z1", size = NULL),
    list(rhs = quote(dyadcov_full("Z1", size = 2)), dyad = "Z1", size = 2),
    list(rhs = quote(dyadcov_full("Z1", size = 2:3)), dyad = "Z1", size = 2:3),
    list(rhs = quote(dyadcov_full("Z2")), dyad = "Z2", size = NULL),
    list(rhs = quote(dyadcov_full("Z2", size = 3)), dyad = "Z2", size = 3)
  )

  for (fixture in dyadic_covariate_fixtures()) {
    nw <- dyadic_covariate_network(fixture$partition, dyads = list(Z1 = fixture$Z1, Z2 = fixture$Z2))

    for (case in cases) {
      observed <- dyadic_covariate_summary(nw, case$rhs)
      expected <- ref_dyadcov_full(fixture$partition, fixture[[case$dyad]], size = case$size)

      #cat("\n")
      #cat("partition:", paste(fixture$partition, collapse = ", "), "\n")
      #cat("dyad:", case$dyad, "\n")
      #cat("term:", paste(deparse(case$rhs), collapse = " "), "\n")
      #cat("observed:", paste(observed, collapse = ", "), "\n")
      #cat("expected:", paste(expected, collapse = ", "), "\n")

      expect_equal(observed, expected)
    }
  }
})

test_that("dyadcov_full validates arguments", {
  withr::local_options(ERPM.dyadcov_full.debug = FALSE)

  partition <- c(1, 1, 2, 2)
  Z <- matrix(1, 4, 4)
  diag(Z) <- 0
  Z_bad <- Z
  Z_bad[1, 2] <- NA
  nw <- dyadic_covariate_network(partition, dyads = list(Z = Z, Z_bad = Z_bad))

  #expect_error(summary(nw ~ dyadcov_full("missing")), "dyadic matrix not found", fixed = TRUE)
  expect_error(summary(nw ~ dyadcov_full("Z_bad")), "NA values are not allowed", fixed = TRUE)
  expect_error(summary(nw ~ dyadcov_full("Z", size = 0)), "'size' must contain positive integers", fixed = TRUE)
  expect_error(summary(nw ~ dyadcov_full("Z", size = 1.5)), "'size' must contain positive integers", fixed = TRUE)
})
