### HELPER FUNCTIONS ###

ref_cliques_group_sizes <- function(partition) {
  as.integer(table(partition))
}

ref_cliques <- function(partition, k = 2L, normalized = FALSE) {
  group_sizes <- ref_cliques_group_sizes(partition)
  k <- as.integer(k)

  values <- vapply(k, function(k_one) {
    if (k_one == 1L) {
      return(sum(group_sizes == 1L))
    }

    raw <- choose(group_sizes, k_one)

    if (isTRUE(normalized)) {
      return(sum(raw / group_sizes))
    }

    sum(raw)
  }, numeric(1))

  unname(values)
}

ref_cliques_GW <- function(partition, lambda = 2) {
  group_sizes <- ref_cliques_group_sizes(partition)

  values <- vapply(lambda, function(lambda_one) {
    r <- (lambda_one - 1) / lambda_one
    sum(lambda_one * (1 - r^group_sizes))
  }, numeric(1))

  unname(values)
}

cliques_network <- function(partition) {
  build_bipartite_from_inputs(partition = partition)$network
}

cliques_summary <- function(nw, rhs) {
  formula <- as.formula(call("~", quote(nw), rhs))
  environment(formula) <- environment()

  as.numeric(summary(formula))
}

### TERM: CLIQUES ###

test_that("cliques summary matches analytic reference values", {
  withr::local_options(ERPM.cliques.debug = FALSE)

  partitions <- list(
    c(1, 2, 2, 3, 3, 3),
    c(1, 1, 2, 3, 3, 4, 4, 4),
    c(1, 1, 1, 2, 2, 3),
    c(1, 2, 3, 4, 5),
    rep(1, 6)
  )

  cases <- list(
    list(rhs = quote(cliques), k = 2L, normalized = FALSE),
    list(rhs = quote(cliques(1)), k = 1L, normalized = FALSE),
    list(rhs = quote(cliques(2)), k = 2L, normalized = FALSE),
    list(rhs = quote(cliques(3)), k = 3L, normalized = FALSE),
    list(rhs = quote(cliques(k = c(1, 2, 3))), k = c(1L, 2L, 3L), normalized = FALSE),
    list(rhs = quote(cliques(clique_size = 3)), k = 3L, normalized = FALSE),
    list(rhs = quote(cliques(k = 2, normalized = TRUE)), k = 2L, normalized = TRUE),
    list(rhs = quote(cliques(k = 3, normalized = TRUE)), k = 3L, normalized = TRUE)
  )

  for (partition in partitions) {
    nw <- cliques_network(partition)

    for (case in cases) {
      observed <- cliques_summary(nw, case$rhs)
      expected <- ref_cliques(partition, k = case$k, normalized = case$normalized)

      expect_equal(observed, expected)
    }
  }
})

test_that("cliques validates arguments", {
  withr::local_options(ERPM.cliques.debug = FALSE)

  nw <- cliques_network(c(1, 1, 2, 3))

  expect_error(summary(nw ~ cliques(size = 2)), "argument 'size' is not supported", fixed = TRUE)
  expect_error(summary(nw ~ cliques(k = numeric(0))), "specify at least one value of k", fixed = TRUE)
  expect_error(summary(nw ~ cliques(k = c(1, NA))), "'k' must not contain NA", fixed = TRUE)
  expect_error(summary(nw ~ cliques(k = 0)), "'k' must contain integers >= 1", fixed = TRUE)
  expect_error(summary(nw ~ cliques(k = 1.5)), "'k' must contain integers >= 1", fixed = TRUE)
  expect_error(summary(nw ~ cliques(k = 2, normalized = c(TRUE, FALSE))), "'normalized' must be a scalar boolean", fixed = TRUE)
  expect_error(summary(nw ~ cliques(k = 2, normalized = NA)), "'normalized' must be a scalar boolean", fixed = TRUE)
})

### TERM: CLIQUES_GW ###

test_that("cliques_GW summary matches analytic reference values", {
  withr::local_options(ERPM.cliques_GW.debug = FALSE)

  partitions <- list(
    c(1, 2, 2, 3, 3, 3),
    c(1, 1, 2, 3, 3, 4, 4, 4),
    c(1, 1, 1, 2, 2, 3),
    c(1, 2, 3, 4, 5),
    rep(1, 6)
  )

  cases <- list(
    list(rhs = quote(cliques_GW), lambda = 2),
    list(rhs = quote(cliques_GW(lambda = 1)), lambda = 1),
    list(rhs = quote(cliques_GW(lambda = 1.5)), lambda = 1.5),
    list(rhs = quote(cliques_GW(lambda = 3)), lambda = 3),
    list(rhs = quote(cliques_GW(lambda = c(1, 2, 4))), lambda = c(1, 2, 4))
  )

  for (partition in partitions) {
    nw <- cliques_network(partition)

    for (case in cases) {
      observed <- cliques_summary(nw, case$rhs)
      expected <- ref_cliques_GW(partition, lambda = case$lambda)

      expect_equal(observed, expected)
    }
  }
})

test_that("cliques_GW validates arguments", {
  withr::local_options(ERPM.cliques_GW.debug = FALSE)

  nw <- cliques_network(c(1, 1, 2, 3))

  expect_error(summary(nw ~ cliques_GW(lambda = NA)), "'lambda' argument is not of any of the expected types: 'numeric'", fixed = TRUE)
  expect_error(summary(nw ~ cliques_GW(lambda = NaN)), "'lambda' must be finite", fixed = TRUE)
  expect_error(summary(nw ~ cliques_GW(lambda = Inf)), "'lambda' must be finite", fixed = TRUE)
  expect_error(summary(nw ~ cliques_GW(lambda = 0.5)), "'lambda' must be >= 1", fixed = TRUE)
})
