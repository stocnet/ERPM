erpm_long_api_partitions <- function() {
  list(
    c(1, 1, 2, 2),
    c(1, 2, 2, 1),
    c(1, 2, 1, 2)
  )
}

erpm_long_api_nodes <- function(partitions) {
  lapply(seq_along(partitions), function(t) {
    data.frame(
      label = paste0("actor", seq_along(partitions[[t]]), "_t", t),
      gender = c("F", "M", "F", "M")[seq_along(partitions[[t]])],
      age = seq_along(partitions[[t]]) + 20 + t,
      stringsAsFactors = FALSE
    )
  })
}

erpm_long_api_dyads <- function(partitions) {
  lapply(seq_along(partitions), function(t) {
    n <- length(partitions[[t]])
    Z <- matrix(seq_len(n * n) + t, n, n)
    diag(Z) <- 0
    list(Z1 = Z)
  })
}

erpm_long_api_call <- function(partitions,
                               rhs = quote(cliques(k = 2)),
                               nodes = NULL,
                               dyads = NULL,
                               estimate = NULL,
                               constraints = NULL) {
  formula <- as.formula(call("~", quote(partitions), rhs))
  environment(formula) <- environment()

  erpm_long(
    formula,
    eval.call = TRUE,
    verbose = FALSE,
    debug = FALSE,
    nodes = nodes,
    dyads = dyads,
    estimate = estimate,
    constraints = constraints
  )
}

erpm_long_api_formula_network <- function(call) {
  formula <- call[[2]]
  eval(formula[[2]], envir = environment(formula))
}

test_that("erpm_long returns an unevaluated erpm call in dry-run mode", {
  partitions <- erpm_long_api_partitions()
  call <- erpm_long_api_call(partitions, rhs = quote(cliques(k = 2)))

  expect_true(is.call(call))
  expect_identical(call[[1]], as.name("erpm"))
  expect_true(inherits(call[[2]], "formula"))
  expect_false(call$verbose)
  expect_true(call$eval.call)

  meta_nw <- erpm_long_api_formula_network(call)
  expect_s3_class(meta_nw, "network")
  expect_equal(network::get.network.attribute(meta_nw, "erpm_mode"), "empile")
  expect_equal(network::get.network.attribute(meta_nw, "erpm_long.selected_partition_indices"), seq_along(partitions))
  expect_equal(network::get.network.attribute(meta_nw, "bipartite"), sum(lengths(partitions)))
})

test_that("erpm_long dry-run call preserves node and dyadic inputs in the meta-network", {
  partitions <- erpm_long_api_partitions()
  nodes <- erpm_long_api_nodes(partitions)
  dyads <- erpm_long_api_dyads(partitions)

  call <- erpm_long_api_call(partitions, rhs = quote(dyadcov("Z1") + cov_match("gender")), nodes = nodes, dyads = dyads)
  meta_nw <- erpm_long_api_formula_network(call)
  stored_dyads <- network::get.network.attribute(meta_nw, "dyads")

  expect_equal(network::get.vertex.attribute(meta_nw, "gender")[seq_len(sum(lengths(partitions)))], unlist(lapply(nodes, `[[`, "gender"), use.names = FALSE))
  expect_named(stored_dyads, "Z1")
  expect_equal(dim(stored_dyads$Z1), rep(sum(lengths(partitions)), 2))
})

test_that("erpm_long forwards fitting options into the composed erpm call", {
  partitions <- erpm_long_api_partitions()
  call <- erpm_long_api_call(partitions, rhs = quote(groups), estimate = "MPLE", constraints = ~ b1part)

  expect_equal(call$estimate, "MPLE")
  expect_equal(deparse(call$constraints), deparse(~ b1part))
})

test_that("erpm_long public examples are stable", {
  skip("TODO: expand after erpm_long documentation and vignettes are finalized.")
})
