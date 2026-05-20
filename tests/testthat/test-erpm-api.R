### HELPER FUNCTIONS ###

erpm_api_partition <- function() {
  c(1, 1, 2, 3)
}

erpm_api_nodes <- function() {
  data.frame(
    label = c("alice", "bob", "chloe", "david"),
    gender = c("F", "M", "F", "M"),
    age = c(21, 22, 31, 41),
    stringsAsFactors = FALSE
  )
}

erpm_api_dyad <- function() {
  matrix(c(
    0, 1, 2, 3,
    1, 0, 4, 5,
    2, 4, 0, 6,
    3, 5, 6, 0
  ), 4, 4, byrow = TRUE)
}

erpm_api_call <- function(partition,
                          rhs = quote(groups),
                          nodes = NULL,
                          dyads = list(),
                          group_labels = NULL,
                          constraints = NULL,
                          estimate = NULL) {
  formula <- as.formula(call("~", quote(partition), rhs))
  environment(formula) <- environment()

  erpm(
    formula,
    eval.call = FALSE,
    verbose = FALSE,
    nodes = nodes,
    dyads = dyads,
    group_labels = group_labels,
    constraints = constraints,
    estimate = estimate
  )
}

erpm_api_call_network <- function(nw, rhs = quote(groups)) {
  formula <- as.formula(call("~", quote(nw), rhs))
  environment(formula) <- environment()

  erpm(formula, eval.call = FALSE, verbose = FALSE)
}

erpm_api_formula_network <- function(call) {
  formula <- call[[2]]
  eval(formula[[2]], envir = environment(formula))
}

### ERPM PUBLIC API / WRAPPER ###

test_that("erpm returns an unevaluated ergm call for partition LHS when eval.call is FALSE", {
  partition <- erpm_api_partition()
  call <- erpm_api_call(partition, rhs = quote(groups))

  expect_true(is.call(call))
  expect_identical(call[[1]], as.name("ergm"))
  expect_true(inherits(call[[2]], "formula"))
  expect_equal(deparse(call$constraints), deparse(~ b1part))
  expect_false(call$verbose)

  nw <- erpm_api_formula_network(call)
  expect_s3_class(nw, "network")
  expect_equal(network::get.network.attribute(nw, "bipartite"), length(partition))
  expect_equal(network::network.edgecount(nw), length(partition))
})

test_that("erpm accepts a pre-built bipartite network LHS", {
  partition <- erpm_api_partition()
  built <- build_bipartite_from_inputs(partition)
  nw <- built$network

  call <- erpm_api_call_network(nw, rhs = quote(groups))
  call_nw <- erpm_api_formula_network(call)

  expect_true(is.call(call))
  expect_s3_class(call_nw, "network")
  expect_equal(network::get.network.attribute(call_nw, "bipartite"), length(partition))
  expect_equal(as.matrix(call_nw, matrix.type = "adjacency"), as.matrix(nw, matrix.type = "adjacency"))
})

test_that("erpm passes nodes and group labels through partition LHS construction", {
  partition <- erpm_api_partition()
  nodes <- erpm_api_nodes()
  group_labels <- c("G_A", "G_B", "G_C", "G_D")

  call <- erpm_api_call(partition, rhs = quote(cov_match("gender")), nodes = nodes, group_labels = group_labels)
  nw <- erpm_api_formula_network(call)

  expect_equal(network::network.vertex.names(nw), c(nodes$label, group_labels))
  expect_equal(network::get.vertex.attribute(nw, "gender"), c(nodes$gender, rep(NA, length(group_labels))))
  expect_equal(network::get.vertex.attribute(nw, "age"), c(nodes$age, rep(NA, length(group_labels))))
})

test_that("erpm passes dyadic covariates through partition LHS construction", {
  partition <- erpm_api_partition()
  nodes <- erpm_api_nodes()
  Z1 <- erpm_api_dyad()

  call <- erpm_api_call(partition, rhs = quote(dyadcov("Z1")), nodes = nodes, dyads = list(Z1 = Z1))
  nw <- erpm_api_formula_network(call)
  stored_dyads <- network::get.network.attribute(nw, "dyads")

  expect_named(stored_dyads, "Z1")
  expect_equal(stored_dyads$Z1, Z1, ignore_attr = TRUE)
  expect_equal(rownames(stored_dyads$Z1), nodes$label)
  expect_equal(colnames(stored_dyads$Z1), nodes$label)
})

test_that("erpm normalizes a single dyadic matrix when RHS names one dyadic covariate", {
  partition <- erpm_api_partition()
  Z1 <- erpm_api_dyad()

  call <- erpm_api_call(partition, rhs = quote(dyadcov("Z1")), dyads = Z1)
  nw <- erpm_api_formula_network(call)
  stored_dyads <- network::get.network.attribute(nw, "dyads")

  expect_named(stored_dyads, "Z1")
  expect_equal(stored_dyads$Z1, Z1, ignore_attr = TRUE)
})

test_that("erpm forwards explicit constraints and estimate settings into the ergm call", {
  partition <- erpm_api_partition()
  call <- erpm_api_call(partition, rhs = quote(groups), constraints = ~ b1part, estimate = "MPLE")

  expect_equal(deparse(call$constraints), deparse(~ b1part))
  expect_equal(call$estimate, "MPLE")
})

test_that("erpm validates wrapper inputs before fitting", {
  partition <- erpm_api_partition()

  expect_error(
    erpm("not a formula", eval.call = FALSE),
    "Expected a `lhs ~ ...` formula",
    fixed = TRUE
  )
  expect_error(
    erpm(partition ~ groups, constraints = "not a formula", eval.call = FALSE),
    "`constraints` must be a formula",
    fixed = TRUE
  )
  expect_error(
    erpm(partition ~ groups, seed = 1.5, eval.call = FALSE),
    "`seed` must be integer-valued",
    fixed = TRUE
  )
  #expect_error(
  #  erpm(partition ~ groups, mh_moves = "toggle", eval.call = FALSE),
  #  "`mh_moves` and `mh_weights` must be provided together",
  #  fixed = TRUE
  #)
  #expect_error(
  #  erpm(partition ~ groups, mh_moves = c("toggle", "bad"), mh_weights = c(1, 1), eval.call = FALSE),
  #  "`mh_moves` contains unsupported move",
  #  fixed = TRUE
  #)
})

test_that("erpm public examples are stable", {
  skip("TODO: expand after erpm documentation and vignettes are finalized.")
})
