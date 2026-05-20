### HELPER FUNCTIONS ###

valid_longitudinal_partitions <- function() {
  list(
    c(1, 1, 2, 2),
    c(1, 2, 2, 1),
    c(1, 2, 1, 2)
  )
}

valid_longitudinal_nodes <- function(partitions) {
  lapply(seq_along(partitions), function(t) {
    data.frame(
      label = paste0("A", seq_along(partitions[[t]]), "_t", t),
      gender = c("F", "M", "F", "M")[seq_along(partitions[[t]])],
      age = seq_along(partitions[[t]]) + 20 + t,
      stringsAsFactors = FALSE
    )
  })
}

valid_longitudinal_dyads <- function(partitions) {
  lapply(partitions, function(partition) {
    n <- length(partition)
    Z <- matrix(seq_len(n * n), n, n)
    diag(Z) <- 0
    list(Z1 = Z)
  })
}

longitudinal_call <- function(partitions, nodes = NULL, dyads = NULL, mode = "empile", rhs = quote(cliques(k = 2))) {
  formula <- as.formula(call("~", quote(partitions), rhs))
  environment(formula) <- environment()

  erpm_long(
    formula,
    mode = mode,
    eval.call = TRUE,
    verbose = FALSE,
    debug = FALSE,
    nodes = nodes,
    dyads = dyads
  )
}

### ERPM_LONG VALIDATORS ###

test_that("erpm_long validates formula and partition inputs", {
  partitions <- valid_longitudinal_partitions()

  expect_error(
    erpm_long("not a formula"),
    "formula must be a formula",
    fixed = TRUE
  )
  expect_error(
    erpm_long(~ cliques(k = 2)),
    "formula must be of the form: partitions_list ~ terms",
    fixed = TRUE
  )
  expect_error(
    erpm_long(partitions ~ 1),
    "RHS must contain at least one term",
    fixed = TRUE
  )

  bad_partitions <- c(1, 1, 2)
  expect_error(
    erpm_long(bad_partitions ~ cliques(k = 2)),
    "LHS must evaluate to a non-empty list of partitions",
    fixed = TRUE
  )

  one_partition <- list(c(1, 1, 2))
  expect_error(
    erpm_long(one_partition ~ cliques(k = 2)),
    "Only one partition provided",
    fixed = TRUE
  )

  empty_partition <- list(c(1, 1), integer(0))
  expect_error(
    erpm_long(empty_partition ~ cliques(k = 2)),
    "partitions[[2]] is empty",
    fixed = TRUE
  )

  na_partition <- list(c(1, 1), c(1, NA))
  expect_error(
    erpm_long(na_partition ~ cliques(k = 2)),
    "partitions[[2]] contains NA values",
    fixed = TRUE
  )
})

test_that("erpm_long validates mode and simple scalar arguments", {
  partitions <- valid_longitudinal_partitions()

  expect_error(
    longitudinal_call(partitions, mode = "bad"),
    "Invalid mode='bad'",
    fixed = TRUE
  )
  expect_error(
    longitudinal_call(partitions, mode = "sequential"),
    "PLS/sequential mode is not implemented yet",
    fixed = TRUE
  )
  expect_error(
    erpm_long(partitions ~ cliques(k = 2), eval.call = NA),
    "eval.call must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    erpm_long(partitions ~ cliques(k = 2), verbose = NA),
    "verbose must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    erpm_long(partitions ~ cliques(k = 2), debug = NA),
    "debug must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    erpm_long(partitions ~ cliques(k = 2), group_labels = c("G", "H")),
    "group_labels must be NULL or a single character string",
    fixed = TRUE
  )
})

test_that("erpm_long validates nodes inputs", {
  partitions <- valid_longitudinal_partitions()
  nodes <- valid_longitudinal_nodes(partitions)

  expect_true(is.call(longitudinal_call(partitions, nodes = nodes)))

  expect_error(
    longitudinal_call(partitions, nodes = nodes[1:2]),
    "Invalid nodes",
    fixed = TRUE
  )

  bad_nodes <- nodes
  bad_nodes[[1]] <- as.list(bad_nodes[[1]])
  expect_error(
    longitudinal_call(partitions, nodes = bad_nodes),
    "nodes[[1]] must be a data.frame",
    fixed = TRUE
  )

  bad_nodes <- nodes
  bad_nodes[[2]] <- bad_nodes[[2]][-1, , drop = FALSE]
  expect_error(
    longitudinal_call(partitions, nodes = bad_nodes),
    "nodes[[2]] has 3 rows but partitions[[2]] has 4 actors",
    fixed = TRUE
  )

  bad_nodes <- nodes
  bad_nodes[[1]]$label <- NULL
  expect_error(
    longitudinal_call(partitions, nodes = bad_nodes),
    "nodes[[1]] must contain a 'label' column",
    fixed = TRUE
  )

  bad_nodes <- nodes
  bad_nodes[[1]] <- bad_nodes[[1]]["label"]
  expect_error(
    longitudinal_call(partitions, nodes = bad_nodes),
    "nodes[[1]] must have at least one covariate besides 'label'",
    fixed = TRUE
  )

  bad_nodes <- nodes
  bad_nodes[[2]]$extra <- 1
  expect_error(
    longitudinal_call(partitions, nodes = bad_nodes),
    "column schema differs at t=2",
    fixed = TRUE
  )

  bad_nodes <- nodes
  bad_nodes[[3]]$label[1] <- NA
  expect_error(
    longitudinal_call(partitions, nodes = bad_nodes),
    "nodes[[3]]$label contains NA/empty values",
    fixed = TRUE
  )
})

test_that("erpm_long flags reordered node labels across time", {
  skip("Known bug: erpm_long currently matches actors by row position but does not flag reordered labels.")

  partitions <- valid_longitudinal_partitions()
  nodes <- valid_longitudinal_nodes(partitions)
  nodes[[2]] <- nodes[[2]][c(2, 1, 3, 4), , drop = FALSE]

  expect_error(
    longitudinal_call(partitions, nodes = nodes),
    "labels must preserve row-position identity across time",
    fixed = TRUE
  )
})

test_that("erpm_long validates dyads inputs", {
  partitions <- valid_longitudinal_partitions()
  dyads <- valid_longitudinal_dyads(partitions)

  expect_true(is.call(longitudinal_call(partitions, dyads = dyads, rhs = quote(dyadcov("Z1")))))

  expect_error(
    longitudinal_call(partitions, dyads = dyads[1:2]),
    "Invalid dyads",
    fixed = TRUE
  )

  bad_dyads <- dyads
  bad_dyads[[1]] <- list(matrix(0, 4, 4))
  expect_error(
    longitudinal_call(partitions, dyads = bad_dyads),
    "dyads[[1]] must be a NAMED list of matrices",
    fixed = TRUE
  )

  bad_dyads <- dyads
  bad_dyads[[1]]$Z1 <- as.data.frame(bad_dyads[[1]]$Z1)
  expect_error(
    longitudinal_call(partitions, dyads = bad_dyads),
    "dyads[[1]][['Z1']] must be a matrix",
    fixed = TRUE
  )

  bad_dyads <- dyads
  bad_dyads[[2]]$Z1 <- matrix("x", 4, 4)
  expect_error(
    longitudinal_call(partitions, dyads = bad_dyads),
    "dyads[[2]][['Z1']] must be numeric",
    fixed = TRUE
  )

  bad_dyads <- dyads
  bad_dyads[[2]]$Z1 <- matrix(0, 3, 3)
  expect_error(
    longitudinal_call(partitions, dyads = bad_dyads),
    "dyads[[2]][['Z1']] has dim 3x3; expected 4x4",
    fixed = TRUE
  )

  bad_dyads <- dyads
  bad_dyads[[3]]$Z1[1, 2] <- Inf
  expect_error(
    longitudinal_call(partitions, dyads = bad_dyads),
    "dyads[[3]][['Z1']] contains non-finite values",
    fixed = TRUE
  )
})
