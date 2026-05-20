### HELPER FUNCTIONS ###

empile_partitions <- function() {
  list(
    c(1, 1, 2, 2),
    c(1, 2, 2, 1),
    c(1, 2, 1, 2),
    c(1, 1, 2, 3)
  )
}

empile_nodes <- function(partitions) {
  lapply(seq_along(partitions), function(t) {
    data.frame(
      label = paste0("actor", seq_along(partitions[[t]]), "_t", t),
      age = seq_along(partitions[[t]]) + 10 * t,
      team = rep(c("A", "B"), length.out = length(partitions[[t]])),
      stringsAsFactors = FALSE
    )
  })
}

empile_dyads <- function(partitions) {
  lapply(seq_along(partitions), function(t) {
    n <- length(partitions[[t]])
    Z <- matrix(t + seq_len(n * n), n, n)
    diag(Z) <- 0
    list(Z1 = Z)
  })
}

build_empile_meta <- function(partitions,
                              rhs = quote(cliques(k = 2)),
                              inertial_present = FALSE,
                              past_influence = 0L,
                              nodes = NULL,
                              dyads = NULL) {
  getFromNamespace(".erpm_long_empile_build_meta_nw", "ERPM")(
    partitions = partitions,
    rhs = rhs,
    inertial_present = inertial_present,
    past_influence = past_influence,
    nodes = nodes,
    dyads = dyads,
    group_labels = NULL,
    directed = FALSE,
    verbose = FALSE,
    debug = FALSE
  )
}

### ERPM_LONG EMPILÉ ENGINE ###

test_that("erpm_long empile engine selects expected time blocks", {
  partitions <- empile_partitions()

  non_inertial <- build_empile_meta(partitions)
  inertial_d1 <- build_empile_meta(
    partitions,
    rhs = quote(inertia_groups(past_influence = 1)),
    inertial_present = TRUE,
    past_influence = 1
  )
  inertial_d2 <- build_empile_meta(
    partitions,
    rhs = quote(inertia_groups(past_influence = 2)),
    inertial_present = TRUE,
    past_influence = 2
  )

  expect_equal(non_inertial$selected_partition_indices, 1:4)
  expect_equal(inertial_d1$selected_partition_indices, 2:4)
  expect_equal(inertial_d2$selected_partition_indices, 3:4)

  expect_equal(network::get.network.attribute(non_inertial$meta_nw, "erpm_long.selected_partition_indices"), 1:4)
  expect_equal(network::get.network.attribute(inertial_d1$meta_nw, "erpm_long.selected_partition_indices"), 2:4)
  expect_equal(network::get.network.attribute(inertial_d2$meta_nw, "erpm_long.selected_partition_indices"), 3:4)
})

test_that("erpm_long empile engine attaches core meta-network attributes", {
  partitions <- empile_partitions()
  built <- build_empile_meta(partitions)
  nw <- built$meta_nw

  expect_s3_class(nw, "network")
  expect_equal(network::get.network.attribute(nw, "erpm_mode"), "empile")
  expect_equal(network::get.network.attribute(nw, "erpm_long.mode"), "PLE")
  expect_equal(network::get.network.attribute(nw, "erpm_long.T"), length(partitions))
  expect_equal(network::get.network.attribute(nw, "erpm_long.d"), 0L)
  expect_equal(network::get.network.attribute(nw, "erpm_long.nbr_actors_meta"), sum(lengths(partitions)))
  expect_equal(network::get.network.attribute(nw, "erpm_long.nbr_actors_by_t"), lengths(partitions))
  expect_equal(network::get.network.attribute(nw, "erpm_long.actor_offsets"), c(0L, 4L, 8L, 12L))

  expect_equal(network::get.network.attribute(nw, "bipartite"), sum(lengths(partitions)))
  expect_equal(network::network.size(nw), 2 * sum(lengths(partitions)))
  expect_equal(network::network.edgecount(nw), sum(lengths(partitions)))
})

test_that("erpm_long empile engine attaches timeblock vertex attributes", {
  partitions <- empile_partitions()
  built <- build_empile_meta(partitions)
  nw <- built$meta_nw

  n_actor <- network::get.network.attribute(nw, "bipartite")
  timeblock <- network::get.vertex.attribute(nw, "timeblock")
  expected_actor_timeblock <- rep(seq_along(partitions), lengths(partitions))
  expected_group_timeblock <- expected_actor_timeblock

  expect_equal(length(timeblock), network::network.size(nw))
  expect_equal(timeblock[seq_len(n_actor)], expected_actor_timeblock)
  expect_equal(timeblock[n_actor + seq_len(n_actor)], expected_group_timeblock)
})

test_that("erpm_long empile engine attaches inertial past partitions", {
  partitions <- empile_partitions()
  built <- build_empile_meta(
    partitions,
    rhs = quote(inertia_groups(past_influence = 2)),
    inertial_present = TRUE,
    past_influence = 2
  )
  nw <- built$meta_nw

  expect_equal(network::get.network.attribute(nw, "erpm_long.d"), 2L)

  past_by_block <- network::get.network.attribute(nw, "erpm_block_past_partitions")
  expect_length(past_by_block, 2)
  expect_equal(past_by_block[[1]][[1]], partitions[[2]])
  expect_equal(past_by_block[[1]][[2]], partitions[[1]])
  expect_equal(past_by_block[[2]][[1]], partitions[[3]])
  expect_equal(past_by_block[[2]][[2]], partitions[[2]])
})

test_that("erpm_long empile engine carries node and dyad data into the meta-network", {
  partitions <- empile_partitions()
  nodes <- empile_nodes(partitions)
  dyads <- empile_dyads(partitions)
  built <- build_empile_meta(partitions, nodes = nodes, dyads = dyads)
  nw <- built$meta_nw

  n_actor <- network::get.network.attribute(nw, "bipartite")
  expected_age <- unlist(lapply(nodes, `[[`, "age"), use.names = FALSE)
  expected_team <- unlist(lapply(nodes, `[[`, "team"), use.names = FALSE)
  expected_label_raw <- unlist(lapply(nodes, `[[`, "label"), use.names = FALSE)

  expect_equal(network::get.network.attribute(nw, "erpm_long.meta_nodes_names"), c("label_raw", "age", "team"))
  expect_equal(network::get.vertex.attribute(nw, "age")[seq_len(n_actor)], expected_age)
  expect_equal(network::get.vertex.attribute(nw, "team")[seq_len(n_actor)], expected_team)
  expect_equal(network::get.vertex.attribute(nw, "label_raw")[seq_len(n_actor)], expected_label_raw)

  stored_dyads <- network::get.network.attribute(nw, "dyads")
  expect_named(stored_dyads, "Z1")
  expect_equal(network::get.network.attribute(nw, "erpm_long.meta_dyads_names"), "Z1")
  expect_equal(network::get.network.attribute(nw, "erpm_long.dyads_mode"), "timeline")
  expect_equal(dim(stored_dyads$Z1), c(n_actor, n_actor))
  expect_equal(stored_dyads$Z1[1:4, 1:4], dyads[[1]]$Z1, ignore_attr = TRUE)
  expect_equal(stored_dyads$Z1[5:8, 5:8], dyads[[2]]$Z1, ignore_attr = TRUE)
  expect_equal(stored_dyads$Z1[1:4, 5:8], matrix(0, 4, 4), ignore_attr = TRUE)
})
