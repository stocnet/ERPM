test_that("build_bipartite_from_inputs builds a padded bipartite membership network", {
  partition <- c(1L, 1L, 2L, 3L)
  built <- build_bipartite_from_inputs(partition)
  nw <- built$network

  expect_s3_class(nw, "network")
  expect_equal(built$partition, partition)
  expect_equal(built$actor_labels, paste0("A", 1:4))
  expect_equal(built$group_labels, paste0("G", 1:4))

  expect_equal(network::network.size(nw), 8)
  expect_equal(network::network.edgecount(nw), 4)
  expect_equal(network::get.network.attribute(nw, "bipartite"), 4)
  expect_equal(network::network.vertex.names(nw), c(paste0("A", 1:4), paste0("G", 1:4)))

  adjacency <- as.matrix(nw, matrix.type = "adjacency")
  expect_equal(adjacency["A1", "G1"], 1)
  expect_equal(adjacency["A2", "G1"], 1)
  expect_equal(adjacency["A3", "G2"], 1)
  expect_equal(adjacency["A4", "G3"], 1)
  expect_equal(sum(adjacency[paste0("A", 1:4), paste0("G", 1:4)]), 4)
  expect_equal(sum(adjacency), 4)
})

test_that("build_bipartite_from_inputs attaches node data.frame attributes to actor vertices", {
  partition <- c(1L, 1L, 2L, 3L)
  nodes <- data.frame(
    label = c("alice", "bob", "chloe", "david"),
    age = c(20, 21, 30, 40),
    type = c("x", "x", "y", "z"),
    stringsAsFactors = FALSE
  )

  built <- build_bipartite_from_inputs(partition, nodes = nodes)
  nw <- built$network

  expect_equal(built$actor_labels, nodes$label)
  expect_equal(network::network.vertex.names(nw), c(nodes$label, paste0("G", 1:4)))
  expect_equal(network::get.vertex.attribute(nw, "age"), c(nodes$age, rep(NA, 4)))
  expect_equal(network::get.vertex.attribute(nw, "type"), c(nodes$type, rep(NA, 4)))
})

test_that("build_bipartite_from_inputs accepts named list node attributes", {
  partition <- c(1L, 2L, 2L)
  nodes <- list(age = c(20, 30, 40), category = c("a", "b", "b"))

  built <- build_bipartite_from_inputs(partition, nodes = nodes)
  nw <- built$network

  expect_equal(built$actor_labels, paste0("A", 1:3))
  expect_equal(network::network.vertex.names(nw), c(paste0("A", 1:3), paste0("G", 1:3)))
  expect_equal(network::get.vertex.attribute(nw, "age"), c(20, 30, 40, rep(NA, 3)))
  expect_equal(network::get.vertex.attribute(nw, "category"), c("a", "b", "b", rep(NA, 3)))
})

test_that("build_bipartite_from_inputs stores dyads with actor label dimnames", {
  partition <- c(1L, 1L, 2L)
  nodes <- data.frame(label = c("a", "b", "c"), stringsAsFactors = FALSE)
  friendship <- matrix(c(
    0, 1, 0,
    1, 0, 1,
    0, 1, 0
  ), 3, 3, byrow = TRUE)

  built <- build_bipartite_from_inputs(partition, nodes = nodes, dyads = list(friendship = friendship))
  stored_dyads <- network::get.network.attribute(built$network, "dyads")

  expect_named(stored_dyads, "friendship")
  expect_equal(stored_dyads$friendship, friendship, ignore_attr = TRUE)
  expect_equal(rownames(stored_dyads$friendship), nodes$label)
  expect_equal(colnames(stored_dyads$friendship), nodes$label)
})


test_that("build_bipartite_from_inputs rejects invalid inputs", {
  expect_error(build_bipartite_from_inputs(NULL), "partition must be a non-empty atomic vector", fixed = TRUE)
  expect_error(build_bipartite_from_inputs(c(1, 0, 2)), "partition must contain finite positive integers", fixed = TRUE)
  expect_error(build_bipartite_from_inputs(c(1, 1.5, 2)), "partition must contain integer-valued group ids", fixed = TRUE)
  expect_error(build_bipartite_from_inputs(c(1, 4, 1)), "max(partition) cannot exceed n", fixed = TRUE)

  bad_nodes <- data.frame(label = c("a", "b"), stringsAsFactors = FALSE)
  expect_error(build_bipartite_from_inputs(c(1, 1, 2), nodes = bad_nodes), "nrow(nodes) must equal length(partition)", fixed = TRUE)

  duplicate_labels <- data.frame(label = c("a", "a", "b"), stringsAsFactors = FALSE)
  expect_error(build_bipartite_from_inputs(c(1, 1, 2), nodes = duplicate_labels), "Invalid `nodes` data.frame: nodes$label contains duplicates.", fixed = TRUE)

  bad_dyad <- matrix(0, 2, 2)
  expect_error(build_bipartite_from_inputs(c(1, 1, 2), dyads = list(friendship = bad_dyad)), "Invalid `dyads` list", fixed = TRUE)

  expect_error(build_bipartite_from_inputs(c(1, 1, 2), group_labels = c("g1", "g2")), "group_labels must be an atomic vector of length G", fixed = TRUE)
  expect_error(build_bipartite_from_inputs(c(1, 1, 2), group_labels = c("g1", "g1", "g3")), "duplicate group_labels are not allowed", fixed = TRUE)
})
