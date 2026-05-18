partition_isolates <- 1:6
partition_one_group <- rep(1L, 6)
partition_mixed <- c(1L, 2L, 2L, 3L, 3L, 4L)
partition_unordered_labels <- c(10L, 5L, 5L, 10L, 20L)

attribute_binary <- c(1L, 0L, 1L, 0L, 0L, 1L)
attribute_numeric <- c(1, 1.5, 2, 0, 3, 0)

dyad_binary_6 <- matrix(c(
  0, 0, 0, 0, 0, 0,
  1, 0, 0, 1, 0, 0,
  0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 0, 1,
  0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0
), 6, 6, byrow = TRUE)
