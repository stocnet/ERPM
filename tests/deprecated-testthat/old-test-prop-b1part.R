n <- 20
k <- 20

nw0 <- network::network.initialize(n + k, FALSE, bipartite = k)
nw0[, k + 1] <- 1 # Dump everything into the first partition.

test_that("simulating from null partition distribution", {
  ## Simulate from the null distribution.
  expect_silent(sim0 <- simulate(nw0 ~ b2sociality(nodes = TRUE), coef = numeric(k), constraints = ~b1part, nsim = 100, output = "stats"))

  ### TODO:
  ## apply(sim0, 1, tabulate, nbins = k) obtains the simulate distribution of partition sizes. What should it be?
})
