
# ISOLATES ----

test_that("Same pairs isolate sum", {
  p1 <- 1:6

  at1 <- rep(1,6)
  at2 <- 1:6
  at3 <- c(1,0,1,0,0,1)

  sta <- 'sum_pergroup'

  expect_equal(same_pairs(p1,at1,sta), 0)
  expect_equal(same_pairs(p1,at2,sta), 0)
  expect_equal(same_pairs(p1,at3,sta), 0)

})

test_that("Same pairs isolate avg", {
  p1 <- 1:6

  at1 <- rep(1,6)
  at2 <- 1:6
  at3 <- c(1,0,1,0,0,1)

  sta <- 'avg_pergroup'

  expect_equal(same_pairs(p1,at1,sta), 0)
  expect_equal(same_pairs(p1,at2,sta), 0)
  expect_equal(same_pairs(p1,at3,sta), 0)

})

# ONE GROUP -----

test_that("Same pairs same group sum", {
  p1 <- rep(1,6)

  at1 <- rep(1,6)
  at2 <- 1:6
  at3 <- c(1,0,1,0,0,1)

  sta <- 'sum_pergroup'

  expect_equal(same_pairs(p1,at1,sta), 15)
  expect_equal(same_pairs(p1,at2,sta), 0)
  expect_equal(same_pairs(p1,at3,sta), 6)

})

test_that("Same pairs same group avg", {
  p1 <- rep(1,6)

  at1 <- rep(1,6)
  at2 <- 1:6
  at3 <- c(1,0,1,0,0,1)

  sta <- 'avg_pergroup'

  expect_equal(same_pairs(p1,at1,sta), 15)
  expect_equal(same_pairs(p1,at2,sta), 0)
  expect_equal(same_pairs(p1,at3,sta), 6)

})


# Random ----

test_that("Same pairs random sum", {
  p1 <- c(1,2,2,3,3,4)

  at1 <- rep(1,6)
  at2 <- 1:6
  at3 <- c(1,0,1,0,0,1)

  sta <- 'sum_pergroup'

  expect_equal(same_pairs(p1,at1,sta), 2)
  expect_equal(same_pairs(p1,at2,sta), 0)
  expect_equal(same_pairs(p1,at3,sta), 1)

})

test_that("Same pairs random avg", {
  p1 <- c(1,2,2,3,3,4)

  at1 <- rep(1,6)
  at2 <- 1:6
  at3 <- c(1,0,1,0,0,1)

  sta <- 'avg_pergroup'

  expect_equal(same_pairs(p1,at1,sta), 2/4)
  expect_equal(same_pairs(p1,at2,sta), 0)
  expect_equal(same_pairs(p1,at3,sta), 1/4)

})


