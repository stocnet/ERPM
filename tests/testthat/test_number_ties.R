
# ISOLATES ----


test_that("same ties isolates", {
  p1 <- 1:6

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)

  sta <- 'sum_pergroup'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 0)

})

test_that("same ties isolates", {
  p1 <- 1:6

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)

  sta <- 'avg_pergroup'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 0)

})


test_that("same ties isolates", {
  p1 <- 1:6

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)

  sta <- 'sum_perind'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 0)

})

test_that("numer ties isolates", {
  p1 <- 1:6

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)

  sta <- 'avg_perind'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 0)

})


# SAME GROUP -----

test_that("numebr ties  same group", {
  p1 <- rep(1,6)

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)
  
  sta <- 'sum_pergroup'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 2)

})

test_that("numebr ties  same group", {
  p1 <- rep(1,6)

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)

  sta <- 'avg_pergroup'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 2)

})


test_that("numebr ties  same group", {
  p1 <- rep(1,6)

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)

  sta <- 'sum_perind'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 4)

})

test_that("numebr ties  same group", {
  p1 <- rep(1,6)

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)

  sta <- 'avg_perind'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 4/6)

})


# RANDOM PARTITION


# SAME GROUP -----

test_that("numebr ties RANDOM", {
  p1 <- c(1,2,2,3,3,4)

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)
  
  sta <- 'sum_pergroup'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 1)

})

test_that("numebr ties RANDOM", {
  p1 <- c(1,2,2,3,3,4)

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)

  sta <- 'avg_pergroup'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 1/4)

})


test_that("number ties RANDOM", {
  p1 <- c(1,2,2,3,3,4)

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)

  sta <- 'sum_perind'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 2)

})

test_that("number ties RANDOM", {
  p1 <- c(1,2,2,3,3,4)

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)

  sta <- 'avg_perind'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 2/6)

})
