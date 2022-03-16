
# ISOLATES ----

test_that("similar pairs isolate sum t= 0", {
  p1 <- 1:6

  at1 <- rep(1,6)
  at2 <- 1:6
  at3 <- c(1,0,1,0,0,1)

  sta <- 'sum_pergroup'
  t <- 0
  expect_equal(similar_pairs(p1,at1,sta,t), 0)
  expect_equal(similar_pairs(p1,at2,sta,t), 0)
  expect_equal(similar_pairs(p1,at3,sta,t), 0)

})

test_that("similar pairs isolate sum t= 0", {
  p1 <- 1:6

  at1 <- rep(1,6)
  at2 <- 1:6
  at3 <- c(1,0,1,0,0,1)

  sta <- 'avg_pergroup'
  t <- 0
  expect_equal(similar_pairs(p1,at1,sta,t), 0)
  expect_equal(similar_pairs(p1,at2,sta,t), 0)
  expect_equal(similar_pairs(p1,at3,sta,t), 0)

})

test_that("similar pairs isolate sum t= 8", {
  p1 <- 1:6

  at1 <- rep(1,6)
  at2 <- 1:6
  at3 <- c(1,0,1,0,0,1)

  sta <- 'sum_pergroup'
  t <- 8
  expect_equal(similar_pairs(p1,at1,sta,t), 0)
  expect_equal(similar_pairs(p1,at2,sta,t), 0)
  expect_equal(similar_pairs(p1,at3,sta,t), 0)

})

test_that("similar pairs isolate sum t= 8", {
  p1 <- 1:6

  at1 <- rep(1,6)
  at2 <- 1:6
  at3 <- c(1,0,1,0,0,1)

  sta <- 'avg_pergroup'
  t <- 8
  expect_equal(similar_pairs(p1,at1,sta,t), 0)
  expect_equal(similar_pairs(p1,at2,sta,t), 0)
  expect_equal(similar_pairs(p1,at3,sta,t), 0)

})

test_that("similar pairs isolate sum t= 1", {
  p1 <- 1:6

  at1 <- rep(1,6)
  at2 <- 1:6
  at3 <- c(1,0,1,0,0,1)

  sta <- 'sum_pergroup'
  t <- 1
  expect_equal(similar_pairs(p1,at1,sta,t), 0)
  expect_equal(similar_pairs(p1,at2,sta,t), 0)
  expect_equal(similar_pairs(p1,at3,sta,t), 0)

})

test_that("similar pairs isolate sum t= 1", {
  p1 <- 1:6

  at1 <- rep(1,6)
  at2 <- 1:6
  at3 <- c(1,0,1,0,0,1)

  sta <- 'avg_pergroup'
  t <- 1
  expect_equal(similar_pairs(p1,at1,sta,t), 0)
  expect_equal(similar_pairs(p1,at2,sta,t), 0)
  expect_equal(similar_pairs(p1,at3,sta,t), 0)

})


# ONE GROUP -----

test_that("similar pairs same group sum t= 0", {
  p1 <- rep(1,6)

  at1 <- c(1,1.5,2,0,3,0)
  at2 <- rep(1,6)
  at3 <- c(1,2,3,4,5,6)

  sta <- 'sum_pergroup'
  t <- 0
  expect_equal(similar_pairs(p1,at1,sta,t), 1)
  expect_equal(similar_pairs(p1,at2,sta,t), 15)
  expect_equal(similar_pairs(p1,at3,sta,t), 0)

})

test_that("similar pairs same group avg t= 0", {
  p1 <- rep(1,6)

  at1 <- c(1,1.5,2,0,3,0)
  at2 <- rep(1,6)
  at3 <- c(1,2,3,4,5,6)

  sta <- 'avg_pergroup'
  t <- 0
  expect_equal(similar_pairs(p1,at1,sta,t), 1)
  expect_equal(similar_pairs(p1,at2,sta,t), 15)
  expect_equal(similar_pairs(p1,at3,sta,t), 0)

})



test_that("similar pairs same group sum t= 8", {
  p1 <- rep(1,6)

  at1 <- c(1,1.5,2,0,3,0)
  at2 <- rep(1,6)
  at3 <- c(1,2,3,4,5,6)

  sta <- 'sum_pergroup'
  t <- 8
  expect_equal(similar_pairs(p1,at1,sta,t), 15)
  expect_equal(similar_pairs(p1,at2,sta,t), 15)
  expect_equal(similar_pairs(p1,at3,sta,t), 15)

})

test_that("similar pairs same group avg t= 8", {
  p1 <- rep(1,6)

  at1 <- c(1,1.5,2,0,3,0)
  at2 <- rep(1,6)
  at3 <- c(1,2,3,4,5,6)

  sta <- 'avg_pergroup'
  t <- 8
  expect_equal(similar_pairs(p1,at1,sta,t), 15)
  expect_equal(similar_pairs(p1,at2,sta,t), 15)
  expect_equal(similar_pairs(p1,at3,sta,t), 15)

})



test_that("similar pairs same group sum t= 1", {
  p1 <- rep(1,6)

  at1 <- c(1,1.5,2,0,3,0)
  at2 <- rep(1,6)
  at3 <- c(1,2,3,4,5,6)

  sta <- 'sum_pergroup'
  t <- 1
  expect_equal(similar_pairs(p1,at1,sta,t), 7)
  expect_equal(similar_pairs(p1,at2,sta,t), 15)
  expect_equal(similar_pairs(p1,at3,sta,t), 5)

})

test_that("similar pairs same group avg t= 0", {
  p1 <- rep(1,6)

  at1 <- c(1,1.5,2,0,3,0)
  at2 <- rep(1,6)
  at3 <- c(1,2,3,4,5,6)

  sta <- 'avg_pergroup'
  t <- 1
  expect_equal(similar_pairs(p1,at1,sta,t), 7)
  expect_equal(similar_pairs(p1,at2,sta,t), 15)
  expect_equal(similar_pairs(p1,at3,sta,t), 5)

})


# Random ----

test_that("similar pairs random sum t= 0", {
  p1 <- c(1,2,2,3,3,4)

  at1 <- c(1,1.5,2,0,3,0)
  at2 <- rep(1,6)
  at3 <- c(1,2,3,4,5,6)

  sta <- 'sum_pergroup'
  t <- 0
  expect_equal(similar_pairs(p1,at1,sta,t), 0)
  expect_equal(similar_pairs(p1,at2,sta,t), 2)
  expect_equal(similar_pairs(p1,at3,sta,t), 0)

})

test_that("similar pairs same group avg t= 0", {
  p1 <- c(1,2,2,3,3,4)

  at1 <- c(1,1.5,2,0,3,0)
  at2 <- rep(1,6)
  at3 <- c(1,2,3,4,5,6)

  sta <- 'avg_pergroup'
  t <- 0
  expect_equal(similar_pairs(p1,at1,sta,t), 0)
  expect_equal(similar_pairs(p1,at2,sta,t), 2/4)
  expect_equal(similar_pairs(p1,at3,sta,t), 0)

})



test_that("similar pairs random sum t= 8", {
  p1 <- c(1,2,2,3,3,4)

  at1 <- c(1,1.5,2,0,3,0)
  at2 <- rep(1,6)
  at3 <- c(1,2,3,4,5,6)

  sta <- 'sum_pergroup'
  t <- 8
  expect_equal(similar_pairs(p1,at1,sta,t), 2)
  expect_equal(similar_pairs(p1,at2,sta,t), 2)
  expect_equal(similar_pairs(p1,at3,sta,t), 2)

})

test_that("similar pairs RANDOM avg t= 8", {
  p1 <- c(1,2,2,3,3,4)
  at1 <- c(1,1.5,2,0,3,0)
  at2 <- rep(1,6)
  at3 <- c(1,2,3,4,5,6)

  sta <- 'avg_pergroup'
  t <- 8
  expect_equal(similar_pairs(p1,at1,sta,t), 2/4)
  expect_equal(similar_pairs(p1,at2,sta,t), 2/4)
  expect_equal(similar_pairs(p1,at3,sta,t), 2/4)

})



test_that("similar pairs RANDOM sum t= 1", {
  p1 <- c(1,2,2,3,3,4)
  at1 <- c(1,1.5,2,0,3,0)
  at2 <- rep(1,6)
  at3 <- c(1,2,3,4,5,6)

  sta <- 'sum_pergroup'
  t <- 1
  expect_equal(similar_pairs(p1,at1,sta,t), 1)
  expect_equal(similar_pairs(p1,at2,sta,t), 2)
  expect_equal(similar_pairs(p1,at3,sta,t), 2)

})

test_that("similar pairs RANDOM avg t= 1", {
  p1 <- c(1,2,2,3,3,4)
  at1 <- c(1,1.5,2,0,3,0)
  at2 <- rep(1,6)
  at3 <- c(1,2,3,4,5,6)

  sta <- 'avg_pergroup'
  t <- 1
  expect_equal(similar_pairs(p1,at1,sta,t), 1/4)
  expect_equal(similar_pairs(p1,at2,sta,t), 2/4)
  expect_equal(similar_pairs(p1,at3,sta,t), 2/4)

})

