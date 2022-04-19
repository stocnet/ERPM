
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
# PB pour at3

partition <- c(1,2,2,3,3,4)

dyadic_attribute = at3
dist<- as.matrix(dist(partition, diag = T, upper = T))
net_partition <- dist
net_partition[dist == 0] <- 1
net_partition[dist != 0] <- 0
net_partition[is.na(dist)] <- 0
diag(net_partition) <- 0

numties <- sum(dyadic_attribute * net_partition) / 2
avgroupnumties <- sum(dyadic_attribute * net_partition) / (2*max(partition, na.rm = T))
indnumties <- sum( rowSums(dyadic_attribute * net_partition) > 0 )
avgind <- indnumties/length(partition)


test_that("same ties isolates", {
  p1 <- 1:6

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)

  sta <- 'sum_perind'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 4)

})
#SAME PB
test_that("numer ties isolates", {
  p1 <- 1:6

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)

  sta <- 'avg_perind'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 4/6)

})


# SAME GROUP -----

test_that("numebr ties  same group", {
  p1 <- rep(1,6)

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,a),6,6, byrow = T)

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
  at3 <- matrix(c(0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,a),6,6, byrow = T)

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

# 2 au lieu de 4 , definition de sum_in , sum indiv , par group ou dans toute la partition?

test_that("numebr ties RANDOM", {
  p1 <- c(1,2,2,3,3,4)

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)

  sta <- 'sum_perind'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 4)

})

test_that("numebr ties RANDOM", {
  p1 <- c(1,2,2,3,3,4)

  at1 <- diag(1,6)
  at2 <- diag(0,6)
  a <- rep(0,12)
  at3 <- matrix(c(0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,a,1,0,0,0,0,0),6,6, byrow = T)

  sta <- 'avg_perind'
  expect_equal(number_ties(p1,at1,sta), 0)
  expect_equal(number_ties(p1,at2,sta), 0)
  expect_equal(number_ties(p1,at3,sta), 4/6)

})
