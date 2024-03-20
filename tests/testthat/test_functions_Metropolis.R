##### Neighborhoods size tests: special cases ----
n <- 5

# extreme case: all isolates
partition <- 1:n

test_that("NS isolates p1", {
  s1 <- compute_size_neighborhood_p1(partition)
  s2 <- compute_size_neighborhood_p1_restricted(partition,1:n,1:n)

  expect_equal(s1$nums.swaps , matrix(0,n,n))
  expect_equal(s1$num.swaps , 0)

  expect_equal(s2$nums.swaps , matrix(0,n,n))
  expect_equal(s2$num.swaps , 0)
})



test_that("NS isolates p2", {
  s1 <- compute_size_neighborhood_p2(partition)
  s2 <- compute_size_neighborhood_p2_restricted(partition,1:n,1:n)

  expect_equal(s1$num.merges , 10)

  expect_equal(s2$merges , matrix(c(0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,1,1,1,0),nrow=5,ncol=5))
  expect_equal(s2$num.merges ,10)

  expect_equal(s1$nums.divisions,c(0,0,0,0,0))
  expect_equal(s1$num.divisions,0)

  expect_equal(s2$nums.divisions,c(0,0,0,0,0))
  expect_equal(s2$num.divisions,0)
})

test_that("Neigh Size isolates p3", {
  s1 <- compute_size_neighborhood_p3(partition)
  s2 <- compute_size_neighborhood_p3_restricted(partition,1:n,1:n)

  expect_equal(s1$nums.swaps , c(4,3,2,1,0))
  expect_equal(s1$num.swaps,10)

  expect_equal(s2$nums.swaps , c(4,3,2,1,0))
  expect_equal(s2$num.swaps,10)
})


# extreme case: all in the same group
partition <- rep(1,n)

test_that("NS same group  p1", {
  s1 <- compute_size_neighborhood_p1(partition)
  s2 <- compute_size_neighborhood_p1_restricted(partition,1:n,1:n)

  expect_equal(s1$nums.swaps , matrix(0,n,n))
  expect_equal(s1$num.swaps , 0)

  expect_equal(s2$nums.swaps , matrix(0,n,n))
  expect_equal(s2$num.swaps , 0)
})



test_that("NS same group p2", {
  s1 <- compute_size_neighborhood_p2(partition)
  s2 <- compute_size_neighborhood_p2_restricted(partition,1:n,1:n)

  expect_equal(s1$num.merges , 0)

  expect_equal(s2$merges,matrix(0,1,1))
  expect_equal(s2$num.merges ,0)

  expect_equal(s1$nums.divisions,15)
  expect_equal(s1$num.divisions, 15)

  expect_equal(s2$nums.divisions, 15)
  expect_equal(s2$num.divisions,15)
})

test_that("Neigh Size same group p3", {
  s1 <- compute_size_neighborhood_p3(partition)
  s2 <- compute_size_neighborhood_p3_restricted(partition,1:n,1:n)

  expect_equal(s1$nums.swaps , c(1,1,1,1,1))
  expect_equal(s1$num.swaps,5)

  expect_equal(s2$nums.swaps , c(1,1,1,1,1))
  expect_equal(s2$num.swaps,5)
})


# random case: c(1,1,2,2,3)
partition <- c(1,1,2,2,3)

test_that("NS Random case  p1", {
  s1 <- compute_size_neighborhood_p1(partition)
  s2 <- compute_size_neighborhood_p1_restricted(partition,1:n,1:n)

  expect_equal(s1$nums.swaps , matrix(c(0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,1,1,1,0),nrow=5,ncol=5))
  expect_equal(s1$num.swaps , 8)

  expect_equal(s2$nums.swaps , matrix(c(0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,1,1,1,0),nrow=5,ncol=5))
  expect_equal(s2$num.swaps , 8)
})

test_that("NS same group p2", {
  s1 <- compute_size_neighborhood_p2(partition)
  s2 <- compute_size_neighborhood_p2_restricted(partition,1:n,1:n)

  expect_equal(s1$num.merges , 3)

  expect_equal(s2$merges,matrix(c(0,0,0,1,0,0,1,1,0),nrow=3,ncol=3))
  expect_equal(s2$num.merges ,3)

  expect_equal(s1$nums.divisions,c(1,1,0))
  expect_equal(s1$num.divisions, 2)

  expect_equal(s2$nums.divisions, c(1,1,0))
  expect_equal(s2$num.divisions,2)
})

test_that("NS random case p3", {
  s1 <- compute_size_neighborhood_p3(partition)
  s2 <- compute_size_neighborhood_p3_restricted(partition,1:n,1:n)

  expect_equal(s1$nums.swaps , c(3,2,3,2,2))
  expect_equal(s1$num.swaps,12)

  expect_equal(s2$nums.swaps , c(3,2,3,2,2))
  expect_equal(s2$num.swaps,12)
})


# random case: c(1,1,1,2,2)
partition <- c(1,1,1,2,2)

test_that("NS Random case  p1", {
  s1 <- compute_size_neighborhood_p1(partition)
  s2 <- compute_size_neighborhood_p1_restricted(partition,1:n,1:n)

  expect_equal(s1$nums.swaps , matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,0),nrow=5,ncol=5))
  expect_equal(s1$num.swaps , 6)

  expect_equal(s2$nums.swaps , matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,0),nrow=5,ncol=5))
  expect_equal(s2$num.swaps , 6)
})

test_that("NS same group p2", {
  s1 <- compute_size_neighborhood_p2(partition)
  s2 <- compute_size_neighborhood_p2_restricted(partition,1:n,1:n)

  expect_equal(s1$num.merges , 1)

  expect_equal(s2$merges,matrix(c(0,0,1,0),nrow=2,ncol=2))
  expect_equal(s2$num.merges ,1)

  expect_equal(s1$nums.divisions,c(3,1))
  expect_equal(s1$num.divisions, 4)

  expect_equal(s2$nums.divisions, c(3,1))
  expect_equal(s2$num.divisions,4)
})

test_that("NS random case p3", {
  s1 <- compute_size_neighborhood_p3(partition)
  s2 <- compute_size_neighborhood_p3_restricted(partition,1:n,1:n)

  expect_equal(s1$nums.swaps , c(2,2,2,2,1))
  expect_equal(s1$num.swaps, 9)

  expect_equal(s2$nums.swaps , c(2,2,2,2,1))
  expect_equal(s2$num.swaps, 9)
})

##### Neighborhoods size tests: check restricted sizes ----

# random case: c(1,1,2,2,3) and restricted sizes 1 to 3
partition <- c(1,1,2,2,3)

test_that('NS restricted sizes (3) p1',{

  s2 <- compute_size_neighborhood_p1_restricted(partition,2:6,1:3)
  expect_equal(s2$nums.swaps , matrix(c(0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,1,1,1,0),nrow=5,ncol=5))
  expect_equal(s2$num.swaps , 8)

})

test_that("NS restricted sizes (3) p2", {
  s2 <- compute_size_neighborhood_p2_restricted(partition,2:6,1:3)
  expect_equal(s2$merges,matrix(c(0,0,0,0,0,0,1,1,0),nrow=3,ncol=3))
  expect_equal(s2$num.merges ,2)

  expect_equal(s2$nums.divisions, c(1,1,0))
  expect_equal(s2$num.divisions, 2)

})

test_that("NS random case p3", {
  s2 <- compute_size_neighborhood_p3_restricted(partition,2:6,1:3)
  expect_equal(s2$nums.swaps , c(3,2,3,2,2))
  expect_equal(s2$num.swaps, 12)

})


# random case: c(1,1,1,2,2) with restricted sizes 2 to 4
partition <- c(1,1,1,2,2)

test_that('NS restricted sizes (2-4) p1',{
  s2 <- compute_size_neighborhood_p1_restricted(partition,2:3,2:4)
  expect_equal(s2$nums.swaps , matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,0),nrow=5,ncol=5))
  expect_equal(s2$num.swaps , 6)

})

test_that("NS restricted sizes (3) p2", {
  s2 <- compute_size_neighborhood_p2_restricted(partition,2:3,2:4)
  expect_equal(s2$merges,matrix(c(0,0,0,0),nrow=2,ncol=2))
  expect_equal(s2$num.merges ,0)

  expect_equal(s2$nums.divisions, c(0,0))
  expect_equal(s2$num.divisions, 0)

})

test_that("NS random case p3", {
  s2 <- compute_size_neighborhood_p3_restricted(partition,2:3,2:4)
  expect_equal(s2$nums.swaps ,c(1,1,1,0,0))
  expect_equal(s2$num.swaps, 3)

})


##### Neighborhoods size tests: same results with restricted on all sizes -----------
n <- 5
S <- 100

test_that('NS restricted', {
  for(i in 1:S){

  partition <- sample(1:n,n,replace=T)
  partition <- order_groupids(partition)

  s1 <- compute_size_neighborhood_p1(partition)
  s2 <- compute_size_neighborhood_p1_restricted(partition,1:n,1:n)
  expect_equal(s1$nums.swaps,s2$nums.swaps)
  expect_equal(s1$num.swaps,s2$num.swaps)

  s1 <- compute_size_neighborhood_p2(partition)
  s2 <- compute_size_neighborhood_p2_restricted(partition,1:n,1:n)
  expect_equal(s1$num.merges,s2$num.merges)
  expect_equal(s1$num.divisions,s2$num.divisions)
  expect_equal(s1$nums.divisions,s2$nums.divisions)

  s1 <- compute_size_neighborhood_p3(partition)
  s2 <- compute_size_neighborhood_p3_restricted(partition,1:n,1:n)

  expect_equal(s1$nums.swaps,s2$nums.swaps)
  expect_equal(s1$num.swaps,s2$num.swaps)

}
})


##### Neighborhoods sample tests: check some samplings of certain partitions ---


# extreme case: all isolates
test_that('NS sampling : isolates p1', {
  S <- 500
  partition <- 1:n
  s1 <- compute_size_neighborhood_p1(partition)
  s2 <- compute_size_neighborhood_p1_restricted(partition,1:n,1:n)
  allsamples1 <- matrix(0,S,n)
  allsamples2 <- matrix(0,S,n)
  for(i in 1:S) {
    allsamples1[i,] <- sample_new_partition_p1(partition,s1)
    allsamples2[i,] <- sample_new_partition_p1_restricted(partition,s2,1:n,1:n)
  }

  for(i in 1:S){
    expect_equal(1:n,allsamples1[i,])
    expect_equal(1:n,allsamples2[i,])
  }

})


# Redo this test: right now it doesn't work because there are two many options possible
# test_that('NS sampling : isolates p2',{
#   S <- 500
#   partition <- 1:n
#   s1 <- compute_size_neighborhood_p2(partition)
#   s2 <- compute_size_neighborhood_p2_restricted(partition,1:n,1:n)
#   allsamples1 <- matrix(0,S,n)
#   allsamples2 <- matrix(0,S,n)
#   for(i in 1:S) {
#     allsamples1[i,] <- sample_new_partition_p2(partition,s1)
#     allsamples2[i,] <- sample_new_partition_p2_restricted(partition,s2,1:n,1:n)
#   }
# 
#   for(i in 1:S){
#     expect_equal(allsamples1[i,],c(1,2,1,3,4))
#     expect_equal(allsamples2[i,],c(1,2,1,3,4))
#   }
# })

# Redo this test: right now it doesn't work because there are two many options possible
# test_that('NS sampling : isolates p3',{
#   S <- 500
#   partition <- 1:n
#   s1 <- compute_size_neighborhood_p3(partition)
#   s2 <- compute_size_neighborhood_p3_restricted(partition,1:n,1:n)
#   allsamples1 <- matrix(0,S,n)
#   allsamples2 <- matrix(0,S,n)
#   for(i in 1:S) {
#     allsamples1[i,] <- sample_new_partition_p3(partition,s1)
#     allsamples2[i,] <- sample_new_partition_p3_restricted(partition,s2,1:n,1:n)
#   }
# 
#   for(i in 1:S){
#     expect_equal(allsamples1[i,],c(1,2,3,4,1))
#     expect_equal(allsamples2[i,],c(1,2,3,4,1))
#   }
# })


# extreme case: all in the same group

test_that('NS sampling : same group p1', {
  S <- 500
  partition <- rep(1,n)
  s1 <- compute_size_neighborhood_p1(partition)
  s2 <- compute_size_neighborhood_p1_restricted(partition,1:n,1:n)
  allsamples1 <- matrix(0,S,n)
  allsamples2 <- matrix(0,S,n)
  for(i in 1:S) {
    allsamples1[i,] <- sample_new_partition_p1(partition,s1)
    allsamples2[i,] <- sample_new_partition_p1_restricted(partition,s2,1:n,1:n)
  }

  for(i in 1:S){
    expect_equal(rep(1,n),allsamples1[i,])
    expect_equal(rep(1,n),allsamples2[i,])
  }

})

# Redo this test: right now it doesn't work because there are two many options possible
# test_that('NS sampling : same group p2',{
#   S <- 500
#   partition <- rep(1,n)
#   s1 <- compute_size_neighborhood_p2(partition)
#   s2 <- compute_size_neighborhood_p2_restricted(partition,1:n,1:n)
#   allsamples1 <- matrix(0,S,n)
#   allsamples2 <- matrix(0,S,n)
#   for(i in 1:S) {
#     allsamples1[i,] <- sample_new_partition_p2(partition,s1)
#     allsamples2[i,] <- sample_new_partition_p2_restricted(partition,s2,1:n,1:n)
#   }
# 
#   for(i in 1:S){
#     expect_equal(allsamples1[i,],c(1,1,2,2,2))
#     expect_equal(allsamples2[i,],c(1,1,2,2,2))
#   }
# })

# Redo this test: right now it doesn't work because there are two many options possible
# test_that('NS sampling : same group p3',{
#   S <- 500
#   partition <- rep(1,n)
#   s1 <- compute_size_neighborhood_p3(partition)
#   s2 <- compute_size_neighborhood_p3_restricted(partition,1:n,1:n)
#   allsamples1 <- matrix(0,S,n)
#   allsamples2 <- matrix(0,S,n)
#   for(i in 1:S) {
#     allsamples1[i,] <- sample_new_partition_p3(partition,s1)
#     allsamples2[i,] <- sample_new_partition_p3_restricted(partition,s2,1:n,1:n)
#   }
# 
#   for(i in 1:S){
#     expect_equal(allsamples1[i,],c(1,2,2,2,2))
#     expect_equal(allsamples2[i,],c(1,2,2,2,2))
#   }
# })


# random case: c(1,1,2,2,3)
partition <- c(1,1,2,2,3)


# Redo this test: right now it doesn't work because there are two many options possible
# test_that('NS sampling : random case p1', {
#   S <- 500
#   partition <- c(1,1,2,2,3)
#   s1 <- compute_size_neighborhood_p1(partition)
#   s2 <- compute_size_neighborhood_p1_restricted(partition,1:n,1:n)
#   allsamples1 <- matrix(0,S,n)
#   allsamples2 <- matrix(0,S,n)
# 
#   for(i in 1:S) {
#     allsamples1[i,] <- sample_new_partition_p1(partition,s1)
#     allsamples2[i,] <- sample_new_partition_p1_restricted(partition,s2,1:n,1:n)
#   }
# 
#   for(i in 1:S){
#     expect_true(all(c(1,2,1,2,3) == allsamples1[i,]))
#     expect_true(all(c(1,2,1,2,3) == allsamples2[i,]))
#   }
# 
# })

test_that('NS sampling : random case p2',{
  S <- 500
  partition <- c(1,1,2,2,3)
  s1 <- compute_size_neighborhood_p2(partition)
  s2 <- compute_size_neighborhood_p2_restricted(partition,1:n,1:n)
  allsamples1 <- matrix(0,S,n)
  allsamples2 <- matrix(0,S,n)
  for(i in 1:S) {
    allsamples1[i,] <- sample_new_partition_p2(partition,s1)
    allsamples2[i,] <- sample_new_partition_p2_restricted(partition,s2,1:n,1:n)
  }
  found1 <- F
  found2 <- F
  for(i in 1:S){
    if(all(c(1,1,2,2,2) == allsamples1[i,])) { found1 = T}
    if(all(c(1,1,2,2,2) == allsamples2[i,])) { found2 = T}
  }
  expect_true(found1)
  expect_true(found2)
})


# Redo this test: right now it doesn't work because there are two many options possible
# test_that('NS sampling : random case p3',{
#   S <- 500
#   partition <- c(1,1,2,2,3)
#   s1 <- compute_size_neighborhood_p3(partition)
#   s2 <- compute_size_neighborhood_p3_restricted(partition,1:n,1:n)
#   allsamples1 <- matrix(0,S,n)
#   allsamples2 <- matrix(0,S,n)
#   for(i in 1:S) {
#     allsamples1[i,] <- sample_new_partition_p3(partition,s1)
#     allsamples2[i,] <- sample_new_partition_p3_restricted(partition,s2,1:n,1:n)
#   }
# 
#   for(i in 1:S){
#     expect_true(all(c(1,1,2,3,3) == allsamples1[i,]))
#     expect_true(all(c(1,1,2,3,3) == allsamples2[i,]))
#   }
# })



##### Neighborhoods reachable tests: check that reachability is correctly tested
# for some extreme cases

test_that('reachable random',{

  partition1 <- c(1,2,3,4,5)
  partition2 <- c(1,2,3,4,4)
  expect_false(reachable_p1(partition1,partition2))

  partition1 <- c(1,1,1,1,1)
  partition2 <- c(1,2,1,1,1)
  expect_false(reachable_p1(partition1,partition2))

  partition1 <- c(1,2,3,4,5)
  partition2 <- c(1,2,3,4,4)
  expect_true(reachable_p2(partition1,partition2))

  partition1 <- c(1,1,1,1,1)
  partition2 <- c(1,2,1,1,1)
  expect_true(reachable_p2(partition1,partition2))

  partition1 <- c(1,2,3,4,5)
  partition2 <- c(1,2,3,4,4)
  expect_true(reachable_p3(partition1,partition2))

  partition1 <- c(1,1,1,1,1)
  partition2 <- c(1,2,1,1,1)
  expect_true(reachable_p3(partition1,partition2))

})


###################################################################################

##### Neighborhoods reachable tests: check that reachability is correctly tested
# for randomly sampled cases

test_that('Reachable random sampled cases',{

  for(i in 1:10){

    if(i == 1) partition <- 1:n
    if(i == 2) partition <- rep(1,n)
    if(i >= 3) partition <- order_groupids(sample(n,n,replace=T))
    S <- 30

    # Phase 1

    s1 <- compute_size_neighborhood_p1(partition)
    s2 <- compute_size_neighborhood_p1_restricted(partition,1:n,1:n)
    allsamples1 <- matrix(0,S,n)
    allsamples2 <- matrix(0,S,n)
    for(i in 1:S) {
      sample1 <- sample_new_partition_p1(partition,s1)
      sample2 <- sample_new_partition_p1_restricted(partition,s2,1:n,1:n)
      if(s1$total != 0 ) {expect_true(reachable_p1(partition,sample1))}
      if(s2$total != 0 ) { expect_true(reachable_p1(partition,sample2))}
    }

    # Phase 2

    s1 <- compute_size_neighborhood_p2(partition)
    s2 <- compute_size_neighborhood_p2_restricted(partition,1:n,1:n)
    allsamples1 <- matrix(0,S,n)
    allsamples2 <- matrix(0,S,n)
    for(i in 1:S) {
      sample1 <- sample_new_partition_p2(partition,s1)
      sample2 <- sample_new_partition_p2_restricted(partition,s2,1:n,1:n)
      expect_true(reachable_p2(partition,sample1))
      expect_true(reachable_p2(partition,sample2))
    }

    # Phase 3

    s1 <- compute_size_neighborhood_p3(partition)
    s2 <- compute_size_neighborhood_p3_restricted(partition,1:n,1:n)
    allsamples1 <- matrix(0,S,n)
    allsamples2 <- matrix(0,S,n)
    for(i in 1:S) {
      sample1 <- sample_new_partition_p3(partition,s1)
      sample2 <- sample_new_partition_p3_restricted(partition,s2,1:n,1:n)
      expect_true(reachable_p3(partition,sample1))
      expect_true(reachable_p3(partition,sample2))
      }
  }

})


