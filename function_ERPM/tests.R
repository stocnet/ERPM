##### Neighborhoods size tests: special cases
n <- 5

# extreme case: all isolates
partition <- 1:n

s1 <- compute_size_neighborhood_p1(partition)
s2 <- compute_size_neighborhood_p1_restricted(partition,1:n)
if(!all(s1$nums.swaps == matrix(0,n,n)) || !(s1$num.swaps == 0)){
  print("problem on swaps with p1 with all isolates")
}
if(!all(s2$nums.swaps == matrix(0,n,n)) || !(s2$num.swaps == 0)){
  print("problem on swaps with p1 restricted with all isolates")
}

s1 <- compute_size_neighborhood_p2(partition)
s2 <- compute_size_neighborhood_p2_restricted(partition,1:n)
if(!(s1$num.merges == 10)){
  print("problem on merges with p2 with all isolates")
}
if(!all(s2$merges == matrix(c(0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,1,1,1,0),nrow=5,ncol=5)) || !(s2$num.merges == 10)){
  print("problem on merges with p2 restricted with all isolates")
}
if(!all(s1$nums.divisions == c(0,0,0,0,0)) || !(s1$num.divisions == 0)){
  print("problem on divisions with p2 with all isolates")
}
if(!all(s2$nums.divisions == c(0,0,0,0,0)) || !(s2$num.divisions == 0)){
  print("problem on divisions with p2 restricted with all isolates")
}

s1 <- compute_size_neighborhood_p3(partition)
s2 <- compute_size_neighborhood_p3_restricted(partition,1:n)
if(!all(s1$nums.swaps == c(4,3,2,1,0)) || !(s1$num.swaps == 10)){
  print("problem on swaps with p3 with all isolates")
}
if(!all(s2$nums.swaps == c(4,3,2,1,0)) || !(s2$num.swaps == 10)){
  print("problem on swaps with p3 restricted with all isolates")
}
  
# s1 <- compute_size_neighborhood_p4(partition)
# s2 <- compute_size_neighborhood_p4_restricted(partition,1:n)
# if(!all(s1$nums.swaps == c(0,0,0,0,0)) || !(s1$num.swaps == 0)){
#   print("problem on swaps with p4 with all isolates")
# }
# if(!all(s2$nums.swaps == c(0,0,0,0,0)) || !(s2$num.swaps == 0)){
#   print("problem on swaps with p4 restricted with all isolates")
# }
# 
# s1 <- compute_size_neighborhood_p5(partition)
# s2 <- compute_size_neighborhood_p5_restricted(partition,1:n)
# if(!(s1$num.swaps == 0)){
#   print("problem on swaps with p5 with all isolates")
# }
# if(!(s2$num.swaps == 0)){
#   print("problem on swaps with p5 restricted with all isolates")
# }
# 
# s1 <- compute_size_neighborhood_p6(partition)
# s2 <- compute_size_neighborhood_p6_restricted(partition,1:n)
# if(!all(s1$nums.swaps == matrix(0,5,5)) || !(s1$num.swaps == 0)){
#   print("problem on swaps with p6 with all isolates")
# }
# if(!all(s2$nums.swaps == matrix(0,5,5)) || !(s2$num.swaps == 0)){
#   print("problem on swaps with p6 restricted with all isolates")
# }
# 
# partition2 <- 1:n
# s1 <- compute_size_neighborhood_p7(partition, partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition, partition2, 1:n)
# if(!all(s1$nums.swaps == matrix(0,5,5)) || !(s1$num.swaps == 0)){
#   print("problem on swaps with p7 with all isolates & all isolates")
# }
# if(!all(s2$nums.swaps == matrix(0,5,5)) || !(s2$num.swaps == 0)){
#   print("problem on swaps with p7 restricted with all isolates & all isolates")
# }
# 
# partition2 <- rep(1,n)
# s1 <- compute_size_neighborhood_p7(partition, partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition, partition2, 1:n)
# if(!all(s1$nums.swaps == matrix(c(0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,1,1,1,0),nrow=5,ncol=5)) || !(s1$num.swaps == 10)){
#   print("problem on swaps with p7 with all isolates & all in the same group")
# }
# if(!all(s2$nums.swaps == matrix(c(0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,1,1,1,0),nrow=5,ncol=5)) || !(s2$num.swaps == 10)){
#   print("problem on swaps with p7 restricted with all isolates & all in the same group")
# }


# extreme case: all in the same group
partition <- rep(1,n)

s1 <- compute_size_neighborhood_p1(partition)
s2 <- compute_size_neighborhood_p1_restricted(partition,1:n)
if(!all(s1$nums.swaps == matrix(0,n,n)) || !(s1$num.swaps == 0)){
  print("problem on swaps with p1 with all in the same group")
}
if(!all(s2$nums.swaps == matrix(0,n,n)) || !(s2$num.swaps == 0)){
  print("problem on swaps with p1 restricted with all in the same group")
}

s1 <- compute_size_neighborhood_p2(partition)
s2 <- compute_size_neighborhood_p2_restricted(partition,1:n)
if(!(s1$num.merges == 0)){
  print("problem on merges with p2 with all in the same group")
}
if(!all(s2$merges == matrix(0,1,1)) || !(s2$num.merges == 0)){
  print("problem on merges with p2 restricted with all in the same group")
}
if(!all(s1$nums.divisions == 15) || !(s1$num.divisions == 15)){
  print("problem on divisions with p2 with all in the same group")
}
if(!all(s2$nums.divisions == 15) || !(s2$num.divisions == 15)){
  print("problem on divisions with p2 restricted with all in the same group")
}

s1 <- compute_size_neighborhood_p3(partition)
s2 <- compute_size_neighborhood_p3_restricted(partition,1:n)
if(!all(s1$nums.swaps == c(1,1,1,1,1)) || !(s1$num.swaps == 5)){
  print("problem on swaps with p3 with all in the same group")
}
if(!all(s2$nums.swaps == c(1,1,1,1,1)) || !(s2$num.swaps == 5)){
  print("problem on swaps with p3 restricted with all in the same group")
}

# s1 <- compute_size_neighborhood_p4(partition)
# s2 <- compute_size_neighborhood_p4_restricted(partition,1:n)
# if(!(s1$nums.swaps == 10) || !(s1$num.swaps == 10)){
#   print("problem on swaps with p4 with all in the same group")
# }
# if(!(s2$nums.swaps == 10) || !(s2$num.swaps == 10)){
#   print("problem on swaps with p4 restricted with all in the same group")
# }
# 
# s1 <- compute_size_neighborhood_p5(partition)
# s2 <- compute_size_neighborhood_p5_restricted(partition,1:n)
# if(!(s1$num.swaps == 0)){
#   print("problem on swaps with p5 with all in the same group")
# }
# if(!(s2$num.swaps == 0)){
#   print("problem on swaps with p5 restricted with all in the same group")
# }
# 
# s1 <- compute_size_neighborhood_p6(partition)
# s2 <- compute_size_neighborhood_p6_restricted(partition,1:n)
# if(!all(s1$nums.swaps == matrix(0,1,1)) || !(s1$num.swaps == 0)){
#   print("problem on swaps with p6 with all in the same group")
# }
# if(!all(s2$nums.swaps == matrix(0,1,1)) || !(s2$num.swaps == 0)){
#   print("problem on swaps with p6 restricted with all in the same group")
# }
# 
# partition2 <- 1:n
# s1 <- compute_size_neighborhood_p7(partition, partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition, partition2, 1:n)
# if(!all(s1$nums.swaps == matrix(0,5,5)) || !(s1$num.swaps == 0)){
#   print("problem on swaps with p7 with all in the same group & all isolates")
# }
# if(!all(s2$nums.swaps == matrix(0,5,5)) || !(s2$num.swaps == 0)){
#   print("problem on swaps with p7 restricted with all in the same group & all isolates")
# }
# 
# partition2 <- rep(1,n)
# s1 <- compute_size_neighborhood_p7(partition, partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition, partition2, 1:n)
# if(!all(s1$nums.swaps == matrix(c(0,1,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,0),nrow=5,ncol=5)) || !(s1$num.swaps == 20)){
#   print("problem on swaps with p7 with all in the same group & all in the same group")
# }
# if(!all(s2$nums.swaps == matrix(c(0,1,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,0),nrow=5,ncol=5)) || !(s2$num.swaps == 20)){
#   print("problem on swaps with p7 restricted with all in the same group & all in the same group")
# }


# random case: c(1,1,2,2,3)
partition <- c(1,1,2,2,3)

s1 <- compute_size_neighborhood_p1(partition)
s2 <- compute_size_neighborhood_p1_restricted(partition,1:n)
if(!all(s1$nums.swaps == matrix(c(0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,1,1,1,0),nrow=5,ncol=5)) || !(s1$num.swaps == 8)){
  print("problem on swaps with p1 with (1,1,2,2,3)")
}
if(!all(s2$nums.swaps == matrix(c(0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,1,1,1,0),nrow=5,ncol=5)) || !(s2$num.swaps == 8)){
  print("problem on swaps with p1 restricted with (1,1,2,2,3)")
}

s1 <- compute_size_neighborhood_p2(partition)
s2 <- compute_size_neighborhood_p2_restricted(partition,1:n)
if(!(s1$num.merges == 3)){
  print("problem on merges with p2 with (1,1,2,2,3)")
}
if(!all(s2$merges == matrix(c(0,0,0,1,0,0,1,1,0),nrow=3,ncol=3)) || !(s2$num.merges == 3)){
  print("problem on merges with p2 restricted with (1,1,2,2,3)")
}
if(!all(s1$nums.divisions == c(1,1,0)) || !(s1$num.divisions == 2)){
  print("problem on divisions with p2 with (1,1,2,2,3)")
}
if(!all(s2$nums.divisions == c(1,1,0)) || !(s2$num.divisions == 2)){
  print("problem on divisions with p2 restricted with (1,1,2,2,3)")
}

s1 <- compute_size_neighborhood_p3(partition)
s2 <- compute_size_neighborhood_p3_restricted(partition,1:n)
if(!all(s1$nums.swaps == c(3,2,3,2,2)) || !(s1$num.swaps == 12)){
  print("problem on swaps with p3 with (1,1,2,2,3)")
}
if(!all(s2$nums.swaps == c(3,2,3,2,2)) || !(s2$num.swaps == 12)){
  print("problem on swaps with p3 restricted with (1,1,2,2,3)")
}

# s1 <- compute_size_neighborhood_p4(partition)
# s2 <- compute_size_neighborhood_p4_restricted(partition,1:n)
# if(!all(s1$nums.swaps == c(2,1,0)) || !(s1$num.swaps == 3)){
#   print("problem on swaps with p4 with (1,1,2,2,3)")
# }
# if(!all(s2$nums.swaps == c(2,1,0)) || !(s2$num.swaps == 3)){
#   print("problem on swaps with p4 restricted with (1,1,2,2,3)")
# }
# 
# s1 <- compute_size_neighborhood_p5(partition)
# s2 <- compute_size_neighborhood_p5_restricted(partition,1:n)
# if(!(s1$num.swaps == 0)){
#   print("problem on swaps with p5 with (1,1,2,2,3)")
# }
# if(!(s2$num.swaps == 0)){
#   print("problem on swaps with p5 restricted with (1,1,2,2,3)")
# }
# 
# s1 <- compute_size_neighborhood_p6(partition)
# s2 <- compute_size_neighborhood_p6_restricted(partition,1:n)
# if(!all(s1$nums.swaps == matrix(c(0,0,0,2,0,0,2,2,0),nrow=3,ncol=3)) || !(s1$num.swaps == 6)){
#   print("problem on swaps with p6 with (1,1,2,2,3)")
# }
# if(!all(s2$nums.swaps == matrix(c(0,0,0,2,0,0,2,2,0),nrow=3,ncol=3)) || !(s2$num.swaps == 6)){
#   print("problem on swaps with p6 restricted with (1,1,2,2,3)")
# }
# 
# partition2 <- 1:n
# s1 <- compute_size_neighborhood_p7(partition, partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition, partition2, 1:n)
# if(!all(s1$nums.swaps == matrix(0,5,5)) || !(s1$num.swaps == 0)){
#   print("problem on swaps with p7 with (1,1,2,2,3) & all isolates")
# }
# if(!all(s2$nums.swaps == matrix(0,5,5)) || !(s2$num.swaps == 0)){
#   print("problem on swaps with p7 restricted with (1,1,2,2,3) & all isolates")
# }
# 
# partition2 <- rep(1,n)
# s1 <- compute_size_neighborhood_p7(partition, partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition, partition2, 1:n)
# if(!all(s1$nums.swaps == matrix(c(0,2,1,1,1,3,0,1,1,1,1,1,0,2,1,1,1,3,0,1,1,1,1,1,0),nrow=5,ncol=5)) || !(s1$num.swaps == 26)){
#   print("problem on swaps with p7 with (1,1,2,2,3) & all in the same group")
# }
# if(!all(s2$nums.swaps == matrix(c(0,2,1,1,1,3,0,1,1,1,1,1,0,2,1,1,1,3,0,1,1,1,1,1,0),nrow=5,ncol=5)) || !(s2$num.swaps == 26)){
#   print("problem on swaps with p7 restricted with (1,1,2,2,3) & all in the same group")
# }

# random case: c(1,1,1,2,2)
partition <- c(1,1,1,2,2)

s1 <- compute_size_neighborhood_p1(partition)
s2 <- compute_size_neighborhood_p1_restricted(partition,1:n)
if(!all(s1$nums.swaps == matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,0),nrow=5,ncol=5)) || !(s1$num.swaps == 6)){
  print("problem on swaps with p1 with (1,1,1,2,2)")
}
if(!all(s2$nums.swaps == matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,0),nrow=5,ncol=5)) || !(s2$num.swaps == 6)){
  print("problem on swaps with p1 restricted with (1,1,1,2,2)")
}

s1 <- compute_size_neighborhood_p2(partition)
s2 <- compute_size_neighborhood_p2_restricted(partition,1:n)
if(!(s1$num.merges == 1)){
  print("problem on merges with p2 with (1,1,1,2,2)")
}
if(!all(s2$merges == matrix(c(0,0,1,0),nrow=2,ncol=2)) || !(s2$num.merges == 1)){
  print("problem on merges with p2 restricted with (1,1,1,2,2)")
}
if(!all(s1$nums.divisions == c(3,1)) || !(s1$num.divisions == 4)){
  print("problem on divisions with p2 with (1,1,1,2,2)")
}
if(!all(s2$nums.divisions == c(3,1)) || !(s2$num.divisions == 4)){
  print("problem on divisions with p2 restricted with (1,1,1,2,2)")
}

s1 <- compute_size_neighborhood_p3(partition)
s2 <- compute_size_neighborhood_p3_restricted(partition,1:n)
if(!all(s1$nums.swaps == c(2,2,2,2,1)) || !(s1$num.swaps == 9)){
  print("problem on swaps with p3 with (1,1,1,2,2)")
}
if(!all(s2$nums.swaps == c(2,2,2,2,1)) || !(s2$num.swaps == 9)){
  print("problem on swaps with p3 restricted with (1,1,1,2,2)")
}

# s1 <- compute_size_neighborhood_p4(partition)
# s2 <- compute_size_neighborhood_p4_restricted(partition,1:n)
# if(!all(s1$nums.swaps == c(6,1)) || !(s1$num.swaps == 7)){
#   print("problem on swaps with p4 with (1,1,1,2,2)")
# }
# if(!all(s2$nums.swaps == c(6,1)) || !(s2$num.swaps == 7)){
#   print("problem on swaps with p4 restricted with (1,1,1,2,2)")
# }
# 
# s1 <- compute_size_neighborhood_p5(partition)
# s2 <- compute_size_neighborhood_p5_restricted(partition,1:n)
# if(!(s1$num.swaps == 3)){
#   print("problem on swaps with p5 with (1,1,1,2,2)")
# }
# if(!(s2$num.swaps == 3)){
#   print("problem on swaps with p5 restricted with (1,1,1,2,2)")
# }
# 
# s1 <- compute_size_neighborhood_p6(partition)
# s2 <- compute_size_neighborhood_p6_restricted(partition,1:n)
# if(!all(s1$nums.swaps == matrix(c(0,0,9,0),nrow=2,ncol=2)) || !(s1$num.swaps == 9)){
#   print("problem on swaps with p6 with (1,1,1,2,2)")
# }
# if(!all(s2$nums.swaps == matrix(c(0,0,9,0),nrow=2,ncol=2)) || !(s2$num.swaps == 9)){
#   print("problem on swaps with p6 restricted with (1,1,1,2,2)")
# }
# 
# partition2 <- 1:n
# s1 <- compute_size_neighborhood_p7(partition, partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition, partition2, 1:n)
# if(!all(s1$nums.swaps == matrix(0,5,5)) || !(s1$num.swaps == 0)){
#   print("problem on swaps with p7 with (1,1,1,2,2) & all isolates")
# }
# if(!all(s2$nums.swaps == matrix(0,5,5)) || !(s2$num.swaps == 0)){
#   print("problem on swaps with p7 restricted with (1,1,1,2,2) & all isolates")
# }
# 
# partition2 <- rep(1,n)
# s1 <- compute_size_neighborhood_p7(partition, partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition, partition2, 1:n)
# if(!all(s1$nums.swaps == matrix(c(0,2,2,1,1,2,0,2,1,1,2,2,0,1,1,1,1,1,0,1,1,1,1,2,0),nrow=5,ncol=5)) || !(s1$num.swaps == 27)){
#   print("problem on swaps with p7 with (1,1,1,2,2) & all in the same group")
# }
# if(!all(s2$nums.swaps == matrix(c(0,2,2,1,1,2,0,2,1,1,2,2,0,1,1,1,1,1,0,1,1,1,1,2,0),nrow=5,ncol=5)) || !(s2$num.swaps == 27)){
#   print("problem on swaps with p7 restricted with (1,1,1,2,2) & all in the same group")
# }


###################################################################################

##### Neighborhoods size tests: check restricted sizes

# random case: c(1,1,2,2,3) and restricted sizes 1 to 3
partition <- c(1,1,2,2,3)

s2 <- compute_size_neighborhood_p1_restricted(partition,1:3)
if(!all(s2$nums.swaps == matrix(c(0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,1,1,1,0),nrow=5,ncol=5)) || !(s2$num.swaps == 8)){
  print("problem on swaps with p1 restricted with (1,1,2,2,3) with restricted sizes 1 to 3")
}

s2 <- compute_size_neighborhood_p2_restricted(partition,1:3)
if(!all(s2$merges == matrix(c(0,0,0,0,0,0,1,1,0),nrow=3,ncol=3)) || !(s2$num.merges == 2)){
  print("problem on merges with p2 restricted with (1,1,2,2,3) with restricted sizes 1 to 3")
}
if(!all(s2$nums.divisions == c(1,1,0)) || !(s2$num.divisions == 2)){
  print("problem on divisions with p2 restricted with (1,1,2,2,3) with restricted sizes 1 to 3")
}

s2 <- compute_size_neighborhood_p3_restricted(partition,1:3)
if(!all(s2$nums.swaps == c(3,2,3,2,2)) || !(s2$num.swaps == 12)){
  print("problem on swaps with p3 restricted with (1,1,2,2,3) with restricted sizes 1 to 3")
}

# s2 <- compute_size_neighborhood_p4_restricted(partition,1:3)
# if(!all(s2$nums.swaps == c(1,1,0)) || !(s2$num.swaps == 2)){
#   print("problem on swaps with p4 restricted with (1,1,2,2,3) with restricted sizes 1 to 3")
# }
# 
# s2 <- compute_size_neighborhood_p5_restricted(partition,1:3)
# if(!(s2$num.swaps == 0)){
#   print("problem on swaps with p5 restricted with (1,1,2,2,3) with restricted sizes 1 to 3")
# }
# 
# s2 <- compute_size_neighborhood_p6_restricted(partition,1:3)
# if(!all(s2$nums.swaps == matrix(c(0,0,0,2,0,0,2,2,0),nrow=3,ncol=3)) || !(s2$num.swaps == 6)){
#   print("problem on swaps with p6 restricted with (1,1,2,2,3) with restricted sizes 1 to 3")
# }
# 
# s2 <- compute_size_neighborhood_p6_restricted(partition,1:3)
# if(!all(s2$nums.swaps == matrix(c(0,0,0,2,0,0,2,2,0),nrow=3,ncol=3)) || !(s2$num.swaps == 6)){
#   print("problem on swaps with p6 restricted with (1,1,2,2,3) with restricted sizes 1 to 3")
# }
# 
# partition2 <- 1:n
# s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,1:3)
# if(!all(s2$nums.swaps == matrix(0,n,n)) || !(s2$num.swaps == 0)){
#   print("problem on swaps with p7 restricted with (1,1,2,2,3) with restricted sizes 1 to 3")
# }
# 
# partition2 <- rep(1,n)
# s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,1:3)
# if(!all(s2$nums.swaps == matrix(c(0,2,1,1,1,3,0,1,1,1,1,1,0,2,1,1,1,3,0,1,1,1,1,1,0),nrow=5,ncol=5)) || !(s2$num.swaps == 26)){
#   print("problem on swaps with p7 restricted with (1,1,2,2,3) with restricted sizes 1 to 3")
# }


# random case: c(1,1,1,2,2) with restricted sizes 2 to 4
partition <- c(1,1,1,2,2)

s2 <- compute_size_neighborhood_p1_restricted(partition,2:4)
if(!all(s2$nums.swaps == matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,0,0),nrow=5,ncol=5)) || !(s2$num.swaps == 6)){
  print("problem on swaps with p1 restricted with (1,1,1,2,2) with restricted sizes 2 to 4")
}

s2 <- compute_size_neighborhood_p2_restricted(partition,2:4)
if(!all(s2$merges == matrix(c(0,0,0,0),nrow=2,ncol=2)) || !(s2$num.merges == 0)){
  print("problem on merges with p2 restricted with (1,1,1,2,2) with restricted sizes 2 to 4")
}
if(!all(s2$nums.divisions == c(0,0)) || !(s2$num.divisions == 0)){
  print("problem on divisions with p2 restricted with (1,1,1,2,2) with restricted sizes 2 to 4")
}

s2 <- compute_size_neighborhood_p3_restricted(partition,2:4)
if(!all(s2$nums.swaps == c(1,1,1,0,0)) || !(s2$num.swaps == 3)){
  print("problem on swaps with p3 restricted with (1,1,1,2,2) with restricted sizes 2 to 4")
}

# s2 <- compute_size_neighborhood_p4_restricted(partition,2:4)
# if(!all(s2$nums.swaps == c(0,0)) || !(s2$num.swaps == 0)){
#   print("problem on swaps with p4 restricted with (1,1,1,2,2) with restricted sizes 2 to 4")
# }
# 
# s2 <- compute_size_neighborhood_p5_restricted(partition,2:4)
# if(!(s2$num.swaps == 3)){
#   print("problem on swaps with p5 restricted with (1,1,1,2,2) with restricted sizes 2 to 4")
# }
# 
# partition2 <- 1:n
# s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,2:4)
# if(!all(s2$nums.swaps == matrix(0,n,n)) || !(s2$num.swaps == 0)){
#   print("problem on swaps with p7 restricted with (1,1,1,2,2) with restricted sizes 2 to 4")
# }
# 
# partition2 <- rep(1,n)
# s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,2:4)
# if(!all(s2$nums.swaps == matrix(c(0,1,1,0,0,1,0,1,0,0,1,1,0,0,0,1,1,1,0,0,1,1,1,0,0),nrow=5,ncol=5)) || !(s2$num.swaps == 12)){
#   print("problem on swaps with p7 restricted with (1,1,1,2,2) with restricted sizes 2 to 4")
# }


###################################################################################


##### Neighborhoods size tests: same results with restricted on all sizes
n <- 5
S <- 100
for(i in 1:S){

  partition <- sample(1:n,n,replace=T)
  partition <- order_groupids(partition)
  
  s1 <- compute_size_neighborhood_p1(partition)
  s2 <- compute_size_neighborhood_p1_restricted(partition,1:n)
  if(!all(s1$nums.swaps == s2$nums.swaps)){
    print("problem on all swaps for p1 with ")
    print(partition)
  }
  if(!(s1$num.swaps == s2$num.swaps)){
    print("problem on total swaps for p1 with ")
    print(partition)
  }
  
  s1 <- compute_size_neighborhood_p2(partition)
  s2 <- compute_size_neighborhood_p2_restricted(partition,1:n)
  if(!(s1$num.merges == s2$num.merges)){
    print("problem on total merges for p2 with ")
    print(partition)
  }
  if(!(s1$num.divisions == s2$num.divisions)){
    print("problem on total divisions for p2 with ")
    print(partition)
  }
  if(!all(s1$nums.divisions == s2$nums.divisions)){
    print("problem on all divisions for p2 with ")
    print(partition)
  }
  
  s1 <- compute_size_neighborhood_p3(partition)
  s2 <- compute_size_neighborhood_p3_restricted(partition,1:n)
  if(!(s1$num.swaps == s2$num.swaps)){
    print("problem on total swaps for p3 with ")
    print(partition)
  } 
  if(!all(s1$nums.swaps == s2$nums.swaps)){
    print("problem on all swaps for p3 with ")
    print(partition)
  } 
  
#   s1 <- compute_size_neighborhood_p4(partition)
#   s2 <- compute_size_neighborhood_p4_restricted(partition,1:n)
#   if(!(s1$num.swaps == s2$num.swaps)){
#     print("problem on total swaps for p4 with ")
#     print(partition)
#   } 
#   if(!all(s1$nums.swaps == s2$nums.swaps)){
#     print("problem on all swaps for p4 with ")
#     print(partition)
#   }
#   
#   s1 <- compute_size_neighborhood_p5(partition)
#   s2 <- compute_size_neighborhood_p5_restricted(partition,1:n)
#   if(!(s1$num.swaps == s2$num.swaps)){
#     print("problem on total swaps for p5 with ")
#     print(partition)
#   } 
#   if(!all(s1$groups.paired.actors == s2$groups.paired.actors)){
#     print("problem on all swaps for p5 with ")
#     print(partition)
#   }
#   
#   s1 <- compute_size_neighborhood_p6(partition)
#   s2 <- compute_size_neighborhood_p6_restricted(partition,1:n)
#   if(!(s1$num.swaps == s2$num.swaps)){
#     print("problem on total swaps for p6 with ")
#     print(partition)
#   } 
#   if(!all(s1$nums.swaps == s2$nums.swaps)){
#     print("problem on all swaps for p6 with ")
#     print(partition)
#   } 
#   
#   partition2 <- 1:n
#   s1 <- compute_size_neighborhood_p7(partition,partition2)
#   s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,1:n)
#   if(!(s1$num.swaps == s2$num.swaps)){
#     print("problem on total swaps for p7 from all isolates with ")
#     print(partition)
#   } 
#   if(!all(s1$nums.swaps == s2$nums.swaps)){
#     print("problem on all swaps for p7 from all isolates with ")
#     print(partition)
#   } 
#   
#   partition2 <- rep(1,n)
#   s1 <- compute_size_neighborhood_p7(partition,partition2)
#   s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,1:n)
#   if(!(s1$num.swaps == s2$num.swaps)){
#     print("problem on total swaps for p7 from all in the same group with ")
#     print(partition)
#   } 
#   if(!all(s1$nums.swaps == s2$nums.swaps)){
#     print("problem on all swaps for p7 from all in the same group with ")
#     print(partition)
#   } 
}


###################################################################################


##### Neighborhoods sample tests: check some samplings of certain partitions
S <- 500

# extreme case: all isolates
partition <- 1:n

s1 <- compute_size_neighborhood_p1(partition)
s2 <- compute_size_neighborhood_p1_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p1(partition,s1)
  allsamples2[i,] <- sample_new_partition_p1_restricted(partition,s2,1:n)
}
found1 <- F
found2 <- F
for(i in 1:S){
 if(all(1:n == allsamples1[i,])) { found1 = T}
 if(all(1:n == allsamples2[i,])) { found2 = T}
}
if(!found1){
  print("problem on sampling with p1 with all isolates")
}
if(!found2){
  print("problem on sampling with p1 restricted with all isolates")
}

s1 <- compute_size_neighborhood_p2(partition)
s2 <- compute_size_neighborhood_p2_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p2(partition,s1)
  allsamples2[i,] <- sample_new_partition_p2_restricted(partition,s2,1:n,)
}
found1 <- F
found2 <- F
for(i in 1:S){
  if(all(c(1,2,1,3,4) == allsamples1[i,])) { found1 = T}
  if(all(c(1,2,1,3,4) == allsamples2[i,])) { found2 = T}
}
if(!found1){
  print("problem on sampling with p2 with all isolates")
}
if(!found2){
  print("problem on sampling with p2 restricted with all isolates")
}

s1 <- compute_size_neighborhood_p3(partition)
s2 <- compute_size_neighborhood_p3_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p3(partition,s1)
  allsamples2[i,] <- sample_new_partition_p3_restricted(partition,s2,1:n)
}
found1 <- F
found2 <- F
for(i in 1:S){
  if(all(c(1,2,3,4,1) == allsamples1[i,])) { found1 = T}
  if(all(c(1,2,3,4,1) == allsamples2[i,])) { found2 = T}
}
if(!found1){
  print("problem on sampling with p3 with all isolates")
}
if(!found2){
  print("problem on sampling with p3 restricted with all isolates")
}

# s1 <- compute_size_neighborhood_p4(partition)
# s2 <- compute_size_neighborhood_p4_restricted(partition,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p4(partition,s1)
#   allsamples2[i,] <- sample_new_partition_p4_restricted(partition,1:n,s2)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(1:n == allsamples1[i,])) { found1 = T}
#   if(all(1:n == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p4 with all isolates")
# }
# if(!found2){
#   print("problem on sampling with p4 restricted with all isolates")
# }
# 
# s1 <- compute_size_neighborhood_p5(partition)
# s2 <- compute_size_neighborhood_p5_restricted(partition,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p5(partition,s1)
#   allsamples2[i,] <- sample_new_partition_p5_restricted(partition,1:n,s2)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(1:n == allsamples1[i,])) { found1 = T}
#   if(all(1:n == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p5 with all isolates")
# }
# if(!found2){
#   print("problem on sampling with p5 restricted with all isolates")
# }
# 
# s1 <- compute_size_neighborhood_p6(partition)
# s2 <- compute_size_neighborhood_p6_restricted(partition,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p6(partition,s1)
#   allsamples2[i,] <- sample_new_partition_p6_restricted(partition,1:n,s2)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(1:n == allsamples1[i,])) { found1 = T}
#   if(all(1:n == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p6 with all isolates")
# }
# if(!found2){
#   print("problem on sampling with p6 restricted with all isolates")
# }
# 
# partition2 <- 1:n
# s1 <- compute_size_neighborhood_p7(partition,partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p7(partition,partition2,s1)
#   allsamples2[i,] <- sample_new_partition_p7_restricted(partition,partition2,s2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(1:n == allsamples1[i,])) { found1 = T}
#   if(all(1:n == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p7 with all isolates & all isolates")
# }
# if(!found2){
#   print("problem on sampling with p7 restricted with all isolates & all isolates")
# }
# 
# partition2 <- rep(1,n)
# s1 <- compute_size_neighborhood_p7(partition,partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p7(partition,partition2,s1)
#   allsamples2[i,] <- sample_new_partition_p7_restricted(partition,partition2,s2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(c(1,2,3,4,1) == allsamples1[i,])) { found1 = T}
#   if(all(c(1,2,3,4,1) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p7 with all isolates & all in the same group")
# }
# if(!found2){
#   print("problem on sampling with p7 restricted with all isolates & all in the same group")
# }


# extreme case: all in the same group
partition <- rep(1,n)

s1 <- compute_size_neighborhood_p1(partition)
s2 <- compute_size_neighborhood_p1_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p1(partition,s1)
  allsamples2[i,] <- sample_new_partition_p1_restricted(partition,s2,1:n)
}
found1 <- F
found2 <- F
for(i in 1:S){
  if(all(rep(1,n) == allsamples1[i,])) { found1 = T}
  if(all(rep(1,n) == allsamples2[i,])) { found2 = T}
}
if(!found1){
  print("problem on sampling with p1 with all in the same group")
}
if(!found2){
  print("problem on sampling with p1 restricted with all in the same group")
}

s1 <- compute_size_neighborhood_p2(partition)
s2 <- compute_size_neighborhood_p2_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p2(partition,s1)
  allsamples2[i,] <- sample_new_partition_p2_restricted(partition,s2,1:n)
}
found1 <- F
found2 <- F
for(i in 1:S){
  if(all(c(1,1,2,2,2) == allsamples1[i,])) { found1 = T}
  if(all(c(1,1,2,2,2) == allsamples2[i,])) { found2 = T}
}
if(!found1){
  print("problem on sampling with p2 with all in the same group")
}
if(!found2){
  print("problem on sampling with p2 restricted with all in the same group")
}

s1 <- compute_size_neighborhood_p3(partition)
s2 <- compute_size_neighborhood_p3_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p3(partition,s1)
  allsamples2[i,] <- sample_new_partition_p3_restricted(partition,s2,1:n)
}
found1 <- F
found2 <- F
for(i in 1:S){
  if(all(c(1,2,2,2,2) == allsamples1[i,])) { found1 = T}
  if(all(c(1,2,2,2,2) == allsamples2[i,])) { found2 = T}
}
if(!found1){
  print("problem on sampling with p3 with all in the same group")
}
if(!found2){
  print("problem on sampling with p3 restricted with all in the same group")
}

# s1 <- compute_size_neighborhood_p4(partition)
# s2 <- compute_size_neighborhood_p4_restricted(partition,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p4(partition,s1)
#   allsamples2[i,] <- sample_new_partition_p4_restricted(partition,s2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(c(1,2,1,2,2) == allsamples1[i,])) { found1 = T}
#   if(all(c(1,2,1,2,2) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p4 with all in the same group")
# }
# if(!found2){
#   print("problem on sampling with p4 restricted with all in the same group")
# }
# 
# s1 <- compute_size_neighborhood_p5(partition)
# s2 <- compute_size_neighborhood_p5_restricted(partition,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p5(partition,s1)
#   allsamples2[i,] <- sample_new_partition_p5_restricted(partition,s2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(rep(1,n) == allsamples1[i,])) { found1 = T}
#   if(all(rep(1,n) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p5 with all in the same group")
# }
# if(!found2){
#   print("problem on sampling with p5 restricted with all in the same group")
# }
# 
# s1 <- compute_size_neighborhood_p6(partition)
# s2 <- compute_size_neighborhood_p6_restricted(partition,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p6(partition,s1)
#   allsamples2[i,] <- sample_new_partition_p6_restricted(partition,s2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(rep(1,n) == allsamples1[i,])) { found1 = T}
#   if(all(rep(1,n) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p6 with all in the same group")
# }
# if(!found2){
#   print("problem on sampling with p6 restricted with all in the same group")
# }
# 
# partition2 <- 1:n
# s1 <- compute_size_neighborhood_p7(partition,partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p7(partition,partition2,s1)
#   allsamples2[i,] <- sample_new_partition_p7_restricted(partition,s2,partition2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(rep(1,n) == allsamples1[i,])) { found1 = T}
#   if(all(rep(1,n) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p7 with all in the same group & all isolates")
# }
# if(!found2){
#   print("problem on sampling with p7 restricted with all in the same group & all isolates")
# }
# 
# partition2 <- rep(1,n)
# s1 <- compute_size_neighborhood_p7(partition,partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p7(partition,partition2,s1)
#   allsamples2[i,] <- sample_new_partition_p7_restricted(partition,s2,partition2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(c(1,1,2,1,1) == allsamples1[i,])) { found1 = T}
#   if(all(c(1,1,2,1,1) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p7 with all in the same group & all in the same group")
# }
# if(!found2){
#   print("problem on sampling with p7 restricted with all in the same group & all in the same group")
# }


# random case: c(1,1,2,2,3)
partition <- c(1,1,2,2,3)

s1 <- compute_size_neighborhood_p1(partition)
s2 <- compute_size_neighborhood_p1_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p1(partition,s1)
  allsamples2[i,] <- sample_new_partition_p1_restricted(partition,s2,1:n)
}
found1 <- F
found2 <- F
for(i in 1:S){
  if(all(c(1,2,1,2,3) == allsamples1[i,])) { found1 = T}
  if(all(c(1,2,1,2,3) == allsamples2[i,])) { found2 = T}
}
if(!found1){
  print("problem on sampling with p1 with (1,1,2,2,3)")
}
if(!found2){
  print("problem on sampling with p1 restricted with (1,1,2,2,3)")
}

s1 <- compute_size_neighborhood_p2(partition)
s2 <- compute_size_neighborhood_p2_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p2(partition,s1)
  allsamples2[i,] <- sample_new_partition_p2_restricted(partition,s2,1:n)
}
found1 <- F
found2 <- F
for(i in 1:S){
  if(all(c(1,1,1,1,2) == allsamples1[i,])) { found1 = T}
  if(all(c(1,1,1,1,2) == allsamples2[i,])) { found2 = T}
}
if(!found1){
  print("problem on sampling with p2 with (1,1,2,2,3)")
}
if(!found2){
  print("problem on sampling with p2 restricted with (1,1,2,2,3)")
}

s1 <- compute_size_neighborhood_p3(partition)
s2 <- compute_size_neighborhood_p3_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p3(partition,s1)
  allsamples2[i,] <- sample_new_partition_p3_restricted(partition,s2,1:n)
}
found1 <- F
found2 <- F
for(i in 1:S){
  if(all(c(1,1,2,3,3) == allsamples1[i,])) { found1 = T}
  if(all(c(1,1,2,3,3) == allsamples2[i,])) { found2 = T}
}
if(!found1){
  print("problem on sampling with p3 with (1,1,2,2,3)")
}
if(!found2){
  print("problem on sampling with p3 restricted with (1,1,2,2,3)")
}

# s1 <- compute_size_neighborhood_p4(partition)
# s2 <- compute_size_neighborhood_p4_restricted(partition,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p4(partition,s1)
#   allsamples2[i,] <- sample_new_partition_p4_restricted(partition,s2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(c(1,1,2,2,2) == allsamples1[i,])) { found1 = T}
#   if(all(c(1,1,2,2,2) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p4 with (1,1,2,2,3)")
# }
# if(!found2){
#   print("problem on sampling with p4 restricted with (1,1,2,2,3)")
# }
# 
# s1 <- compute_size_neighborhood_p5(partition)
# s2 <- compute_size_neighborhood_p5_restricted(partition,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p5(partition,s1)
#   allsamples2[i,] <- sample_new_partition_p5_restricted(partition,s2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(c(1,1,2,2,3) == allsamples1[i,])) { found1 = T}
#   if(all(c(1,1,2,2,3) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p5 with (1,1,2,2,3)")
# }
# if(!found2){
#   print("problem on sampling with p5 restricted with (1,1,2,2,3)")
# }
# 
# s1 <- compute_size_neighborhood_p6(partition)
# s2 <- compute_size_neighborhood_p6_restricted(partition,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p6(partition,s1)
#   allsamples2[i,] <- sample_new_partition_p6_restricted(partition,s2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(c(1,2,3,3,1) == allsamples1[i,])) { found1 = T}
#   if(all(c(1,2,3,3,1) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p6 with (1,1,2,2,3)")
# }
# if(!found2){
#   print("problem on sampling with p6 restricted with (1,1,2,2,3)")
# }
# 
# partition2 <- 1:n
# s1 <- compute_size_neighborhood_p7(partition,partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p7(partition,partition2,s1)
#   allsamples2[i,] <- sample_new_partition_p7_restricted(partition,s2,partition2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(c(1,1,2,2,3) == allsamples1[i,])) { found1 = T}
#   if(all(c(1,1,2,2,3) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p7 with (1,1,2,2,3) & all isolates")
# }
# if(!found2){
#   print("problem on sampling with p7 restricted with (1,1,2,2,3) & all isolates")
# }
# 
# partition2 <- rep(1,n)
# s1 <- compute_size_neighborhood_p7(partition,partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p7(partition,partition2,s1)
#   allsamples2[i,] <- sample_new_partition_p7_restricted(partition,s2,partition2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(c(1,1,2,3,3) == allsamples1[i,])) { found1 = T}
#   if(all(c(1,1,2,3,3) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p7 with (1,1,2,2,3) & all in the same group")
# }
# if(!found2){
#   print("problem on sampling with p7 restricted with (1,1,2,2,3) & all in the same group")
# }


# random case: (1,1,1,2,2)
partition <- c(1,1,1,2,2)

s1 <- compute_size_neighborhood_p1(partition)
s2 <- compute_size_neighborhood_p1_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p1(partition,s1)
  allsamples2[i,] <- sample_new_partition_p1_restricted(partition,s2,1:n)
}
found1 <- F
found2 <- F
for(i in 1:S){
  if(all(c(1,2,1,1,2) == allsamples1[i,])) { found1 = T}
  if(all(c(1,2,1,1,2) == allsamples2[i,])) { found2 = T}
}
if(!found1){
  print("problem on sampling with p1 with (1,1,1,2,2)")
}
if(!found2){
  print("problem on sampling with p1 restricted with (1,1,1,2,2)")
}

s1 <- compute_size_neighborhood_p2(partition)
s2 <- compute_size_neighborhood_p2_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p2(partition,s1)
  allsamples2[i,] <- sample_new_partition_p2_restricted(partition,s2,1:n)
}
found1 <- F
found2 <- F
for(i in 1:S){
  if(all(c(1,2,2,3,3) == allsamples1[i,])) { found1 = T}
  if(all(c(1,2,2,3,3) == allsamples2[i,])) { found2 = T}
}
if(!found1){
  print("problem on sampling with p2 with (1,1,1,2,2)")
}
if(!found2){
  print("problem on sampling with p2 restricted with (1,1,1,2,2)")
}

s1 <- compute_size_neighborhood_p3(partition)
s2 <- compute_size_neighborhood_p3_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p3(partition,s1)
  allsamples2[i,] <- sample_new_partition_p3_restricted(partition,s2,1:n)
}
found1 <- F
found2 <- F
for(i in 1:S){
  if(all(c(1,1,2,2,2) == allsamples1[i,])) { found1 = T}
  if(all(c(1,1,2,2,2) == allsamples2[i,])) { found2 = T}
}
if(!found1){
  print("problem on sampling with p3 with (1,1,1,2,2)")
}
if(!found2){
  print("problem on sampling with p3 restricted with (1,1,1,2,2)")
}

# s1 <- compute_size_neighborhood_p4(partition)
# s2 <- compute_size_neighborhood_p4_restricted(partition,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p4(partition,s1)
#   allsamples2[i,] <- sample_new_partition_p4_restricted(partition,s2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(c(1,1,2,1,1) == allsamples1[i,])) { found1 = T}
#   if(all(c(1,1,2,1,1) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p4 with (1,1,1,2,2)")
# }
# if(!found2){
#   print("problem on sampling with p4 restricted with (1,1,1,2,2)")
# }
# 
# s1 <- compute_size_neighborhood_p5(partition)
# s2 <- compute_size_neighborhood_p5_restricted(partition,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p5(partition,s1)
#   allsamples2[i,] <- sample_new_partition_p5_restricted(partition,s2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(c(1,1,2,2,2) == allsamples1[i,])) { found1 = T}
#   if(all(c(1,1,2,2,2) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p5 with (1,1,1,2,2)")
# }
# if(!found2){
#   print("problem on sampling with p5 restricted with (1,1,1,2,2)")
# }
# 
# s1 <- compute_size_neighborhood_p6(partition)
# s2 <- compute_size_neighborhood_p6_restricted(partition,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p6(partition,s1)
#   allsamples2[i,] <- sample_new_partition_p6_restricted(partition,s2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(c(1,2,1,2,1) == allsamples1[i,])) { found1 = T}
#   if(all(c(1,2,1,2,1) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p6 with (1,1,1,2,2)")
# }
# if(!found2){
#   print("problem on sampling with p6 restricted with (1,1,1,2,2)")
# }
# 
# partition2 <- 1:n
# s1 <- compute_size_neighborhood_p7(partition,partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p7(partition,partition2,s1)
#   allsamples2[i,] <- sample_new_partition_p7_restricted(partition,s2,partition2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(c(1,1,1,2,2) == allsamples1[i,])) { found1 = T}
#   if(all(c(1,1,1,2,2) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p7 with (1,1,1,2,2) & all isolates")
# }
# if(!found2){
#   print("problem on sampling with p7 restricted with (1,1,1,2,2) & all isolates")
# }
# 
# partition2 <- rep(1,n)
# s1 <- compute_size_neighborhood_p7(partition,partition2)
# s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,1:n)
# allsamples1 <- matrix(0,S,n)
# allsamples2 <- matrix(0,S,n)
# for(i in 1:S) {
#   allsamples1[i,] <- sample_new_partition_p7(partition,partition2,s1)
#   allsamples2[i,] <- sample_new_partition_p7_restricted(partition,s2,partition2,1:n)
# }
# found1 <- F
# found2 <- F
# for(i in 1:S){
#   if(all(c(1,1,2,3,3) == allsamples1[i,])) { found1 = T}
#   if(all(c(1,1,2,3,3) == allsamples2[i,])) { found2 = T}
# }
# if(!found1){
#   print("problem on sampling with p7 with (1,1,1,2,2) & all in the same group")
# }
# if(!found2){
#   print("problem on sampling with p7 restricted with (1,1,1,2,2) & all in the same group")
# }


###################################################################################


##### Neighborhoods sample tests: check that both restricted and non restricted 
# versions sample the same partitions

for(i in 1:50){
  
  if(i == 1) partition <- 1:n
  if(i == 2) partition <- rep(1,n)
  if(i >= 3) partition <- order_groupids(sample(n,n,replace=T))
  S <- 500

  s1 <- compute_size_neighborhood_p1(partition)
  s2 <- compute_size_neighborhood_p1_restricted(partition,1:n)
  allsamples1 <- matrix(0,S,n)
  allsamples2 <- matrix(0,S,n)
  for(i in 1:S) {
    allsamples1[i,] <- sample_new_partition_p1(partition,s1)
    allsamples2[i,] <- sample_new_partition_p1_restricted(partition,s2,1:n)
  } 
  notin1 <- c()
  notin2 <- c()
  for(i in 1:S){
    p1 <- allsamples1[i,]
    p2 <- allsamples2[i,]
    found1 <- F
    found2 <- F
    for(j in 1:S) {
      if(all(p1 == allsamples2[j,])) { found1 = T}
      if(all(p2 == allsamples1[j,])) { found2 = T}
    }
    if(!found1) {notin1 <- rbind(notin1,p1)}
    if(!found2) {notin2 <- rbind(notin2,p2)}
  }

  if(!is.null(notin1)) {
    print(partition)
    print("partitions sampled normally but not sampled in the restricted version with p1")
    print(notin1)
  }
  if(!is.null(notin2)) {
    print(partition)
    print("partitions sampled in the restricted version but not normally with p1")
    print(notin2)
  }
  
  s1 <- compute_size_neighborhood_p2(partition)
  s2 <- compute_size_neighborhood_p2_restricted(partition,1:n)
  allsamples1 <- matrix(0,S,n)
  allsamples2 <- matrix(0,S,n)
  for(i in 1:S) {
    allsamples1[i,] <- sample_new_partition_p2(partition,s1)
    allsamples2[i,] <- sample_new_partition_p2_restricted(partition,s2,1:n)
  } 
  notin1 <- c()
  notin2 <- c()
  for(i in 1:S){
    p1 <- allsamples1[i,]
    p2 <- allsamples2[i,]
    found1 <- F
    found2 <- F
    for(j in 1:S) {
      if(all(p1 == allsamples2[j,])) { found1 = T}
      if(all(p2 == allsamples1[j,])) { found2 = T}
    }
    if(!found1) {notin1 <- rbind(notin1,p1)}
    if(!found2) {notin2 <- rbind(notin2,p2)}
  }
  if(!is.null(notin1)) {
    print(partition)
    print("partitions sampled normally but not sampled in the restricted version with p2")
    print(notin1)
  }
  if(!is.null(notin2)) {
    print(partition)
    print("partitions sampled in the restricted version but not normally with p2")
    print(notin2)
  }

  s1 <- compute_size_neighborhood_p3(partition)
  s2 <- compute_size_neighborhood_p3_restricted(partition,1:n)
  allsamples1 <- matrix(0,S,n)
  allsamples2 <- matrix(0,S,n)
  for(i in 1:S) {
    allsamples1[i,] <- sample_new_partition_p3(partition,s1)
    allsamples2[i,] <- sample_new_partition_p3_restricted(partition,s2,1:n)
  } 
  notin1 <- c()
  notin2 <- c()
  for(i in 1:S){
    p1 <- allsamples1[i,]
    p2 <- allsamples2[i,]
    found1 <- F
    found2 <- F
    for(j in 1:S) {
      if(all(p1 == allsamples2[j,])) { found1 = T}
      if(all(p2 == allsamples1[j,])) { found2 = T}
    }
    if(!found1) {notin1 <- rbind(notin1,p1)}
    if(!found2) {notin2 <- rbind(notin2,p2)}
  }
  if(!is.null(notin1)) {
    print(partition)
    print("partitions sampled normally but not sampled in the restricted version with p3")
    print(notin1)
  }
  if(!is.null(notin2)) {
    print(partition)
    print("partitions sampled in the restricted version but not normally with p3")
    print(notin2)
  }
  
#   s1 <- compute_size_neighborhood_p4(partition)
#   s2 <- compute_size_neighborhood_p4_restricted(partition,1:n)
#   allsamples1 <- matrix(0,S,n)
#   allsamples2 <- matrix(0,S,n)
#   for(i in 1:S) {
#     allsamples1[i,] <- sample_new_partition_p4(partition,s1)
#     allsamples2[i,] <- sample_new_partition_p4_restricted(partition,s2,1:n)
#   } 
#   notin1 <- c()
#   notin2 <- c()
#   for(i in 1:S){
#     p1 <- allsamples1[i,]
#     p2 <- allsamples2[i,]
#     found1 <- F
#     found2 <- F
#     for(j in 1:S) {
#       if(all(p1 == allsamples2[j,])) { found1 = T}
#       if(all(p2 == allsamples1[j,])) { found2 = T}
#     }
#     if(!found1) {notin1 <- rbind(notin1,p1)}
#     if(!found2) {notin2 <- rbind(notin2,p2)}
#   }
#   if(!is.null(notin1)) {
#     print(partition)
#     print("partitions sampled normally but not sampled in the restricted version with p4")
#     print(notin1)
#   }
#   if(!is.null(notin2)) {
#     print(partition)
#     print("partitions sampled in the restricted version but not normally with p4")
#     print(notin2)
#   }
# 
#   s1 <- compute_size_neighborhood_p5(partition)
#   s2 <- compute_size_neighborhood_p5_restricted(partition,1:n)
#   allsamples1 <- matrix(0,S,n)
#   allsamples2 <- matrix(0,S,n)
#   for(i in 1:S) {
#     allsamples1[i,] <- sample_new_partition_p5(partition,s1)
#     allsamples2[i,] <- sample_new_partition_p5_restricted(partition,s2,1:n)
#   } 
#   notin1 <- c()
#   notin2 <- c()
#   for(i in 1:S){
#     p1 <- allsamples1[i,]
#     p2 <- allsamples2[i,]
#     found1 <- F
#     found2 <- F
#     for(j in 1:S) {
#       if(all(p1 == allsamples2[j,])) { found1 = T}
#       if(all(p2 == allsamples1[j,])) { found2 = T}
#     }
#     if(!found1) {notin1 <- rbind(notin1,p1)}
#     if(!found2) {notin2 <- rbind(notin2,p2)}
#   }
#   if(!is.null(notin1)) {
#     print(partition)
#     print("partitions sampled normally but not sampled in the restricted version with p5")
#     print(notin1)
#   }
#   if(!is.null(notin2)) {
#     print(partition)
#     print("partitions sampled in the restricted version but not normally with p5")
#     print(notin2)
#   }
#   
#   s1 <- compute_size_neighborhood_p6(partition)
#   s2 <- compute_size_neighborhood_p6_restricted(partition,1:n)
#   allsamples1 <- matrix(0,S,n)
#   allsamples2 <- matrix(0,S,n)
#   for(i in 1:S) {
#     allsamples1[i,] <- sample_new_partition_p6(partition,s1)
#     allsamples2[i,] <- sample_new_partition_p6_restricted(partition,s2,1:n)
#   } 
#   notin1 <- c()
#   notin2 <- c()
#   for(i in 1:S){
#     p1 <- allsamples1[i,]
#     p2 <- allsamples2[i,]
#     found1 <- F
#     found2 <- F
#     for(j in 1:S) {
#       if(all(p1 == allsamples2[j,])) { found1 = T}
#       if(all(p2 == allsamples1[j,])) { found2 = T}
#     }
#     if(!found1) {notin1 <- rbind(notin1,p1)}
#     if(!found2) {notin2 <- rbind(notin2,p2)}
#   }
#   if(!is.null(notin1)) {
#     print(partition)
#     print("partitions sampled normally but not sampled in the restricted version with p6")
#     print(notin1)
#   }
#   if(!is.null(notin2)) {
#     print(partition)
#     print("partitions sampled in the restricted version but not normally with p6")
#     print(notin2)
#   }
#   
#   partition2 <- 1:n
#   s1 <- compute_size_neighborhood_p7(partition,partition2)
#   s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,1:n)
#   allsamples1 <- matrix(0,S,n)
#   allsamples2 <- matrix(0,S,n)
#   for(i in 1:S) {
#     allsamples1[i,] <- sample_new_partition_p7(partition,partition2,s1)
#     allsamples2[i,] <- sample_new_partition_p7_restricted(partition,s2,partition2,1:n)
#   } 
#   notin1 <- c()
#   notin2 <- c()
#   for(i in 1:S){
#     p1 <- allsamples1[i,]
#     p2 <- allsamples2[i,]
#     found1 <- F
#     found2 <- F
#     for(j in 1:S) {
#       if(all(p1 == allsamples2[j,])) { found1 = T}
#       if(all(p2 == allsamples1[j,])) { found2 = T}
#     }
#     if(!found1) {notin1 <- rbind(notin1,p1)}
#     if(!found2) {notin2 <- rbind(notin2,p2)}
#   }
#   if(!is.null(notin1)) {
#     print(partition)
#     print("partitions sampled normally but not sampled in the restricted version with p7 from all isolates")
#     print(notin1)
#   }
#   if(!is.null(notin2)) {
#     print(partition)
#     print("partitions sampled in the restricted version but not normally with p7 from all isolates")
#     print(notin2)
#   }
#   
#   partition2 <- rep(1,n)
#   s1 <- compute_size_neighborhood_p7(partition,partition2)
#   s2 <- compute_size_neighborhood_p7_restricted(partition,partition2,1:n)
#   allsamples1 <- matrix(0,S,n)
#   allsamples2 <- matrix(0,S,n)
#   for(i in 1:S) {
#     allsamples1[i,] <- sample_new_partition_p7(partition,partition2,s1)
#     allsamples2[i,] <- sample_new_partition_p7_restricted(partition,s2,partition2,1:n)
#   } 
#   notin1 <- c()
#   notin2 <- c()
#   for(i in 1:S){
#     p1 <- allsamples1[i,]
#     p2 <- allsamples2[i,]
#     found1 <- F
#     found2 <- F
#     for(j in 1:S) {
#       if(all(p1 == allsamples2[j,])) { found1 = T}
#       if(all(p2 == allsamples1[j,])) { found2 = T}
#     }
#     if(!found1) {notin1 <- rbind(notin1,p1)}
#     if(!found2) {notin2 <- rbind(notin2,p2)}
#   }
#   if(!is.null(notin1)) {
#     print(partition)
#     print("partitions sampled normally but not sampled in the restricted version with p7 from all in the same group")
#     print(notin1)
#   }
#   if(!is.null(notin2)) {
#     print(partition)
#     print("partitions sampled in the restricted version but not normally with p7 from all in the same group")
#     print(notin2)
#   }

}


###################################################################################

##### Neighborhoods reachable tests: check that reachability is correctly tested
# for some extreme cases

partition1 <- c(1,2,3,4,5)
partition2 <- c(1,2,3,4,4)
if(reachable_p1(partition1,partition2)) print("problem with reachable test with p1")
partition1 <- c(1,1,1,1,1)
partition2 <- c(1,2,1,1,1)
if(reachable_p1(partition1,partition2)) print("problem with reachable test with p1")

partition1 <- c(1,2,3,4,5)
partition2 <- c(1,2,3,4,4)
if(!reachable_p2(partition1,partition2)) print("problem with reachable test with p2")
partition1 <- c(1,1,1,1,1)
partition2 <- c(1,2,1,1,1)
if(!reachable_p2(partition1,partition2)) print("problem with reachable test with p2")

partition1 <- c(1,2,3,4,5)
partition2 <- c(1,2,3,4,4)
if(!reachable_p3(partition1,partition2)) print("problem with reachable test with p3")
partition1 <- c(1,1,1,1,1)
partition2 <- c(1,2,1,1,1)
if(!reachable_p3(partition1,partition2)) print("problem with reachable test with p3")


###################################################################################

##### Neighborhoods reachable tests: check that reachability is correctly tested
# for randomly sampled cases

for(i in 1:50){
  
  if(i == 1) partition <- 1:n
  if(i == 2) partition <- rep(1,n)
  if(i >= 3) partition <- order_groupids(sample(n,n,replace=T))
  S <- 500
  
  s1 <- compute_size_neighborhood_p1(partition)
  s2 <- compute_size_neighborhood_p1_restricted(partition,1:n)
  allsamples1 <- matrix(0,S,n)
  allsamples2 <- matrix(0,S,n)
  for(i in 1:S) {
    sample1 <- sample_new_partition_p1(partition,s1)
    sample2 <- sample_new_partition_p1_restricted(partition,s2,1:n)
    if(s1$total != 0 && !reachable_p1(partition,sample1)) {
      print("problem with reachable test with p1 with ")
      print(partition)
      print("and")
      print(sample1)
    }
    if(s2$total != 0 && !reachable_p1(partition,sample2)) {
      print("problem with reachable test with p1 with ")
      print(partition)
      print("and")
      print(sample2)
    } 
  }
  
  s1 <- compute_size_neighborhood_p2(partition)
  s2 <- compute_size_neighborhood_p2_restricted(partition,1:n)
  allsamples1 <- matrix(0,S,n)
  allsamples2 <- matrix(0,S,n)
  for(i in 1:S) {
    sample1 <- sample_new_partition_p2(partition,s1)
    sample2 <- sample_new_partition_p2_restricted(partition,s2,1:n)
    if(!reachable_p2(partition,sample1)) {
      print("problem with reachable test with p2 with ")
      print(partition)
      print("and")
      print(sample1)
    }
    if(!reachable_p2(partition,sample2)) {
      print("problem with reachable test with p2 with ")
      print(partition)
      print("and")
      print(sample2)
    } 
  }
  
  s1 <- compute_size_neighborhood_p3(partition)
  s2 <- compute_size_neighborhood_p3_restricted(partition,1:n)
  allsamples1 <- matrix(0,S,n)
  allsamples2 <- matrix(0,S,n)
  for(i in 1:S) {
    sample1 <- sample_new_partition_p3(partition,s1)
    sample2 <- sample_new_partition_p3_restricted(partition,s2,1:n)
    if(!reachable_p3(partition,sample1)) {
      print("problem with reachable test with p3 with ")
      print(partition)
      print("and")
      print(sample1)
    }
    if(!reachable_p3(partition,sample2)) {
      print("problem with reachable test with p3 with ")
      print(partition)
      print("and")
      print(sample2)
    } 
  }
}

###################################################################################


##### Neighborhoods tests: check that sampling for 10 nodes, with a null model,
# we get an expected number of groups around 4.85 (can be calculated exactly)
# uncomment, and vary the neighborhoods (it's not going to be exactly 4.85)

#draws.p <- draw_Metropolis_single(theta = 0, 
#                                  first.partition = 1:10, 
#                                  nodes = data.frame(id = 1:10), 
#                                  effects = list(names = "num_groups", objects = "partition"),
#                                  objects = list(), 
#                                  burnin = 100, 
#                                  thining = 100, 
#                                  num.steps = 2000,  
#                                  neighborhood = c(0,0,0.2,0,0.8,0), 
#                                  return.all.partitions = T)
#allgroupcounts <- rep(0,10)
#for(d in 1:2000){
#  partition <- draws.p$all.partitions[d,]
#  for(group in 1:max(partition)){
#    size <- sum(partition == group)
#    allgroupcounts[size] <- allgroupcounts[size] + 1
#  }
#}
#sum(allgroupcounts / num.steps) # around 4.85
