# Neighborhoods tests: sizes
n <- 5
S <- 100
for(i in 1:S){
  
  #extreme case: all isolates
  if(i==1) partition <- 1:n
  if(i==2) partition <- rep(1,n)
  if(i>=3) partition <- sample(1:n,n,replace=T)
  
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
}


# Neighborhoods tests: sampling
partition <- 1:n
partition <- rep(1,n)
partition <- c(1,2,2,3,3)
n <- 5
S <- 500

s1 <- compute_size_neighborhood_p1(partition)
s2 <- compute_size_neighborhood_p1_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p1(partition,"normalized",s1)
  allsamples2[i,] <- sample_new_partition_p1_restricted(partition,"normalized",1:n,s2)
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
print("partitions sampled normally but not sampled in the restricted version with p1")
print(notin1)
print("partitions sampled in the restricted version but not normally with p1")
print(notin2)

s1 <- compute_size_neighborhood_p2(partition)
s2 <- compute_size_neighborhood_p2_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p2(partition,"normalized",s1)
  allsamples2[i,] <- sample_new_partition_p2_restricted(partition,"normalized",1:n,s2)
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
print("partitions sampled normally but not sampled in the restricted version with p2")
print(notin1)
print("partitions sampled in the restricted version but not normally with p2")
print(notin2)


s1 <- compute_size_neighborhood_p3(partition)
s2 <- compute_size_neighborhood_p3_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p3(partition,"normalized",s1)
  allsamples2[i,] <- sample_new_partition_p3_restricted(partition,"normalized",1:n,s2)
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
print("partitions sampled normally but not sampled in the restricted version with p3")
print(notin1)
print("partitions sampled in the restricted version but not normally with p3")
print(notin2)



# Neighborhoods tests: sampling
partition <- 1:n
partition <- rep(1,n)
partition <- c(1,2,2,3,3)
n <- 5
S <- 500

s1 <- compute_size_neighborhood_p1(partition)
s2 <- compute_size_neighborhood_p1_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p1(partition,"normalized",s1)
  allsamples2[i,] <- sample_new_partition_p1_restricted(partition,"normalized",1:n,s2)
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
print("partitions sampled normally but not sampled in the restricted version with p1")
print(notin1)
print("partitions sampled in the restricted version but not normally with p1")
print(notin2)

s1 <- compute_size_neighborhood_p2(partition)
s2 <- compute_size_neighborhood_p2_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p2(partition,"normalized",s1)
  allsamples2[i,] <- sample_new_partition_p2_restricted(partition,"normalized",1:n,s2)
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
print("partitions sampled normally but not sampled in the restricted version with p2")
print(notin1)
print("partitions sampled in the restricted version but not normally with p2")
print(notin2)


s1 <- compute_size_neighborhood_p3(partition)
s2 <- compute_size_neighborhood_p3_restricted(partition,1:n)
allsamples1 <- matrix(0,S,n)
allsamples2 <- matrix(0,S,n)
for(i in 1:S) {
  allsamples1[i,] <- sample_new_partition_p3(partition,"normalized",s1)
  allsamples2[i,] <- sample_new_partition_p3_restricted(partition,"normalized",1:n,s2)
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
print("partitions sampled normally but not sampled in the restricted version with p3")
print(notin1)
print("partitions sampled in the restricted version but not normally with p3")
print(notin2)