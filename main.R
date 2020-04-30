######################################################################
## Simulation and estimation of Exponential Random Partition Models ## 
## Example script with random data                                  ##
##                                                                  ## 
## Author: Marion Hoffman                                           ##
######################################################################

# Libraries required to run the ERPM functions
library(ggplot2)
library(numbers)
library(combinat)
library(gmp)


# Source files for the ERPM functions
current_wd <- getwd()
filesLoad = list.files(paste(current_wd,"/functions_ERPM",sep=""), pattern = "R$")
lapply(filesLoad, function(x) source(paste(current_wd,"/functions_ERPM/",x,sep="")))


############################## SIMULATION ############################

## ------------------------- DEFINE OBJECTS --------------------------
# Here we define an arbitrary set of n=6 nodes with attributes, and an 
# arbitrary covariate matrix
n <- 6
nodes <- data.frame(label = c("A","B","C","D","E","F"),
                    gender = c(1,1,2,1,2,2),
                    age = c(20,22,25,30,30,31)) 
friendship <- matrix(c(0, 1, 1, 1, 0, 0,
                       1, 0, 0, 0, 1, 0,
                       1, 0, 0, 0, 1, 0,
                       1, 0, 0, 0, 0, 0,
                       0, 1, 1, 0, 0, 1,
                       0, 0, 0, 0, 1, 0), 6, 6, TRUE) 


## --------------------- EXACT CALCULATIONS --------------------------
# For small values of n, we can calculate exactly several things,
# first we can compute the number of possible partitions
# (this is the Bell number):
number_partitions <- bell(n)
# or the number f partitions with k=2 groups for example
number_partitions_2groups <- Stirling2(n,2)
# or the number of 

# we can also compute the average size (number of groups) of a partition
# (under a null model)
average_size_nullmodel <- compute_averagesize(n)



##################### ESTIMATION ############################

## --------------------- DEFINE OBJECTS ---------------------

load("data/RFID_interactions.RData")

nodes <- data.frame(label = c("A","B","C","D","E","F"),
                    gender = c(1,1,2,1,2,2),
                    age = c(20,22,25,30,30,31))
friendship <- matrix(c(0, 1, 1, 1, 0, 0,
                       1, 0, 0, 0, 1, 0,
                       1, 0, 0, 0, 1, 0,
                       1, 0, 0, 0, 0, 0,
                       0, 1, 1, 0, 0, 1,
                       0, 0, 0, 0, 1, 0), 6, 6, TRUE)


## ----------------- CROSS-SECTIONAL VERSION -----------------

# Estimation
partition <- c(1,2,2,1,3,3)
effects <- list( names = c("size","diff"),
                 objects = c("partition","age"))
objects <- list(friendship = friendship)
startingestimates <- rep(0,length(effects$names))
multiplicationfactor <- 30
gainfactor <- 0.1 
mini.steps <- "normalized"
burnin <- 100
thining <- 30
length.p1 <- 100
min.iter.p2 <- 100
max.iter.p2 <- 200
num.steps.p2 <- 6
length.p3 <- 1000
neighborhood <- 1

results.erpm <- estimate_ERPM(partition, nodes, objects, effects, 
                              startingestimates, multiplicationfactor, gainfactor, mini.steps, burnin, thining,
                              length.p1, min.iter.p2, max.iter.p2, num.steps.p2, length.p3, neighborhood)


## --------------------- EVENTS VERSION ---------------------

firstPartition <- c(1,2,3,3,4,4)
events <- list()
events[[1]] <- list(partition = c(1, 1, 2, 2, 3, 3), time = 1,  old_g1 = 1, old_g2 = 2, new_g1 = 1, new_g2 = 0)
events[[2]] <- list(partition = c(1, 1, 1, 1, 2, 2), time = 5,  old_g1 = 1, old_g2 = 2, new_g1 = 1, new_g2 = 0)
events[[3]] <- list(partition = c(1, 2, 1, 2, 3, 3), time = 10, old_g1 = 1, old_g2 = 0, new_g1 = 1, new_g2 = 2)
events[[4]] <- list(partition = c(1, 2, 1, 3, 4, 4), time = 12, old_g1 = 2, old_g2 = 0, new_g1 = 2, new_g2 = 3)
events[[5]] <- list(partition = c(1, 2, 1, 3, 2, 2), time = 20, old_g1 = 2, old_g2 = 4, new_g1 = 2, new_g2 = 0)
events[[6]] <- list(partition = c(1, 2, 1, 3, 4, 4), time = 25, old_g1 = 2, old_g2 = 0, new_g1 = 2, new_g2 = 4)
effects <- list( names = c("diff"),
                 objects = c("age"))
objects <- list(friendship = friendship)
intercept <- TRUE

startingestimates <- rep(0,length(effects$names))
if(intercept) startingestimates <- rep(0,length(effects$names)+1)
initialDamping <- 1
maxIterations <- 20 
dampingIncreaseFactor <- 2
dampingDecreaseFactor <- 3
maxScoreStopCriterion <- 0.001

neighborhood <- 1

results.pem <- estimate_PEM(firstPartition, events, nodes, objects, effects, intercept,
                            startingestimates, maxIterations, dampingIncreaseFactor, dampingDecreaseFactor, maxScoreStopCriterion, neighborhood)


# TODO:
# find a better example to avoid collinearity problems


##################### SIMULATION ############################

## --------------------- DEFINE OBJECTS ---------------------

nodes <- data.frame(label = paste("Actor",1:100),
                    att = sample(c(0,1), replace=TRUE, size=100))
nodes <- data.frame(label = paste("Actor",1:10),
                    att = sample(c(0,1), replace=TRUE, size=10))


## ----------------- CROSS-SECTIONAL VERSION -----------------

# Trace plots
effects <- list( names = c("num_groups"),#,"size2","diff","alter"),
                 objects = c("partition"))#,"partition","att","att"))
objects <- list()
minestimates <- c(-5)#,-1,-1,1)
maxestimates <- c()
mini.steps <- "normalized"
burnin <- 0
thining <-500
num.steps <- 500
neighborhood <- 1

traces <- simulate_ERPM(nodes, objects, effects, minestimates, maxestimates, mini.steps, burnin, thining, num.steps, neighborhood) 


traces <- data.frame(x=traces[,1], y1=traces[,2])#, y2=traces[,3], y3=traces[,4], y4=traces[,5])
ggplot(traces, aes(x)) +                    
  geom_line(aes(y=y1), colour="red")# +  
#  geom_line(aes(y=y2), colour="green") +
#  geom_line(aes(y=y3), colour="blue") +  
#  geom_line(aes(y=y4), colour="orange") 

# Parameters evolution: size
effects <- list( names = c("num_groups"),
                 objects = c("partition"))
objects <- list()
minestimates <- c(-15)
maxestimates <- c(10)
burnin <- 500
thining <- 1000
num.steps <- 500
neighborhood <- 1

plots <- simulate_ERPM(nodes, objects, effects, minestimates, maxestimates, mini.steps, burnin, thining, num.steps, neighborhood) 

par(mfrow=c(1,1))
plot(plots[[1]]$x,plots[[1]]$y,type="o",main="number of groups")

par(mfrow=c(2,2))
plot(plots[[1]]$x,plots[[1]]$y,type="o",main="size")
plot(plots[[2]]$x,plots[[2]]$y,type="o",main="size2")
plot(plots[[3]]$x,plots[[3]]$y,type="o",main="diff")
plot(plots[[4]]$x,plots[[4]]$y,type="o",main="alter")

# estimate parameter with one fixed
# Estimation
partition <- c(1,1,1,2,2,3)
effects <- list( names = c("num_groups","sizes_squared"),
                 objects = c("partition","partition"))
objects <- list()
startingestimates <- c(0,0.4)
multiplicationfactor <- 30
gainfactor <- 0.1 
mini.steps <- "normalized"
burnin <- 100
thining <- 100
length.p1 <- 500
min.iter.p2 <- 100
max.iter.p2 <- 200
num.steps.p2 <- 15
length.p3 <- 500
neighborhood <- 1
fixed.estimates <- list()
fixed.estimates[[1]] <- NULL
fixed.estimates[[2]] <- 0.4

results.erpm <- estimate_ERPM(partition, nodes, objects, effects, 
                              startingestimates, multiplicationfactor, gainfactor, mini.steps, burnin, thining,
                              length.p1, min.iter.p2, max.iter.p2, num.steps.p2, length.p3, neighborhood, fixed.estimates)


## --------------------- EVENTS VERSION ---------------------

firstPartition <- c(1,2,3,3,4,4)
events <- list()
events[[1]] <- list(partition = c(1, 1, 2, 2, 3, 3), time = 1,  old_g1 = 1, old_g2 = 2, new_g1 = 1, new_g2 = 0)
events[[2]] <- list(partition = c(1, 1, 1, 1, 2, 2), time = 5,  old_g1 = 1, old_g2 = 2, new_g1 = 1, new_g2 = 0)
events[[3]] <- list(partition = c(1, 2, 1, 2, 3, 3), time = 10, old_g1 = 1, old_g2 = 0, new_g1 = 1, new_g2 = 2)
events[[4]] <- list(partition = c(1, 2, 1, 3, 4, 4), time = 12, old_g1 = 2, old_g2 = 0, new_g1 = 2, new_g2 = 3)
events[[5]] <- list(partition = c(1, 2, 1, 3, 2, 2), time = 20, old_g1 = 2, old_g2 = 4, new_g1 = 2, new_g2 = 0)
events[[6]] <- list(partition = c(1, 2, 1, 3, 4, 4), time = 25, old_g1 = 2, old_g2 = 0, new_g1 = 2, new_g2 = 4)
effects <- list( names = c("diff"),
                 objects = c("age"))
objects <- list(friendship = friendship)
intercept <- TRUE

startingestimates <- rep(0,length(effects$names))
if(intercept) startingestimates <- rep(0,length(effects$names)+1)
initialDamping <- 1
maxIterations <- 20 
dampingIncreaseFactor <- 2
dampingDecreaseFactor <- 3
maxScoreStopCriterion <- 0.001

neighborhood <- 1

results.pem <- estimate_PEM(firstPartition, events, nodes, objects, effects, intercept,
                            startingestimates, maxIterations, dampingIncreaseFactor, dampingDecreaseFactor, maxScoreStopCriterion, neighborhood)






### Future maybe interface

# define objects
actors <- defineNodes()
partition.events <- defineEvents()
interaction.groups <- definePartition(nodes =, initialPartition =)
interaction.groups <- linkEvents(partition = interaction.groups, events = partition.events)
friendship <- defineNetwork(network =, nodes = )

# PEM estimation
results.pem <- estimate.pem(partition.events ~ inertia + recip + same())

# ERPM estimation
results.erpm <- estimate.erpm(interaction.groups ~ inertia + recip + same())