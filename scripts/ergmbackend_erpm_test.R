dev_mode(TRUE)
library(ergm)


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
partition <- c(1,1,2,2,2,3)

nw <- network.initialize(n*2, dir=FALSE, bip=n)
for(a in names(nodes)) set.vertex.attribute(nw, a, c(nodes[[a]], rep(NA, n)))
nw %n% "friendship" <- friendship
nw[cbind(seq_along(partition), partition+n)] <- 1

as.matrix(nw)

summary(nw~Proj1(~B(~edges + nodematch("gender") + edgecov("friendship"), "nonzero")))

fit <- ergm(nw~b2degrange(1,Inf) + Proj1(~B(~nodematch("gender") + absdiff("age") + edgecov("friendship"), "nonzero")), constraints = ~b1part)

fit2 <- ergm(nw~b2degrange(1,Inf) + Proj1(~B(~nodematch("gender") + absdiff("age"), "nonzero")), constraints = ~b1part)

fit <- ergm(nw~b2degrange(1,Inf) + offset(b2degrange(5,Inf)) + Proj1(~B(~nodematch("gender") + absdiff("age") + edgecov("friendship"), "nonzero")), offset.coef = -Inf, constraints = ~b1part)

