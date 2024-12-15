dev_mode(TRUE)
library(ergm)

covrange <- function(nw, a, which = c("u", "i", "o", "1", "2")){
  which <- match.arg(which)
  m <- as.matrix(nw)
  m[m==0] <- NA
  w <- switch(which,
              u =,
              i =,
              o = nw %v% a,
              `1` = (nw %v% a)[-seq_len(nw%n%"bipartite")],
              `2` = (nw %v% a)[seq_len(nw%n%"bipartite")])

  if(which %in% c("o","1")) m <- t(m)

  r <- apply(m*w, 2, range, na.rm = TRUE) |> diff()
  sum(r[!is.infinite(r)])
}

data(florentine)
covrange(flomarriage, "wealth", "u")
summary(flomarriage~nodecovrange("wealth"))

data(davis, package="latentnet")
davis %v% "w" <- rnorm(network.size(davis))
covrange(davis, "w", "1")
summary(davis~b1covrange("w"))

covrange(davis, "w", "2")
summary(davis~b2covrange("w"))

data(sampson)
samplike %v% "w" <- rnorm(network.size(samplike))
covrange(samplike, "w", "o")
summary(samplike~nodeocovrange("w"))


covrange(samplike, "w", "i")
summary(samplike~nodeicovrange("w"))


flomarriage %v% "c" <- sample.int(5, network.size(flomarriage), replace=TRUE)
summary(flomarriage~nodefactordistinct("c", levels = NULL))

m <- as.matrix(flomarriage)
m[m==0] <- NA
c <- flomarriage %v% "c"
apply(m*c,2, unique, simplify = FALSE) |> lapply(na.omit) |> lengths() |> sum()
