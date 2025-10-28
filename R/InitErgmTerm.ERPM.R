# InitErgmTerm.groups<-function(nw, arglist, ..., version=packageVersion("ergm")) {
#     ### Check the network and arguments to make sure they are appropriate.
#     # message("Start using ERGM \"groups\" Init Ergm Term")
#     termname <- "groups"
#     coefpre <- "groups_size"
#   a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=FALSE,
#                         varnames = c("from", "to", "pow"),
#                         vartypes = c("numeric","numeric", "numeric"),
#                         defaultvalues = list(NULL,Inf, 2),
#                         required = c(TRUE,FALSE, FALSE))
#     ### Process the arguments
#   from <- a$from; to <- a$to; pow <- a$pow;

#   if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
#   else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
#   else if(length(from)!=length(to)) ergm_Init_stop("The arguments of term ", sQuote(termname), " must have arguments either of the same length, or one of them must have length 1.")

#   to <- ifelse(to==Inf, pmax(from, network.size(nw)) + 1, to)
#   if(any(from>=to)) ergm_Init_stop("Term ", sQuote(termname), " must have ", sQuote("from"), "<", sQuote("to"), ".")

#   emptynwstats<-NULL

#   name <- "sgroups"
#   coefpre <- "groups_square"
#   if (any(from==0)) {
#       emptynwstats <- rep(0, length(from))
#       emptynwstats[from==0] <- network.size(nw)
#   }

#   if(length(from)==0){return(NULL)}
#   coef.names <- ifelse(to>=network.size(nw)+1,
#                   paste(coefpre, from,"+",if(a$pow!=1) a$pow else "",sep=""),
#                   paste(coefpre, from,"to",to,"pow",if(a$pow!=1) a$pow else "",sep=""))

#   inputs <- c(rbind(from,to,pow))

#   ### Construct the list to return
#     list(name=name,coef.names=coef.names, inputs=inputs, dependence=TRUE, minval = 0, maxval=network.size(nw), conflicts.constraints=paste0("degreedist"), emptynwstats=emptynwstats)
# }

# InitErgmTerm.squared_sizes<-function(nw, arglist, ..., version=packageVersion("ergm")) {
#     ### Check the network and arguments to make sure they are appropriate.
#     message("Start using CUBI \"groups\" Init Ergm Term")
#     termname <- "squared_sizes"
#     coefpre <- "squared_sizes"
#   a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=TRUE,
#                         varnames = c("from", "to", "pow"),
#                         vartypes = c("numeric","numeric", "numeric"),
#                         defaultvalues = list(NULL,Inf, 2),
#                         required = c(TRUE,FALSE, FALSE))
#     ### Process the arguments
#   from <- a$from; to <- a$to; pow <- a$pow;

#   if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
#   else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
#   else if(length(from)!=length(to)) ergm_Init_stop("The arguments of term ", sQuote(termname), " must have arguments either of the same length, or one of them must have length 1.")

#   to <- ifelse(to==Inf, pmax(from, network.size(nw)) + 1, to)
#   if(any(from>=to)) ergm_Init_stop("Term ", sQuote(termname), " must have ", sQuote("from"), "<", sQuote("to"), ".")

#   emptynwstats<-NULL

#   name <- "squared_sizes"
#   coefpre <- "squared_sizes"
#   if (any(from==0)) {
#       emptynwstats <- rep(0, length(from))
#       emptynwstats[from==0] <- network.size(nw) - nw %n% "bipartite"
#   }

#   if(length(from)==0){return(NULL)}
#   coef.names <- ifelse(to>=network.size(nw)+1,
#                   paste(coefpre, from,"+",if(a$pow!=1) a$pow else "",sep=""),
#                   paste(coefpre, from,"to",to,"pow",if(a$pow!=1) a$pow else "",sep=""))

#   inputs <- c(rbind(from,to,pow))

#   ### Construct the list to return
#     list(name=name,coef.names=coef.names, inputs=inputs, dependence=TRUE, minval = 0, maxval=network.size(nw), conflicts.constraints=paste0("b2","degreedist"), emptynwstats=emptynwstats)
# }

# InitErgmTerm.cliques <- function(nw, arglist, ..., version=packageVersion("ergm")) {
#   ### Check the network and arguments to make sure they are appropriate.
#   message("Start using \"cliques\" Init Ergm Term")
#   a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=TRUE,
#                       varnames = c("nb_cliques"),
#                       vartypes = c("numeric"),
#                       defaultvalues = list(2),
#                       required = c(FALSE))
#   ### Construct the list to return
#   message("Return list for cliques")
#   list(name="cliques",                               #name: required
#        coef.names = "cliques",                       #coef.names: required
#        emptynwstats = network.size(nw), # When nw is empty, isolates=n, not 0,
#        minval = 0,
#        maxval = network.size(nw),
#        conflicts.constraints="degreedist"
#        )
# }