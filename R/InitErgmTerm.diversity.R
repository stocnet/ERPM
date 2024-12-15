#' @templateVar name nodecovrange
#' @title Range of covariate values for neighbors of a node
#' @description This term adds a single network statistic equalling
#'   the sum over the nodes of the difference between the highest
#'   value of a nodal covariate and its lower covariate.
#'
#' @usage
#' # binary: nodecovrange(attr)
#'
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept quantitative nodal attribute
InitErgmTerm.nodecovrange<-function (nw, arglist, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                      varnames = c("attr"),
                      vartypes = c(ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
### Process the arguments
  nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric")
  coef.names <- nodecov_names(nodecov, "nodecovrange")
  list(name="nodecovrange", coef.names=coef.names, inputs=c(nodecov))
}


#' @templateVar name nodeocovrange
#' @title Range of covariate values for out-neighbors of a node
#' @description This term adds a single network statistic equalling
#'   the sum over the nodes of the difference between the highest
#'   value of a nodal covariate and its lower covariate.
#'
#' @usage
#' # binary: nodeocovrange(attr)
#'
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept quantitative nodal attribute
InitErgmTerm.nodeocovrange<-function (nw, arglist, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                      varnames = c("attr"),
                      vartypes = c(ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
### Process the arguments
  nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric")
  coef.names <- nodecov_names(nodecov, "nodeocovrange")
  list(name="nodeocovrange", coef.names=coef.names, inputs=c(nodecov))
}


#' @templateVar name nodeicovrange
#' @title Range of covariate values for in-neighbors of a node
#' @description This term adds a single network statistic equalling
#'   the sum over the nodes of the difference between the highest
#'   value of a nodal covariate and its lower covariate.
#'
#' @usage
#' # binary: nodeicovrange(attr)
#'
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept quantitative nodal attribute
InitErgmTerm.nodeicovrange<-function (nw, arglist, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                      varnames = c("attr"),
                      vartypes = c(ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
### Process the arguments
  nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric")
  coef.names <- nodecov_names(nodecov, "nodeicovrange")
  list(name="nodeicovrange", coef.names=coef.names, inputs=c(nodecov))
}


#' @templateVar name b1covrange
#' @title Range of covariate values for neighbors of a mode-1 node
#' @description This term adds a single network statistic equalling
#'   the sum over the nodes of the difference between the highest
#'   value of a nodal covariate and its lower covariate.
#'
#' @usage
#' # binary: nodecovrange(attr)
#'
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept quantitative nodal attribute
InitErgmTerm.b1covrange<-function (nw, arglist, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=TRUE,
                      varnames = c("attr"),
                      vartypes = c(ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
### Process the arguments
  nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric", bip="b2")
  coef.names <- nodecov_names(nodecov, "b1covrange")
  list(name="b1covrange", coef.names=coef.names, inputs=c(nodecov))
}



#' @templateVar name b2covrange
#' @title Range of covariate values for neighbors of a mode-2 node
#' @description This term adds a single network statistic equalling
#'   the sum over the nodes of the difference between the highest
#'   value of a nodal covariate and its lower covariate.
#'
#' @usage
#' # binary: nodecovrange(attr)
#'
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept quantitative nodal attribute
InitErgmTerm.b2covrange<-function (nw, arglist, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=TRUE,
                      varnames = c("attr"),
                      vartypes = c(ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
### Process the arguments
  nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric", bip="b1")
  coef.names <- nodecov_names(nodecov, "b2covrange")
  list(name="nodeicovrange", coef.names=coef.names, inputs=c(nodecov))
}



#' @templateVar name nodefactordistinct
#' @title Number of distinct neighbor types
#' @description This term adds multiple network statistics to the
#'   model, one for each of (a subset of) the unique values of the
#'   `attr` attribute (or each combination of the attributes
#'   given). Each of these statistics gives the number of times a node
#'   with that attribute or those attributes appears in an edge in the
#'   network.
#'   
#' @usage
#' # binary: nodefactordistinct(attr, levels=-1)
#'
#' @template ergmTerm-attr
#' @templateVar explain this optional argument controls which levels of the attribute
#'   should be included and which should be excluded.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-levels-not-first
#'
#' @template ergmTerm-base-dep
#'
#' @concept directed
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.nodefactordistinct<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL),
                        required = c(TRUE, FALSE),
                        dep.inform = list(FALSE, FALSE))
  attrarg <- a$attr                        
  levels <- a$levels    

  nodecov <- ergm_get_vattr(attrarg, nw)
  attrname <- attr(nodecov, "name")
  u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))


  if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
    return()
  } 
  #   Recode to numeric
  nodepos <- match(nodecov,u,nomatch=0)
  ### Construct the list to return
  inputs <- c(max(nodepos), nodepos)
  list(name="nodefactordistinct",                                        #required
       coef.names = paste("nodefactordistinct", paste(attrname,collapse="."), sep="."), #required
       iinputs = inputs,
       minval = 0
       )
}


# ------TEST----------
  
#' @templateVar name b1factordistinct
#' @title Number of distinct neighbor types for second mode nodes in bipartite 
#'  networks
#' @description This term adds multiple network statistics to the
#'   model, one for each of (a subset of) the unique values of the
#'   `attr` attribute (or each combination of the attributes
#'   given). Each of these statistics gives the number of times a node
#'   with that attribute or those attributes appears in an edge in the
#'   network.
#'   
#' @usage
#' # binary: b1factordistinct(attr, levels=-1)
#'
#' @template ergmTerm-attr
#' @templateVar explain this optional argument controls which levels of the attribute
#'   should be included and which should be excluded.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-levels-not-first
#'
#' @template ergmTerm-base-dep
#'
#' @concept directed
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.b1factordistinct<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "levels"),
                      vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE),
                      dep.inform = list(FALSE, FALSE))
  attrarg <- a$attr                        
  levels <- a$levels    
  
  nodecov <- ergm_get_vattr(attrarg, nw, bip="b1")
  attrname <- attr(nodecov, "name")
  u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
  
  if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
    return()
  } 
  #   Recode to numeric
  nodepos <- match(nodecov,u,nomatch=0)
  ### Construct the list to return
  inputs <- c(max(nodepos), nodepos)
  list(name="b1factordistinct",                                        #required
       coef.names = paste("b1factordistinct", paste(attrname,collapse="."), sep="."), #required
       iinputs = inputs,
       minval = 0
  )
}
