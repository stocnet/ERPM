######################################################################
## Simulation and estimation of Exponential Random Partition Models ##
## Visualization of paritions                                       ##
## Author: Alexandra Amani                                          ##
######################################################################


#' Visualization of partition
#'
#'This function plot the groups of a partition
#'
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph make_clusters
#' @importFrom igraph layout.auto
#' @importFrom RColorBrewer brewer.pal
#' @param partition A partition (vector)
#' @param title Character, the title of the plot (default = NULL)
#' @param group.color A vector with the colors of the groups (default = NULL)
#' @param attribute.color A vector, attribute to represent with colors (default = NULL)
#' @param attribute.shape A vector, attribute to represent with shapes (default = NULL)
#' @return A plot of the partition
#' @importFrom grDevices rgb
#' @examples
#' p <- c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4)
#' attr1 <- c(1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 2)
#' attr2 <- c(1, 1, 1, 1, 0, 0, 3, 0, 1, 0, 1, 1, 1, 1, 1, 2)
#' plot_partition(p, attribute.color = attr1, attribute.shape = attr2)
#' @export

plot_partition <- function(partition, title = NULL, group.color = NULL,
                           attribute.color = NULL, attribute.shape = NULL) {

  num.nodes <- length(partition)
  num.group <- table(partition)
  cat.col <- unique(attribute.color)
  possible.color <- brewer.pal(length(cat.col), name = "Reds")
  possible.shape <-  c("circle", "square", "csquare", "rectangle", "crectangle", "vrectangle")

  # Color for Attribute
  if (is.null(attribute.color)) {
    n.color <- "black"
  }else {
    n.color <- attribute.color
    for (i in 1: length(cat.col)) {
      n.color <- replace(n.color, n.color == cat.col[i], possible.color[i])
    }
  }


  # Shape for Attribute
  if (is.null(attribute.shape)) {
    n.shape <- "circle"
  }else {
    n.shape <- attribute.shape
    cat.shape <- unique(attribute.shape)
    for (i in 1:length(cat.shape)) {
      n.shape <- replace(n.shape, n.shape == cat.shape[i], possible.shape[i])
    }
  }


  # set colors for groups
  if (is.null(group.color)) {
    g.colors <- brewer.pal(max(partition), name = "Blues")
  }else {
    g.colors <- group.color
  }

  # make graph and clusters for igraph
  affiliation <- as.matrix(table(data.frame(actor = 1:num.nodes, group = partition)))
  adjacency <- affiliation %*% t(affiliation)
  diag(adjacency) <- 0
  graph <- graph_from_adjacency_matrix(adjacency)
  clusters <- make_clusters(graph, membership = partition, modularity = FALSE)
  layout <- layout.auto(graph)

  # Plot

  plot(clusters, graph,
       layout = layout,
       col = n.color,
       vertex.shape = n.shape,
       vertex.label = NA,
       vertex.size = 5,
       edge.color = rgb(1, 0, 0, alpha = 0), # transparent edges
       mark.border = "white", # color for group border
       mark.col = g.colors,
       main = title) # color for group
}