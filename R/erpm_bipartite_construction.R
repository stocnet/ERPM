build_bipartite_network <- function(labels, partition, attributes = list()) {
    stopifnot(length(labels) == length(partition))

    # Create all possible groups for partition, event the ones not used
    group_labels <- paste0("G", seq_along(labels))
    all_labels <- c(labels, group_labels)
    n <- length(all_labels)
  
    # Create an empty adjacency matrix
    adj <- matrix(0, nrow = n, ncol = n)
    rownames(adj) <- all_labels
    colnames(adj) <- all_labels
  
    # connect objects to their respective groups from the partition
    for (i in seq_along(labels)) {
        object_label <- labels[i]
        group_label <- paste0("G", partition[i])
        adj[object_label, group_label] <- 1
        adj[group_label, object_label] <- 1
    }
  
    # Building the bipartite network
    net <- network(adj, directed = FALSE, matrix.type = "adjacency")
    network::network.vertex.names(net) <- all_labels
    network::set.network.attribute(net, "bipartite", length(labels))
  
    # add the attributes to the network
    for (attr_name in names(attributes)) {
        full_attr <- c(attributes[[attr_name]], rep(NA, length(group_labels)))
        network::set.vertex.attribute(net, attr_name, full_attr)
    }
    return(net)
}

plot_partition_clusters <- function(net) {
  library(intergraph)
  library(igraph)
  library(RColorBrewer)

  # create a 2 column layout for the plots
  par(mfrow=c(1,2))

  # Identify groups and standard nodes
  group_indices <- which(grepl("^G", network.vertex.names(net)))
  standard_indices <- setdiff(1:network.size(net), group_indices)

  # Define coordinates for plotting a bipartite network
  coords <- matrix(NA, nrow = network.size(net), ncol = 2)

  # X position : fix for everyone
  coords[standard_indices, 1] <- seq(-1, 1, length.out = length(standard_indices))
  coords[group_indices, 1] <- seq(-1, 1, length.out = length(group_indices))

  # Y position : 0 for standard nodes, 1 for group nodes
  coords[standard_indices, 2] <- 0
  coords[group_indices, 2] <- 1

  # Plot first graph : bipartite network with edges and nodes
  plot(net,
     coord = coords,
     main = "Exemple ERPM",
     cex.main = 0.8,
     label = network.vertex.names(net),
     label.cex = 0.8,
     pad = 3)

  # Delete the bipartite attribute for the next plot
  if (!is.null(network::get.network.attribute(net, "bipartite"))) {
    network::delete.network.attribute(net, "bipartite")
  }

  # Graph conversion
  g <- intergraph::asIgraph(net)
  V(g)$name <- V(g)$'vertex.names'

  # Labels
  all_labels <- network.vertex.names(net)
  object_labels <- V(g)$name[!grepl("^G", V(g)$name)]
  group_nodes <- V(g)[grepl("^G", V(g)$name)]
  g <- delete_vertices(g, group_nodes)

  # Adjacency matrix
  adj <- as.matrix.network.adjacency(net)

  # Find the partition gourps for each object
  partition <- vapply(object_labels, function(lbl) {
  neighbors <- names(which(adj[lbl, ] == 1))
  group <- neighbors[grepl("^G", neighbors)]
  if (length(group) == 0) {
    warning(paste("No group found for node:", lbl))
    return(NA_real_)
  } else {
    return(as.numeric(sub("G", "", group[1])))
  }
}, FUN.VALUE = numeric(1))

na_nodes <- object_labels[is.na(partition)]
cat("Partition :", partition, "\n")
cat("na_nodes :", na_nodes, "\n")
if (length(na_nodes) > 0) {
  warning("Les objets suivants ne sont connectés à aucun groupe : ", paste(na_nodes, collapse = ", "))
}

  # Subgraph for objects
  g_obj <- induced_subgraph(g, vids = object_labels)

  # Define clusters
  membership <- partition
  names(membership) <- object_labels
  clusters <- make_clusters(g_obj, membership = as.numeric(membership))

  # coloring
  if (all(is.na(partition))) {
  stop("Erreur : aucun groupe n’a été trouvé pour les objets. Vérifie les connexions entre objets et sommets de type 'Gx'.")
}
  group.colors <- brewer.pal(max(partition, na.rm = TRUE), "Paired")

  # Layout
  layout <- layout_with_fr(g)

  # Plot the network as groups
  plot(clusters, g,
       layout = layout,
       vertex.label = V(g)$name,
       vertex.label.cex = 0.6,
       vertex.size = 6,
       edge.color = rgb(0.5, 0.5, 0.5, 0.3),
       mark.border ="white",
       mark.col = group.colors
  )
}
