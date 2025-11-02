# ==============================================================================
# Fichier : functions_erpm_bip_network.R
# Fonction : partition_to_bipartite_network(), plot_partition_clusters()
# Utilité : Initialise, transforme et affiche une partition en graphe bipartie
# ==============================================================================
if(!exists(".__functions_erpm_bip_network_loaded", envir = .GlobalEnv)){
  
  # -- Bloc minimal: dépendances et symboles requis --
  # Garantit l'accès à network() et à as.matrix.network.adjacency() sans
  # nécessiter library() à chaque sourcing.
  stopifnot(requireNamespace("network", quietly = TRUE),
            requireNamespace("sna",     quietly = TRUE))

  if (!exists("network", mode = "function", inherits = TRUE)) {
    assign("network",
           get("network", envir = asNamespace("network")),
           envir = .GlobalEnv)
  }

  if (!exists("as.matrix.network.adjacency", mode = "function", inherits = TRUE)) {
    assign("as.matrix.network.adjacency",
           get("as.matrix.network.adjacency", envir = asNamespace("sna")),
           envir = .GlobalEnv)
  }

  if (!exists("VERBOSE", envir = .GlobalEnv)) assign("VERBOSE", FALSE, envir = .GlobalEnv)

  #' Transforme une partition en réseau bipartite
  #'
  #' Cette fonction prend un vecteur d'objets, un vecteur de partition et éventuellement des attributs
  #' pour construire un objet \code{network} bipartite utilisable avec \pkg{ergm}.
  #'
  #' @param labels Vecteur de caractères représentant les noms des objets.
  #' @param partition Vecteur numérique ou factor indiquant le groupe de chaque objet.
  #' @param attributes Liste optionnelle contenant les attributs des objets (ex: \code{gender}, \code{age}).
  #'
  #' @return Un objet de classe \code{network} représentant le réseau bipartite.
  #' @examples
  #' labels <- c("A", "B", "C", "D")
  #' partition <- c(1, 2, 2, 1)
  #' attrs <- list(gender = c(1, 2, 2, 1))
  #' nw <- partition_to_bipartite_network(labels, partition, attrs)
  partition_to_bipartite_network <- function(labels, partition, attributes = list()) {
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
      for ( i in seq_along(labels) ) {
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

  #' Affiche un réseau bipartite avec les clusters de la partition
  #'
  #' Cette fonction prend un objet \code{network} bipartite et affiche deux graphiques :
  #' \itemize{
  #'   \item le réseau bipartite avec ses arêtes
  #'   \item le réseau des objets coloré selon les clusters
  #' }
  #'
  #' @param net Objet de classe \code{network} construit avec \code{partition_to_bipartite_network}.
  #' @return Aucun retour, la fonction produit uniquement des graphiques.
  #' @examples
  #' plot_partition_clusters(nw)
  plot_partition_clusters <- function(net) {
    library(intergraph)
    library(igraph)
    library(RColorBrewer)

    oldpar <- par(no.readonly = TRUE)               # Sauv. la config. modifiable graphique actuelle
    on.exit(par(oldpar), add = TRUE)                # Action à effectuer à la fin de la fonction
    par(mfrow=c(1,2))

    group_indices <- which(grepl("^G", network.vertex.names(net)))
    standard_indices <- setdiff(1:network.size(net), group_indices)

    coords <- matrix(NA, nrow = network.size(net), ncol = 2)
    coords[standard_indices, 1] <- seq(-1, 1, length.out = length(standard_indices))
    coords[group_indices, 1] <- seq(-1, 1, length.out = length(group_indices))
    coords[standard_indices, 2] <- 0
    coords[group_indices, 2] <- 1

    plot(net,
      coord = coords,
      main = "Exemple ERPM",
      cex.main = 0.8,
      label = network.vertex.names(net),
      label.cex = 0.8,
      pad = 3)

    if (!is.null(network::get.network.attribute(net, "bipartite"))) {
      network::delete.network.attribute(net, "bipartite")
    }

    g <- intergraph::asIgraph(net)
    V(g)$name <- V(g)$'vertex.names'

    all_labels <- network.vertex.names(net)
    object_labels <- V(g)$name[!grepl("^G", V(g)$name)]
    group_nodes <- V(g)[grepl("^G", V(g)$name)]
    g <- delete_vertices(g, group_nodes)

    adj <- as.matrix.network.adjacency(net)

    partition <- vapply(object_labels, function(lbl) {
      neighbors <- names(which(adj[lbl, ] == 1))
      group <- neighbors[grepl("^G", neighbors)]
      if (length(group) == 0) {
        return(NA_real_)
      } else {
        return(as.numeric(sub("G", "", group[1])))
      }
    }, FUN.VALUE = numeric(1))

    na_nodes <- object_labels[is.na(partition)]
    
    # Remplacer l'affichage console par log
    # if (exists("log_msg")) {
    #   log_msg("INFO", paste("Partition :", paste(partition, collapse = ", ")))
    #   if (length(na_nodes) > 0) log_msg("WARN", paste("Objets sans groupe :", paste(na_nodes, collapse = ", ")))
    # }

    g_obj <- induced_subgraph(g, vids = object_labels)

    membership <- partition
    names(membership) <- object_labels
    clusters <- make_clusters(g_obj, membership = as.numeric(membership))

    if (all(is.na(partition))) stop("Erreur : aucun groupe n’a été trouvé pour les objets.")

    group.colors <- brewer.pal(max(partition, na.rm = TRUE), "Paired")
    layout <- layout_with_fr(g)

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


  #' Crée un exemple de réseau ERPM bipartite
  #'
  #' Génère un petit réseau de test avec partitions et attributs prédéfinis.
  #'
  #' @return Un objet \code{network} ERPM bipartite prêt à l'emploi.
  #' @examples
  #' nw <- create_erpm_network()
  create_erpm_network <- function() {
    # Exemple de base
    labels <- c("A", "B", "C", "D", "E", "F")
    partition <- c(1, 2, 2, 3, 3, 3)

    attributes <- list(
      gender = c(1, 1, 1, 2, 2, 1),
      age = c(15, 22, 22, 40, 30, 30)#,
      # friendship_matrices <- c(
      #   matrix(c(0, 1, 1, 1, 0, 0), nrow=1),
      #   matrix(c(1, 0, 0, 0, 1, 0), nrow=1),
      #   matrix(c(1, 0, 0, 0, 1, 0), nrow=1),
      #   matrix(c(1, 0, 0, 0, 0, 0), nrow=1),
      #   matrix(c(0, 1, 1, 0, 0, 1), nrow=1),
      #   matrix(c(0, 0, 0, 0, 1, 0), nrow=1))
    )

    if (VERBOSE) cat("\nConstruction du réseau bipartite...\n")

    nw <- partition_to_bipartite_network(labels, partition, attributes)

    if (VERBOSE) cat(green("Réseau bipartite construit avec succès.\n\n"))

    # Retourne le réseau et la partition
    return(list(
      network = nw,
      partition = partition
    ))
    
  }

  #' Log et affiche les informations d'un réseau ERPM
  #'
  #' @param res Liste contenant `network` et `partition` (retour de create_erpm_network)
  #' @param verbose Booléen, si TRUE affiche les infos sur la console (par défaut : VERBOSE global)
  log_erpm_network <- function(res, verbose = get0("VERBOSE", ifnotfound = FALSE)) {
  
    nw <- res$network
    partition <- res$partition
    na_nodes <- names(which(is.na(partition)))

    # Infos formatées
    partition_info <- paste("Partition :", paste(partition, collapse = ", "))
    na_nodes_info  <- if (length(na_nodes) > 0) paste("Objets sans groupe :", paste(na_nodes, collapse = ", ")) else NULL
    nw_info        <- paste0("Réseau ERPM créé : ", network.size(nw), " noeuds, ", network.edgecount(nw), " arêtes")
    vertex_list    <- paste("Liste des sommets :", paste(network.vertex.names(nw), collapse = ", "))

    # Log et affichage
    if (exists("log_msg")) {
      log_msg("INFO", partition_info)
      log_msg("INFO", nw_info)
      log_msg("INFO", vertex_list)
      if (!is.null(na_nodes_info)) log_msg("WARN", na_nodes_info)
    } else if (verbose) {
      cat(cyan("INFO : Réseau ERPM créé\n"))
      cat(partition_info, "\n")
      cat(nw_info, "\n")
      cat(vertex_list, "\n")
      if (!is.null(na_nodes_info)) cat(na_nodes_info, "\n")
    } else {
      if (!is.null(na_nodes_info)) cat(na_nodes_info, "\n")
    }
    
    invisible(NULL)
  } 

  assign("create_erpm_network",                   create_erpm_network,          envir = .GlobalEnv)
  assign("partition_to_bipartite_network",          partition_to_bipartite_network, envir = .GlobalEnv)
  assign("plot_partition_clusters",               plot_partition_clusters,      envir = .GlobalEnv)
  assign("log_erpm_network",                      log_erpm_network,             envir = .GlobalEnv)
  assign(".__functions_erpm_bip_network_loaded",  TRUE,                         envir = .GlobalEnv)

}
