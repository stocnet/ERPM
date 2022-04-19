library(igraph)
library(RColorBrewer)

# partition
partition <- c(1,1,1,2,2,2,2,3,3,3,4,4,4,4,4,4)
num.nodes <- 16

# attributes for node colors and shapes
colors <- c(0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3)
shapes <- c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1)

# make graph and clusters for igraph
affiliation <- as.matrix(table(data.frame(actor = 1:num.nodes, group= partition)))
adjacency <- affiliation %*% t(affiliation)
diag(adjacency) <- 0
graph <- graph_from_adjacency_matrix(adjacency)
clusters <- make_clusters(graph, membership = partition, modularity = F)

# set colors for groups
group.colors <- brewer.pal(max(partition), name = "Paired")

# Plot
layout <- layout.auto(graph)
plot(clusters, graph, 
     layout=layout, 
     vertex.label=NA, 
     vertex.size=5,
     edge.color = rgb(1,0,0,alpha = 0), # transparent edges
     mark.border = "white", # color for group border
     mark.col = group.colors) # color for group

