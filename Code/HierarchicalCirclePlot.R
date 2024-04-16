##########################################################################################################################################################################
##########################################################################################################################################################################
# Hierarchical Circle Plot Creation
# This script generates a circular dendrogram representing connectivity data with hierarchical clustering of brain regions.

# Load necessary libraries
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)

# Load statistical matrix data (V by V matrix)
z_statistic_matrix <- mean_z_stat_mat_td

##### Prepare parcellation labels
parc_cluster <- parc$meta$cifti$labels[[1]][-1, ] %>% dplyr::select(Key)
parc_cluster$Region <- rownames(parc_cluster)
parc_cluster$Group <- sub("^([^_]+)_([^_]+)_([^_]+)_.*$", "\\3", parc_cluster$Region)

# Define the hierarchical structure
d1_test <- data.frame(from = "origin", to = unique(parc_cluster$Group))
d2_test <- data.frame(from = parc_cluster$Group, to = parc_cluster$Key)
edges_test <- rbind(d1_test, d2_test)

# Convert z-statistic matrix for edge creation
colnames(z_statistic_matrix) <- 1:100
rownames(z_statistic_matrix) <- 1:100
edges_df <- as.data.frame(as.table(z_statistic_matrix))
connections_test <- edges_df
connections_test$Freq <- as.numeric(as.character(connections_test$Freq))
connections_test <- subset(connections_test, connections_test$Freq != 0)

# Need to change the threshold based on specifc purpose, i.e., to keep the "connections" you wanted to show in the plot
#threshold <- quantile(abs(connections_test$Freq), 0.99)
#connections_test <- connections_test[abs(connections_test$Freq) >= threshold, ]

# Define vertices for the graph
vertices_test <- data.frame(
  name = unique(c(as.character(edges_test$from), as.character(edges_test$to)))
)
vertices_test$group <- edges_test$from[match(vertices_test$name, edges_test$to)]

# Ordering groups in vertices
group_order <- unique(parc_cluster$Group)
vertices_test$group <- factor(vertices_test$group, levels = group_order)
vertices_test <- vertices_test[order(vertices_test$group), ]

# Calculate label angles
vertices_test$id <- NA
myleaves <- which(is.na(match(vertices_test$name, edges_test$from)))
nleaves <- length(myleaves)
vertices_test$id[myleaves] <- seq(1:nleaves)
vertices_test$angle <- 90 - 360 * vertices_test$id / nleaves
vertices_test$hjust <- ifelse(vertices_test$angle < -90, 1, 0)
vertices_test$angle <- ifelse(vertices_test$angle < -90, vertices_test$angle + 180, vertices_test$angle)


##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

# Create graph object
my_graph <- igraph::graph_from_data_frame(edges_test, vertices = vertices_test)

# Create a plot
p <- ggraph(my_graph, layout = 'dendrogram', circular = TRUE) +
  geom_node_text(aes(x = x * 1.3, y = y * 1.31, filter = leaf, label = group, colour = group, angle = angle), size = 3, alpha = 1, fontface = "bold") +
  geom_node_point(aes(filter = leaf, x = x * 1.07, y = y * 1.07, colour = group, size = 2, alpha = 0.2)) +
  scale_colour_manual(values = my_colors) +
  scale_size_continuous(range = c(0.2, 10)) +
  theme_void()

# Add connectivity bundles to the plot
p + geom_conn_bundle(data = connections_test, alpha = 0.8, tension = 0.5, aes(color = Freq)) +
  scale_edge_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey50") +
  guides(alpha = FALSE, size = FALSE) +
  labs(colour = "Brain Network")
