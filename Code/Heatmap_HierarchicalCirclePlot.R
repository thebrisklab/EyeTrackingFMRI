##########################################################################################################################################################################
##########################################################################################################################################################################
# Heatmap & Hierarchical Circle Plot Creation
# This script includes the function to generate the circular dendrogram representing connectivity data with hierarchical clustering of brain regions.
# Results AIM 2 in the poster

Heatmap.plot <- function(mat = zstat.mat.1917203, lower.bound, upper.bound){
  # store the Z matrix in new variable "mat", and index the row and column names
  colnames(mat) <- 1:100
  rownames(mat) <- 1:100
  
  #### to make the mat to be a full matrix instead of lower trangle matrix
  makeSymm <- function(m) {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    return(m)
  }
  mat <- makeSymm(mat)
  
  heatmap <- as.data.frame(as.table(mat))
  colnames(heatmap) <- c("Region1", "Region2", "Zstat")
  heatmap$Region1 <- as.numeric(heatmap$Region1)
  heatmap$Region2 <- as.numeric(heatmap$Region2)
  #heatmap.filter <- subset(heatmap, abs(heatmap$Z_stat) > 1.96) # to only plot some Z_stats meet specific condition
  
  # Heatmap 
  heatmap.plot <- ggplot(data = heatmap, aes(x = Region1, y = Region2)) + # Changed color to fill for gradient coloring in points
    geom_point(aes(color = Zstat),shape = 15, size = 1.4) + # Added 'size' to adjust point size for better visualization
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(lower.bound, upper.bound)) + # Ensure the gradient scales properly
    theme_minimal() +
    labs(x = "Node 1", y = "Node 2", title = "Connectivity")
  
  return(heatmap.plot)
}

##########################################################################################################################################################################
##########################################################################################################################################################################
Hierarchical.Circle.plot <- function(z_statistic_matrix, threshold.value) {

    ##### Prepare parcellation labels
    # Need to read the meta data from cifti data.
    parc <- load_parc("Schaefer_100")
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
    threshold <- quantile(abs(connections_test$Freq), threshold.value)
    connections_test <- connections_test[abs(connections_test$Freq) >= threshold, ]
    
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
    
    # Define color 
    # Extract the unique color combinations
    unique_colors <- unique(parc$meta$cifti$labels[[1]][-1,][, c("Red", "Green", "Blue", "Alpha")])
    # Now, convert the RGB color values to hexadecimal format
    my_colors <- rgb(unique_colors)
    
    # Create a plot
    p <- ggraph(my_graph, layout = 'dendrogram', circular = TRUE) +
      geom_node_text(aes(x = x * 1.3, y = y * 1.31, filter = leaf, label = group, colour = group, angle = angle), size = 3, alpha = 1, fontface = "bold") +
      geom_node_point(aes(filter = leaf, x = x * 1.07, y = y * 1.07, colour = group, size = 2, alpha = 0.2)) +
      scale_colour_manual(values = my_colors) +
      scale_size_continuous(range = c(0.2, 10)) +
      theme_void()
    
    # Add connectivity bundles to the plot
    p.final <- p +  geom_conn_bundle(data =  get_con(from = connections_test$Var1, to = connections_test$Var2, values = connections_test$Freq), 
                                     alpha = 0.8, 
                                     tension = 0.5,
                                     aes(color = values)) +
      scale_edge_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey50")+
      guides( alpha = FALSE, size = FALSE) +
      labs(colour = "Brain Network")
    
    return(p.final)
}
