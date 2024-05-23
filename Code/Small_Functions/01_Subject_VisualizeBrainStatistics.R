# Define a function to map statistical data (e.g., mean signals, coefficients, or t-statistics)
# onto a brain graph. This is often used in neuroimaging to visualize data on a brain model.

# Inputs:
# - dtseries_data: Raw .dtseries.nii MRI data.
# - mapping_data: The statistical data to be mapped onto the brain graph.
#                 Must be of the same length as the number of brain regions.
# - brain_parcellation: The parcellation scheme used for brain region delineation.
# - region_count: The number of brain regions (excluding the 19 subcortex regions).

# Output:
# - A brain graph object with the mapped statistical data.

fMRI_xii_pmean <- function(dtseries_data, brain_parcellation = "Schaefer_100", region_count = 100) {
  # Load the brain parcellation data
  parc <- load_parc(brain_parcellation)
  
  # Convert parcellation data into a vector, indicating the associated brain region
  parc_vec <- c(as.matrix(parc))
  
  # Adjust for non-cortical vertices to align with the parcellation scheme
  adjusted_data <- move_from_mwall(dtseries_data, NA)
  
  # Convert the adjusted data into a matrix format (voxel by time)
  xii_mat <- as.matrix(adjusted_data)
  
  # Retrieve and process the labels for subcortex regions
  subcortex_labels <- adjusted_data$meta$subcort$labels
  
  # Combine cortical and subcortex labels
  combined_labels <- c(as.character(parc_vec), as.character(subcortex_labels))
  
  # Generate labels for brain regions including cortex and subcortex
  labels_vector <- c(as.character(1:region_count), levels(subcortex_labels)[3:21])
  
  # Initialize a matrix to store the computed mean fMRI signal values
  regional_means <- matrix(nrow = region_count + 19, ncol = ncol(xii_mat))
  
  # Calculate mean signals for each parcellation region across time points
  for (label in labels_vector) {
    region_data <- xii_mat[combined_labels == label, ]
    regional_means[which(labels_vector == label), ] <- colMeans(region_data, na.rm = TRUE)
  }
  
  return(regional_means)
}

# Example of usage:
# brain_graph_result <- fMRIBrain_Mapping(dtseries_data, statistical_data, 'Schaefer_100', 100)
