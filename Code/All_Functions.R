# Define a function to compute the mean fMRI time series for specified brain regions.
# The function calculates the mean across vertices for each time point.

# Inputs:
# - dtseries_data: the raw .dtseries.nii MRI data to be processed.
# - brain_parcellation: the label of the parcellation scheme, e.g., 'Schaefer_100', 'Schaefer_400', etc.
# - region_count: the total number of cortical brain regions in the parcellation 
#   (this does not include subcortex regions, which are assumed to be 19).

# Output:
# - A matrix with brain regions as rows and time points as columns (region by time).

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

#########################################################################################################
#########################################################################################################


