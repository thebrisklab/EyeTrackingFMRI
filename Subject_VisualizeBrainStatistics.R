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

fMRIBrain_Mapping <- function(dtseries_data, mapping_data, brain_parcellation = "Schaefer_100", region_count = 100) {
  # Load the parcellation data
  parc <- load_parc(brain_parcellation)
  
  # Convert the parcellation to a vector to identify the brain region of each voxel
  parc_vec <- c(as.matrix(parc))
  
  # Adjust subcortex labels to match cortical labeling
  subcortex_labels_adjusted <- as.numeric(dtseries_data$meta$subcort$labels) - 2
  subcortex_labels_adjusted <- region_count + subcortex_labels_adjusted
  
  # Create a combined brain vector of cortical and adjusted subcortex labels
  brain_vec <- c(parc_vec, subcortex_labels_adjusted)
  
  # Adjust the fMRI data for non-cortical vertices to align with parcellation
  adjusted_data <- move_from_mwall(dtseries_data, NA)
  
  # Initialize fMRI data with empty values for subsequent mapping
  prepared_data <- select_xifti(adjusted_data * 0, 1)
  
  # If the data represents p-values, perform -log10 transformation for visualization
  # Uncomment the following line if necessary
  # transformed_data <- -log10(mapping_data)
  
  # Map the statistical data onto the brain vector
  brain_graph <- newdata_xifti(prepared_data, c(NA, mapping_data)[brain_vec + 1])
  
  # Return the brain graph object for visualization
  return(brain_graph)
}

# Example of usage:
# brain_graph_result <- fMRIBrain_Mapping(dtseries_data, statistical_data, 'Schaefer_100', 100)
