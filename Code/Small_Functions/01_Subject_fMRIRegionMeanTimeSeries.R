# Define a function to compute the mean fMRI time series for specified brain regions.
# The function calculates the mean across vertices for each time point.

# Inputs:
# - dtseries_data: the raw .dtseries.nii MRI data to be processed.
# - brain_parcellation: the label of the parcellation scheme, e.g., 'Schaefer_100', 'Schaefer_400', etc.
# - region_count: the total number of cortical brain regions in the parcellation 
#   (this does not include subcortex regions, which are assumed to be 19).

# Output:
# - A matrix with brain regions as rows and time points as columns (region by time).

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
# result_matrix <- fMRI_xii_pmean(dtseries_data, 'Schaefer_100', 100)
