# Define a function to process fMRI data and calculate mean time series signals
# for each brain region based on a given parcellation.

# Inputs:
# - xii: fMRI data in .dtseries.nii format.
# - parcellation: The parcellation scheme to use (e.g., 'Schaefer_100', 'Schaefer_400').
# - num: Number of regions in the parcellation (excluding subcortex regions).

# Output:
# - A matrix of mean fMRI signals with regions as rows and time points as columns (v by t).

fMRI_xii_pmean <- function(xii, parcellation = "Schaefer_100", num = 100) {
  # Load the specified parcellation data
  parc <- load_parc(parcellation)
  
  # Transform parcellation data into a vector, indicating the brain region
  # each fMRI voxel is associated with (left/right cortex)
  parc_vec <- c(as.matrix(parc))
  
  # Adjust the fMRI data to align with parcellation, setting non-applicable
  # vertices (medial wall) to NA
  xii_adj <- adjust_fMRI_data(xii) 
  
  # Create a matrix representation of the adjusted fMRI data
  xii_mat <- as.matrix(xii_adj)
  
  # Retrieve and process subcortex labels
  subcortex_labels <- xii_adj$meta$subcort$labels
  labels_all <- c(as.character(parc_vec), as.character(subcortex_labels))
  
  # Generate labels for both cortex and subcortex regions
  brain_labels <- c(as.character(1:num), levels(subcortex_labels)[3:21])
  
  # Initialize a matrix to store mean fMRI signal values
  xii_pmean <- matrix(nrow = num + 19, ncol = ncol(xii_mat))
  
  # Calculate the mean signal for each parcellation region and time point
  for (label in brain_labels) {
    region_data <- xii_mat[labels_all == label, ]
    xii_pmean[which(brain_labels == label), ] <- colMeans(region_data, na.rm = TRUE)
  }
  
  return(xii_pmean)
}

# Example of usage:
# Assuming 'xii_data' is your loaded fMRI data and 'Schaefer_100' is the parcellation
# mean_signals <- fMRI_xii_pmean(xii_data, 'Schaefer_100', 100)
