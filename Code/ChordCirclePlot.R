##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
# Need to be revised. Updated on 2024, May 15
# Script to generate a circle plot (chord diagram) visualizing relationships between groups.

# Load necessary libraries
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)

# Load your statistical matrix
z_statistic_matrix <- zstat_mat_1887002

# Set the group names for columns and rows
colnames(z_statistic_matrix) <- parc_cluster$Group
rownames(z_statistic_matrix) <- parc_cluster$Group

# Convert the matrix to a data frame for processing
df <- as.data.frame(as.table(z_statistic_matrix))

# Optional filtering based on the significance of frequency
# df_filtered <- subset(df, abs(df$Freq) > 1.96)
df_filtered <- df

# Ensure that the factor levels are consistent for proper plotting
df_filtered$Var1 <- factor(df_filtered$Var1, levels = unique(df_filtered$Var2))
df_filtered$Var2 <- factor(df_filtered$Var2, levels = unique(df_filtered$Var2))

# Initialize the circos plot parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# Define a color palette
my_colors <- viridis(17, alpha = 1, begin = 0, end = 1, option = "D")
my_colors <- my_colors[sample(1:17)]

# Extract and convert the unique RGB colors to hexadecimal
unique_colors <- unique(parc$meta$cifti$labels[[1]][-1,][, c("Red", "Green", "Blue", "Alpha")])
hex_colors <- rgb(unique_colors)

# Assign colors based on the sign of the frequency value
df_filtered$FreqFactor <- ifelse(df_filtered$Freq < 0, "blue", "red")

# Base plot using chord diagram
chordDiagramFromDataFrame(
  df = df_filtered[, c(2, 1, 3:5)],
  col = df_filtered$FreqFactor,
  grid.col = hex_colors,
  transparency = 0.3,
  annotationTrack = "grid",
  annotationTrackHeight = c(0.05, 0.1),
  direction.type = "diffHeight",
  directional = 1,
  diffHeight = -0.04
)

# Add labels and annotations
circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add sector labels
    circos.text(
      x = mean(xlim),
      y = 2,
      labels = sector.index,
      facing = "bending",
      cex = 1,
      font = 2
    )
  }
)

# Add graduation on axis
circos.axis(
  h = "top",
  major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)),
  minor.ticks = 1,
  major.tick.length = 0.5,
  labels.niceFacing = FALSE)
 }
)
