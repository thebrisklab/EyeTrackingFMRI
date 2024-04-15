##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
################################## Code used to calculate the percentage of ET events in each video clip. (Median, Q1, Q3) ###############################################
##########################################################################################################################################################################
##########################################################################################################################################################################

# Function to extract eyeblink and eyefixation events for each video segment across two task sessions.
# This function processes the Eye Tracking data to compute the count and duration percentage of events
# within specified video intervals.

# task session 1
# video 1. 0021SAND 3 - 130
# video 2. 0023DOLL 145 - 252
# video 3. 0024DOLL 267 - 333

# task session 2
# video 4 0022SAND 3 - 241
# video 5 0025DOLL 256 - 390

# Input:
# - path_ET_ses1: Path to the Eye Tracking data file for task session 1.
# - path_ET_ses2: Path to the Eye Tracking data file for task session 2.
# - Tasksession1: Boolean to indicate if session 1 should be processed.
# - Tasksession2: Boolean to indicate if session 2 should be processed.
# - ID: Participant ID for annotation.

# Output:
# - A dataframe summarizing eyeblink and eyefixation events for each video.

Function_summary_ET <- function(path_ET_ses1, path_ET_ses2, Tasksession1 = FALSE, Tasksession2 = FALSE, ID = ID) {
  # Prepare an empty data frame to collect results
  results <- data.frame()
  
  # Helper function to process session data
  process_session <- function(path, session_number, video_intervals) {
    ETasc <- read.asc(fname = path, samples = TRUE, events = TRUE)
    process_type <- function(type, label) {
      ETtypein <- ETasc[[type]] %>% mutate(ts = (stime - min(ETasc$raw$time)) / 1e3,
                                           te = (etime - min(ETasc$raw$time)) / 1e3,
                                           duration = dur / 1e3)
      lapply(video_intervals, function(interval) {
        subset(ETtypein, ts >= interval$start & ts <= interval$end) %>%
          summarise(event = n(), 
                    percent = sum(duration) * 100 / (interval$end - interval$start), 
                    video = interval$name, 
                    Movie = session_number, 
                    ETtype = label,
                    ID = ID)
      })
    }
    # Combine blink and fixation data
    rbind(
      do.call(rbind, process_type("blinks", "Blink")),
      do.call(rbind, process_type("fix", "Fixation"))
    )
  }
  
  # Define video intervals for each session
  video_intervals_ses1 <- list(
    list(name = "0021SAND", start = 3, end = 130),
    list(name = "0023DOLL", start = 145, end = 252),
    list(name = "0024DOLL", start = 267, end = 333)
  )
  video_intervals_ses2 <- list(
    list(name = "0022SAND", start = 3, end = 241),
    list(name = "0025DOLL", start = 256, end = 390)
  )
  
  # Process each session based on input flags
  if (Tasksession1) {
    results <- rbind(results, process_session(path_ET_ses1, 1, video_intervals_ses1))
  }
  if (Tasksession2) {
    results <- rbind(results, process_session(path_ET_ses2, 2, video_intervals_ses2))
  }
  
  return(results)
}

# Example of usage:
# summary_data <- Function_summary_ET(path_ET_ses1 = "path/to/session1.et", 
#                                     path_ET_ses2 = "path/to/session2.et", 
#                                     Tasksession1 = TRUE, 
#                                     Tasksession2 = TRUE, 
#                                     ID = "Participant123")

##########################################################################################################################################################################
##########################################################################################################################################################################

#####################################################################################
# Script to extract eye tracking summary for individual participants across different videos.
# This script uses the `Function_summary_ET` function to process data and then merges the results
# with participant demographic and clinical information.

# Example calls to Function_summary_ET for a few individuals:
# In real time analysis, need to do this for everyone
ETsummary.1898902 <- Function_summary_ET(path_ET_ses2 = "path/to/189890202.asc",
                                         Tasksession2 = TRUE, ID = 1898902 )

ETsummary.1876402 <- Function_summary_ET(path_ET_ses1 = "path/to/1876401.asc",
                                         Tasksession1 = TRUE, ID = 1876402)

ETsummary.1879803 <- Function_summary_ET(path_ET_ses1 = "path/to/1879831.asc",
                                         path_ET_ses2 = "path/to/1879832.asc",
                                         Tasksession1 = TRUE, Tasksession2 = TRUE, ID = 1879803)

# Merge the eye tracking summary, and add ASD/TD label, as well as Sex
# Create a list of all objects in your environment
all_objects <- ls()

# Use a regular expression to find names that match the 'ETsummary.xxxx' pattern
pattern <- "^ETsummary\\.\\d+$"
ETsummary.name  <- grep(pattern, all_objects, value = TRUE)

# Loop through the list and rbind each data frame
ETsummary <- data.frame()
for(df_name in ETsummary.name) {
  # Assume the data frames are in the global environment
  # Use get() to retrieve them by name
  df <- get(df_name)
  # Combine the data frame with the previous ones
  ETsummary <- rbind(ETsummary, df)
}

#### read ASD, TD file for merging together.
asd.td.info <- read.csv("path/to/ASD_TD_label.tsv", header = T)
# remove the "-" in the first column
asd.td.info$src_subject_id <- as.numeric(gsub("-", "", asd.td.info$src_subject_id))
# merge with combined_data_frame
ETsummary.asd <- left_join(ETsummary, asd.td.info, by = c("ID" = "src_subject_id"))

#### read the file with SEX
MRIcontrol <- read.csv("path/to/eachsubject_gender.csv", header = T)
# remove the "-"
MRIcontrol$Individual <- as.numeric(gsub("-", "", MRIcontrol$Individual))

# merge the ET summary and ASD_TD and Gender, need to revise based on specific condition
ETsummary.asd.sex <- left_join(ETsummary.asd, MRIcontrol[!duplicated(MRIcontrol[, 1:4]), ][,1:4], by = c("ID" = "Individual"))

##################################### Table 1.
# Analyze and summarize data for real Table 1.
ET_summary_analysis <- ETsummary.asd.sex %>%
  group_by(phenotype, Sex, video, ETtype) %>%
  summarise(
    Median = median(percent, na.rm = TRUE),
    Q1 = quantile(percent, 0.25, na.rm = TRUE),
    Q3 = quantile(percent, 0.75, na.rm = TRUE)
  )

# Output the summary to a CSV file
write.csv(ET_summary_analysis, "path/to/ET_summary_analysis.csv")

# Print the result for review
print(ET_summary_analysis)

