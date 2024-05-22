# Define a function to process eye-tracking data from an .ASC file and extract
# either blink or fixation data, creating a timestamped dataset for further analysis.

# Inputs:
# - file_path: The path to the eye-tracking .ASC file.
# - blink: Logical flag to indicate if blink data should be extracted.
# - fixation: Logical flag to indicate if fixation data should be extracted.

# Output:
# - A list containing the total time of the eye-tracking task, onsets, durations,
#   and sampling rate accuracy.

ETascDataProcess <- function(file_path, blink = FALSE, fixation = FALSE) {
  # Read the raw .asc data
  ETasc <- read.asc(fname = file_path, samples = TRUE, events = TRUE)
  # Initialize an empty container for the extracted type
  ETtypein <- NULL
  
  # Extract blink or fixation data and calculate real time
  if (blink) {
    # eyeQuality to detect the eye-blink event
    eyeQuality.result <- detectBlinks(data = ETasc$raw,
                                      pupilLeft  = "ps",
                                      recHz = 500)
    # process to get the ts, te, duration
    # get the whole ET blink subset dataframe
    etblink <- subset(eyeQuality.result, eyeQuality.result$pupilLeft.blink == 1)
    # calculate each blink's end time point, and store the index
    et.end.indice <- which(diff(etblink$newtime) > 0.003)
    # we only need to keep the start and end time for each eye blink event
    etblink.se <- etblink[sort(c(1, et.end.indice, et.end.indice + 1, nrow(etblink))),]
    # add lables
    etblink.se$event <- rep(c("ts", "te"), time = nrow(etblink.se)/2)
    etblink.se$blink_seq <- rep(c(1:(nrow(etblink.se)/2)), each = 2)
    # remove rownames
    rownames(etblink.se) <- NULL
    
    # convert to wide format, no need for pupil information
    ETtypein <- etblink.se[,c("newtime","event", "blink_seq")] %>% group_by(blink_seq) %>% 
      pivot_wider(names_from = "event", values_from = "newtime") %>%
      dplyr::mutate(duration = te - ts)

  } else if (fixation) {
    ETtypein <- ETasc$fix %>% 
      dplyr::mutate(ts = (stime - min(ETasc$raw$time)) / 1e3,
                    te = (etime - min(ETasc$raw$time)) / 1e3,
                    duration = dur / 1e3)
  } else {
    stop("Please specify blink or fixation data to process.")
  }
  
  # Convert the results to a data.frame
  ETtype <- as.data.frame(ETtypein)
  # Calculate the total time, onsets, and durations for further analysis
  totaltime <- (max(ETasc$raw$time) - min(ETasc$raw$time)) / 1e3
  onsets <- ETtype$ts
  durations <- ETtype$duration
  sampling_rate <- 0.002
  
  return(list(totaltime, onsets, durations, sampling_rate))
}

# Example of usage:
# et_data <- ETascDataProcess(file_path = "path/to/data.asc", blink = TRUE, fixation = FALSE)
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

# Define a function to perform the convolution of eye-tracking time series data
# with a hemodynamic response function (HRF).

# Inputs:
# - totaltime: Total duration of the stimulus experiment in seconds.
# - onsets: A vector containing the start times for each stimulus event.
# - durations: A vector containing the duration of each stimulus event in seconds.
# - sampling_rate: The sampling rate for the experiment (e.g., 0.002 sec for 500 Hz).

# Output:
# - A list containing the convolution result and the corresponding time vector.

Convolution_function <- function(totaltime, onsets, durations, sampling_rate) {
  # Generate a stimulus boxcar function based on the inputs
  stimulus <- stimfunction(totaltime = totaltime, onsets = onsets, durations = durations, accuracy = sampling_rate)
  # Calculate the time vector for the length of the experiment
  time_vector <- seq(0, totaltime, by = sampling_rate)
  
  # Perform convolution with the stimulus and HRF
  # (Ensure 'stimulus' and 'newhrf' have the same length for later convolution)
  stimulus_padded <- c(stimulus, rep(0, length(time_vector) - length(stimulus)))
  hrf_kernel <- canonicalHRF(time_vector, verbose = FALSE) # Double-gamma Function
  hrf_padded <- c(hrf_kernel, rep(0, length(stimulus)))
  
  # Compute the frequency domain representation via FFT
  fft_stimulus <- fft(stimulus_padded)
  fft_hrf <- fft(hrf_padded)
  # Multiply element-wise in the frequency domain and perform inverse FFT
  convolved_signal <- Re(fft(fft_stimulus * fft_hrf, inverse = TRUE) / length(stimulus))
  # Trim the convolution result to the original time vector length
  convolved_signal <- convolved_signal[1:length(time_vector)]
  
  return(list(convolved_signal, time_vector))
}

# Example of usage:
# convolution_result <- Convolution_function(totaltime = 780, onsets = c(1, 100, 200), durations = c(1, 1, 1),  accuracy = 0.002)
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

# Define a function to extract time series data from convolved eye-tracking (ET) data
# and align it with fMRI data time points.

# Inputs:
# - conv_data: Convolution result from the ET data processed with HRF.
# - fmri_data: Pre-processed fMRI time series data obtained from earlier analysis.
# - tr: Repetition time or the time interval of fMRI data acquisition (in seconds).

# Output:
# - A vector of the standardized extracted ET time series aligned with fMRI data.

Extraction_ETtime <- function(conv_data, fmri_data, tr = 1.127) {
  # Calculate the indices to extract from the convolved data based on fMRI TR
  indices <- round(seq(from = 1, to = length(conv_data[[1]]), by = tr / 0.002), digits = 0)
  
  # Extract the convolved data at the specified indices
  extracted_data <- conv_data[[1]][indices]
  
  # Ensure the extracted ET time series length matches the fMRI time series
  if (ncol(fmri_data) > length(extracted_data)) {
    # Pad the end with zeros if the fMRI data is longer
    extracted_data <- c(extracted_data, rep(0, ncol(fmri_data) - length(extracted_data)))
  } else {
    # Truncate the extracted data if it is longer than the fMRI data
    extracted_data <- extracted_data[1:ncol(fmri_data)]
  }
  
  # Standardize the extracted time series
  standardized_data <- extracted_data / max(extracted_data)
  
  return(standardized_data)
}

# Example of usage:
# Assuming 'convolved_hrf_data' is your convolved ET data and 'fmri_time_series' is the fMRI time series
# standardized_et_series <- Extraction_ETtime(conv_data = convolved_hrf_data, fmri_data = fmri_time_series, tr = 1.127)
