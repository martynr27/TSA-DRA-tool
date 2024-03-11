# Temporary Storage Area Drainage Rate Analysis (TSA-DRA) tool functions

# Martyn Roberts (University of Aberdeen & James Hutton Institute)
# Tue Jul 11 2023 ------------------------------

# This script contains functions that have been created for the TSA-DRA tool. 


# Time split function -----------------------------------------------------

# Identify time breaks in data and splits into a list
time_split.fn <- function(x) {
  breaks <- c(1, diff(x$N_dt)) > time_inc           # Identify breaks in DF when time step is greater than time_inc
  split_id <- which(breaks =="TRUE")                # identify row numbers where condition is met
  split(x, f = cumsum(1:nrow(x) %in% split_id))     # create a list of dataframes based on split_id
}



# Season function ---------------------------------------------------------

# Creates a meteorological seasons column in DF
metseasons.fn <- function(x) {
  metseasons <- c(
    "01" = "Winter", "02" = "Winter",
    "03" = "Spring", "04" = "Spring", "05" = "Spring",
    "06" = "Summer", "07" = "Summer", "08" = "Summer",
    "09" = "Autumn", "10" = "Autumn", "11" = "Autumn",
    "12" = "Winter"
  )
  x$season <- metseasons[format(x$date_time, "%m")]
  return(x)
}



# MRC matching method -----------------------------------------------------

### Duration match ###
duration_match <- function(R_period_DF, min.Rlength) {
  # Create a vector to store the values for R_id2
  R_id2_values <- character()
  
  # Iterate through each unique value in R_id
  for (unique_id in unique(R_period_DF$R_id)) {
    # Get the subset of rows for the current R_id
    subset_rows <- R_period_DF[R_period_DF$R_id == unique_id, ]
    
    # Determine the number of times to repeat the value of R_id
    repeat_times <- nrow(subset_rows)
    
    # Generate the suffix to append to R_id based on the number of repetitions
    suffix <- rep(LETTERS[seq_len(ceiling(repeat_times / min.Rlength))], each = min.Rlength)[seq_len(repeat_times)]
    
    # Append the modified R_id values to the R_id2_values vector
    R_id2_values <- c(R_id2_values, paste0(unique_id, suffix))
  }
  
  # Add the R_id2 column to the dataframe
  R_period_DF$R_id2 <- R_id2_values
  
  # Identify each recession start date
  R_startDates <- R_period_DF %>% 
    group_by(R_id2) %>% 
    filter(date_time == min(date_time)) %>% 
    dplyr::select(R_id2, date_time)
  
  # Order recessions based on each min value (high to low)
  R_minvalues <- R_period_DF %>%      
    group_by(R_id2) %>%                # group by recession id
    slice(which.min(MA_Q)) %>%      # extract min depth/volume values for each recession
    dplyr::select(R_id2, MA_Q)              # select R_id and depth/volume data
  
  x <- R_minvalues[with(R_minvalues, order(-MA_Q)), ]     # order by MA-vol min value (high to low)
  R_id_order <- x$R_id2     # R_id order - recession min values high to low 
  
  # Convert dataframe to wide format
  MRC_DF <- R_period_DF %>% 
    arrange(match(R_id2, R_id_order)) %>%        # match R_id rows with R_id_order 
    transform(t = 1:nrow(R_period_DF)) %>%     # create t (time) column the length of the df
    dplyr::select(t, R_id2, MA_Q)                    # select columns
  
  R_order_chr <- as.character(R_id_order)     # convert R_id_order to character so can rearrange columns
  
  MRC_DF_wide <- dcast(setDT(MRC_DF), t ~ R_id2, value.var = "MA_Q") %>%     # convert df to wide format 
    dplyr::select(all_of(R_order_chr))       # rearrange column order
  
  colnames(MRC_DF_wide)[1:ncol(MRC_DF_wide)] <-
    paste("R", colnames(MRC_DF_wide)[1:ncol(MRC_DF_wide)], sep = "_")     # add "R_" to start of numbers so syntax is correct
  
  # Create empty DF to put at the top of MRC_DF_wide - able to shift values up without losing data
  NA_DF <- MRC_DF_wide      
  NA_DF[] <- NA
  MRC_DF_wide <- rbind(NA_DF, MRC_DF_wide)
  
  # Match lowest point with the corresponding value in the next column
  matching.strip.fn <- function(.x, .y) {
    .y <- .y[!is.na(.y)]
    pos <- which.min(.x) - which.min(abs(min(.x, na.rm = T) - .y))
    c(rep(NA, pos), .y, rep(NA, length(.x) - pos - length(.y)))
  }
  
  MRC_DF_wide <- MRC_DF_wide %>% 
    accumulate(matching.strip.fn) %>%
    set_names(names(MRC_DF_wide)) %>%
    as.data.frame()
  
  # Remove rows with all NAs
  MRC_DF_wide <- MRC_DF_wide[rowSums(is.na(MRC_DF_wide)) != ncol(MRC_DF_wide),]
  
  # Create time(t) column
  MRC_DF_wide <- MRC_DF_wide %>%
    transform(t = 1:nrow(MRC_DF_wide)) %>%
    dplyr::select(t, everything())
  
  # Convert to long format
  MRC_long <- MRC_DF_wide %>% gather(R_id2, Q,-c(t))
  
  # Remove "R_" from recession id and convert to number 
  MRC_long$R_id2 <- gsub("\\R_","", MRC_long$R_id2)
  
  # Create R_id column to use as match for DF merge
  MRC_long$R_id <- as.numeric(gsub("[A-Za-z]", "", MRC_long$R_id2)) 

  return(MRC_long)
}