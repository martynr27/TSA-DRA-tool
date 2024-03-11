# Temporary Storage Area Drainage Rate Analysis (TSA-DRA) tool

# Martyn Roberts (University of Aberdeen & James Hutton Institute)
# Tue Jul 11 2023 ------------------------------

# The Temporary Storage Area Drainage Rate Analysis (TSA-DRA) tool is a novel
# data-based mechanistic method that only requires rainfall and level data to 
# describe individual TSA drainage rates. The tool provides a systematic 
# approach for characterising the functioning of a wide range of TSA types and 
# sizes. The tool can also explore time-variable TSA functioning and provides 
# useful soil infiltration estimates for modelling. 

# Method:
# This script extracts recession periods based on user input and then creates a
# master recession curve (MRC) to provide a general description for TSA drainage 
# over a longer time period. A segmented linear model is then fitted to the MRC 
# to describe the TSA drainage rate. User input is based on site knowledge and
# the desired MRC analysis. 

# Data requirements:
# TSA data | 15 min | Headings: date_time, temp, depth 
# Precipitation | hourly | Headings: date_time, mm_hr
# API | daily | Headings: date_time, mm_daily, api_daily


# Packages ----------------------------------------------------------------

library(tidyverse)
library(zoo)
library(data.table)
library(strucchange)
library(segmented)


# Functions ---------------------------------------------------------------

source("R/TSA-DRA_tool_functions.R")


# ### USER INPUT ### ------------------------------------------------------

### Recession criteria ###
MA <- 5                     # TSA moving average
min.dry <- 2                # Min dry period (Precip = hourly, so 2 = 2 hours ) 
maxP.neg <- 0.2             # max precipitation negligibility within dry period
min.Rlength <- 1 * 4        # Min recession length ( _ hours * 4(15min))
min.temp <- 0               # Min temperature (avoid freezing errors)
max.TSA <- 0.4              # Max TSA water level

# Time-variable
season <- "Spring"
yr <- 2015


### TSA dimensions ### 
# Outlet pipe
pipe_diameter <- 0.22
pipe_area <- pi * (pipe_diameter ^ 2 / 4)
pipe_height_m <- 0.25
pipe_base_m <- pipe_height_m - (pipe_diameter / 2)
pipe_top_m <- pipe_height_m + (pipe_diameter / 2)


# Function to identify recessions 
# Subsets TSA DF that meets the following recession criteria 
RA_criteria.fn <- function(x) {
  subset(x[x$diff_Q < 0 &
             x$cum_precip <= maxP.neg &
             x$MA_Q <= max.TSA &
             x$temp > min.temp 
           
           # # Time-variable
           # x$season == season
           # x$year == yr
           # 
           # # Additional criteria e.g.
           # x$MA_Q <= pipe_base_m           # only soil infiltration
           # x$stage < stream_spill_height   # Extract when no stream spill (Offline pond)
           
           ,])
}



# Recession criteria DF
R_criteria_DF <-
  data.frame(MA,
             maxP.neg,
             min.dry,
             min.Rlength,
             min.temp,
             max.TSA,
             season,
             yr)


# # csv file
# write.table(
#   R_criteria_DF,
#   file = 'outputs/RAcriteria.csv',
#   col.names = TRUE,
#   row.names = FALSE,
#   sep = ","
# )


# Import data -------------------------------------------------------------

## TSA data (can use either level or volume) ##
 
TSA_DF <-
  read.table(
    file = "data/processed_data/Bund/datasets/TarBund_15min_MASTER.csv",   ### USER INPUT ###
    header = TRUE,
    sep = ",",
    dec = "."
  ) %>% 
  transform(date_time = as.POSIXct(date_time,
                                   format = "%Y-%m-%d %H:%M",
                                   tz = "UTC")) %>% 
  transform(N_dt = as.numeric(date_time)) %>%                
  dplyr::select(date_time, N_dt, everything()) %>%                   
  arrange(N_dt)

str(TSA_DF)
summary(TSA_DF)


## precipitation data (hourly) ##
rain_DF <- read.table(
  file = "data/processed_data/precipitation/hourly_precip_MASTER.csv",   ### USER INPUT ###
  header = TRUE,
  sep = ",",
  dec = "."
) %>% 
  transform(date_time = as.POSIXct(date_time,
                                   format = "%Y-%m-%d %H:%M",
                                   tz = "UTC")) %>% 
  transform(N_dt = as.numeric(date_time)) %>%               
  dplyr::select(date_time, N_dt, everything()) %>%                           
  arrange(N_dt)

str(rain_DF)
summary(rain_DF)
table(rain_DF$mm_hr)



# TSA moving average ------------------------------------------------------

# Split TSA DF when the time interval is not continuous
time_inc <- 900     # 15 min interval data (15min * 60 sec = 900)
TSA_list <- time_split.fn(TSA_DF)      

# Moving average (MA) and difference in depth or volume (diff_Q)
TSA_list <- TSA_list %>%
  lapply(transform,             
         MA_Q = rollmean(     
           depth,               ### depth/level or volume ###
           k = MA,           ### rolling mean length ###
           align = "center",    ### alignment ###
           fill = NA            
         )) %>%
  lapply(transform,                              
         diff_Q = c(NA, diff(MA_Q)))         

# Check NA values
TSA_list %>% lapply(function(x)
  sum(is.na(x)))

# Remove NA values for analysis
TSA_list <- TSA_list %>% lapply(na.omit)     

# Create TSA DF
TSA_DF <- bind_rows(TSA_list) %>% arrange(N_dt)

str(TSA_DF)
summary(TSA_DF)

# # Depth and MA_Q plot
# TSA_DF %>%
#   ggplot() +
#   geom_line(aes(x = date_time, y = depth)) +
#   geom_line(aes(x = date_time, y = MA_Q), color = "red") +
#   labs(x = "Date",
#        y = "Depth")




# Precipitation -----------------------------------------------------------

### Calculate cumulative precipitation ###

# Split precipitation DF by time breaks
time_inc <- 60*60     # 60 mins  
precip_list <- time_split.fn(rain_DF)

# Check NA values
precip_list %>% lapply(function(x)
  sum(is.na(x)))

# Cumulative precipitation 
precip_list <- precip_list %>%
  lapply(transform,
         cum_precip = rollapply(
           mm_hr,
           width = min.dry,     ### cumulative sum for the length of min.dry ###
           FUN = sum,
           align = "right",     ### alignment ###
           fill = NA
         ))

# Create precipitation DF
precip_DF <- bind_rows(precip_list) %>% arrange(N_dt)

# Check cumulative precipitation values 
table(precip_DF$cum_precip)


### Create 15-minute precipitation data ###

# Identify minimum and maximum date_time 
min_ts <- min(precip_DF$date_time)     # min date
max_ts <- max(precip_DF$date_time)     # max date

# Create a new dataframe with 15-minute intervals
new_df <- data.frame(date_time = seq(min_ts, max_ts, by = "15 min"))

# Merge the new dataframe with the original dataframe
merged_df <- merge(precip_DF, new_df, by = "date_time", all = TRUE)

# Sort the merged dataframe by date_time
precip_DF <- merged_df[order(merged_df$date_time), ]

precip_DF$N_dt <- as.numeric(precip_DF$date_time)


# Next observation carried backward - fill NA values with hourly data 
precip_DF[, c("mm_hr", "cum_precip")] <-
  na.locf(precip_DF[, c("mm_hr", "cum_precip")], fromLast = TRUE)

head(precip_DF)
summary(precip_DF)


# Create recession analysis DF --------------------------------------------

# Merge TSA and precipitation DFs
RA_DF <-
  left_join(TSA_DF, precip_DF, by = "date_time") %>%          # merge by TSA_DF date_time
  transform(N_dt = as.numeric(date_time)) %>%                 # create numeric date_time column
  arrange(N_dt) %>%  
  metseasons.fn() %>%
  mutate(year = format(date_time, format = "%Y")) %>% 
  dplyr::select(date_time, N_dt, year, season, everything(), -N_dt.x, -N_dt.y)            

# Check RA DF
head(RA_DF)
summary(RA_DF)


# Identify recession periods  ---------------------------------------------

# Subset RA DF based on recession criteria
RA_filtered <- RA_criteria.fn(RA_DF)

# Split RA_DF by time breaks
time_inc <- 900     # 15 min interval data (15min * 60 sec)
RA_list <- time_split.fn(RA_filtered)

# Check time step breaks
tail(RA_list[[1]], 1)
head(RA_list[[2]], 1)

# Filter RA_list (>= min.Rlength)
R_periods_list <-
  Filter(function(x)
    nrow(x) >= min.Rlength, RA_list)

# Label each dataframe with recession id
R_id <- seq(1:length(R_periods_list))

R_periods_list <-
  mapply(cbind,
         R_periods_list,
         "R_id" = R_id,
         SIMPLIFY = FALSE) %>% 
  lapply(transform,
         rel_time = seq(1:length(date_time)))

# Create R_period DF
R_period_DF <- bind_rows(R_periods_list)

# Check R_period DF
head(R_period_DF)

# Add number of recessions column 
R_criteria_DF$n.recessions <- max(R_period_DF$R_id)


# MRC matching method -----------------------------------------------------

# duration_match()
MRC_long <- duration_match(R_period_DF, min.Rlength)


# MRC descriptive data  ---------------------------------------------------

# Descriptive data 
recession_data <- R_period_DF %>%
  distinct(R_id, .keep_all = TRUE) %>% 
  dplyr::select(-rel_time, -N_dt)

MRC_data <- MRC_long %>% 
  distinct(R_id2, .keep_all = TRUE)

MRC_data <- 
  left_join(MRC_data, recession_data, "R_id")


# # Save MRC descriptive data
# # csv file
# write.table(
#   MRC_data,
#   file = 'outputs/MRC_descriptive_data.csv',
#   col.names = TRUE,
#   row.names = FALSE,
#   sep = ","
# )


MRC_data_filtered <- MRC_data %>% 
  dplyr::select(R_id2, year, season)

MRC_long <- MRC_long %>%
  left_join(MRC_data_filtered, by = "R_id2")


# MRC plot ----------------------------------------------------------------

# Convert time (t) to hours 
MRC_long$t <- MRC_long$t*15 / 60

# Extract x values for plot 
min_x <- min(MRC_long$t, na.rm = TRUE)
max_x <- max(MRC_long$t, na.rm = TRUE)

# Order seasons 
MRC_long$season <-
  factor(MRC_long$season,
         levels = c("Autumn", "Winter", "Spring", "Summer"))        

# Seasonal colours
cbp2 <- c("#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "red", "#000000")


### MRC plot ###
plot <- MRC_long %>%
  ggplot(aes(x = t,
             y = Q,
             group = R_id2, 
             colour = season)) +
  geom_line() +
  labs(x = "Time (hours)",
       y = "TSA depth (m)") +
  geom_hline(
    yintercept = c(pipe_base_m),
    linetype = 'dashed',
    color = 'red',
    linewidth = 1
  ) +
  scale_colour_manual(values = cbp2) +   ### seasonal colours ###
  annotate("text", x = 2, y = (pipe_base_m), label = "Outlet pipe", size = 8/.pt) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = "top") 

plot

max(MRC_long$R_id)
max(MRC_long$t)

# # Save MRC plot
# ggsave(
#   filename = "MRC.tiff",
#   plot = plot,
#   device = "tiff",
#   path = "outputs/Figures",
#   width = 160,
#   height = 100,
#   units = "mm",
#   dpi = 600,
#   limitsize = TRUE
# )



# Extract MRC points ------------------------------------------------------

# Duration match method
MRC_fit <- MRC_long %>%
  group_by(t) %>%
  slice_min(order_by = Q) %>%    
  dplyr::select(t, Q) %>%   
  distinct(t, .keep_all = TRUE)

## Filter MRC_fit if sparse data at higher levels ##
# MRC_fit <- MRC_fit[MRC_fit$Q <= 0.3, ]

plot <- MRC_fit %>%
  ggplot(aes(x = t,
             y = Q)) +
  geom_point() +
  geom_hline(yintercept = pipe_base_m)

plot


MRC_fit$MRC <- "2015"   ### USER INPUT - MRC name ###

# ### Save MRC points ###
# # csv file
# write.table(
#   MRC_fit,
#   file = 'outputs/MRC_points.csv',
#   col.names = TRUE,
#   row.names = FALSE,
#   sep = ","
# )

# Append MRC_fit to saved MRC points df
saved_df <-
  read.csv('outputs/MRC_points.csv')

combined_df <- rbind(saved_df, MRC_fit)

# ### Save appended MRC points ###
# # csv file
# write.table(
#   combined_df,
#   file = 'outputs/MRC_points.csv',
#   col.names = TRUE,
#   row.names = FALSE,
#   sep = ","
# )


# Fit linear model to MRC -------------------------------------------------

### Identify optimal breakpoints ###
bpts <- breakpoints(Q ~ t, data = MRC_fit)
MRC_breaks <- bpts$breakpoints

# Q breaks
MRC_Q_breaks <- MRC_fit[MRC_breaks, "Q"][[1]]
MRC_Q_breaks

MRC_Q_breaks <- sort(MRC_Q_breaks, decreasing = TRUE)  
MRC_Q_breaks

# t breaks
closest_indices <- sapply(MRC_Q_breaks, function(x) which.min(abs(MRC_fit$Q - x)))
MRC_t_breaks <- MRC_fit$t[closest_indices]
MRC_t_breaks

### Segmented regression ###
# Fit the initial model
m1 <- lm(Q ~ t, data = MRC_fit)

# Fit the segmented models
s1 <- segmented(m1)  # One breakpoint
s2 <- segmented(m1, psi = c(MRC_t_breaks))  # Optimal breakpoints

# Create data frames based on predictions
plot_data <- data.frame(t = MRC_fit$t, Q = MRC_fit$Q)
s1_data <- data.frame(t = MRC_fit$t, Q = predict(s1))
s2_data <- data.frame(t = MRC_fit$t, Q = predict(s2))


# Plot the data and segmented regression lines
ggplot() +
  geom_point(data = plot_data, aes(x = t, y = Q)) +
  # geom_line(data = s1_data, aes(x = t, y = Q), color = "blue") +
  geom_line(data = s2_data, aes(x = t, y = Q), color = "red", size = 1.5)


### USER INPUT - MRC name ###
s2_data$MRC <- "2015"

# ### Save segmented lm prediction ###
# # csv file
# write.table(
#   s2_data,
#   file = 'outputs/MRC_lm.csv',
#   col.names = TRUE,
#   row.names = FALSE,
#   sep = ","
# )


# Append dfs
saved_df <-
  read.csv('outputs/MRC_lm.csv')

combined_df <- rbind(saved_df, s2_data)

# ### Save appended df ###
# # csv file
# write.table(
#   combined_df,
#   file = 'outputs/MRC_lm.csv',
#   col.names = TRUE,
#   row.names = FALSE,
#   sep = ","
# )


### Extract estimated breakpoints ###
Est_t_bkps <- as.numeric(summary(s2)$psi[, "Est."])

# Create a new data frame with the specific time values and estimated break points
new_data <- data.frame(t = Est_t_bkps)

# Predict the values for Q using the specific times and break points
Est_Q_bkps <- as.numeric(predict(s2, newdata = new_data))
Est_Q_bkps

# Create data frame with estimated breakpoints
Est_bkps <- data.frame(Est_t_bkps, Est_Q_bkps)


### Linear model information ###
## Depth breaks ##
depth_breaks <- c(max(MRC_fit$Q, na.rm = TRUE), Est_Q_bkps, 0)

# Create depth ranges
Depth_range <- character(length(depth_breaks) - 1)

for (i in 1:(length(depth_breaks) - 1)) {
  if (depth_breaks[i+1] == 0) {
    range <- sprintf("%.2f to 0", depth_breaks[i])
  } else {
    range <- sprintf("%.2f to %.2f", depth_breaks[i], depth_breaks[i+1])
  }
  Depth_range[i] <- range
}

Depth_range


## Linear model equation ##
## Extract slope and intercept values ## 
Slope <- as.numeric(slope(s2)$t[, "Est."])
Intercept <- as.numeric(intercept(s2)$t[, "Est."])

## lm (no breakpoints) ## 
# Slope <- as.numeric(coef(m1)["t"])
# Intercept <- as.numeric(coef(m1)["(Intercept)"])

# Create data frame with linear model information 
MRC_lm <- data.frame(Depth_range, Slope, Intercept)

## Linear model equation ##
# Create an empty vector for storing the linear model equations
model_equations <- character()

# Iterate over each row of the MRC_lm data frame
for (i in 1:nrow(MRC_lm)) {
  slope <- MRC_lm$Slope[i]
  intercept <- MRC_lm$Intercept[i]
  
  # Create the linear model equation
  equation <- paste("y =", slope, "* x +", intercept)
  
  # Append the equation to the model_equations vector
  model_equations <- c(model_equations, equation)
}

# Print the linear model equations
print(model_equations)

# Store equations in data frame
MRC_lm$model_equation <- model_equations
MRC_lm


# Categorise MRC_fit by breakpoints
breaks <- c(-Inf, Est_Q_bkps, Inf)
MRC_Q_bkps <- cut(MRC_fit$Q, breaks = breaks, right = FALSE)


# Subset data based on MRC_Q_bkps categories
subset_data <- function(data, categories) {
  subsets <- list()
  
  for (category in categories) {
    subset <- data[MRC_Q_bkps == category, ]
    subsets[[category]] <- subset
  }
  
  return(subsets)
}

categories <- unique(MRC_Q_bkps)
data_subsets <- subset_data(MRC_fit, categories)

## R2 value ##
# Create an empty vector to store the rounded R-squared values
R_squared_values <- numeric()

# Iterate over each subset of data and corresponding model equation
for (i in 1:length(data_subsets)) {
  # Extract the subset of data and model equation
  subset_data <- data_subsets[[i]]
  equation <- model_equations[i]
  
  # Extract the Q and t values from the subset data
  Q <- subset_data$Q
  t <- subset_data$t
  
  # Evaluate the predicted Q values using the model equation
  predicted_Q <- eval(parse(text = gsub("x", "t", equation)))
  
  # Calculate the R-squared value
  R_squared <- 1 - sum((Q - predicted_Q)^2) / sum((Q - mean(Q))^2)
  
  # Round the R-squared value to 2 decimal places
  R_squared_rounded <- round(R_squared, 2)
  
  # Append the rounded R-squared value to the R_squared_values vector
  R_squared_values <- c(R_squared_values, R_squared_rounded)
}

# Add the rounded R-squared values to the MRC_lm data frame
MRC_lm$R_squared <- R_squared_values

# # ## lm (no breakpoints) ## 
# R_squared <- summary(m1)$r.squared
# R_squared_rounded <- round(R_squared, 2)
# MRC_lm$R_squared <- R_squared_rounded

## Add name of MRC ##
MRC_lm$MRC <- "2015"   ### USER INPUT ###

# Reorder the columns in the MRC_lm data frame
MRC_lm <- MRC_lm[, c("MRC", "Depth_range", "model_equation", "Slope", "Intercept", "R_squared")]
MRC_lm

# ### Save lm prediction information ###
# # csv file
# write.table(
#   MRC_lm,
#   file = 'outputs/MRC_lm_info.csv',
#   col.names = TRUE,
#   row.names = FALSE,
#   sep = ","
# )


# Append dfs
saved_df <-
  read.csv('outputs/MRC_lm_info.csv')

combined_df <- rbind(saved_df, MRC_lm)

# ### Save appended df ###
# # csv file
# write.table(
#   combined_df,
#   file = 'outputs/MRC_lm_info.csv',
#   col.names = TRUE,
#   row.names = FALSE,
#   sep = ","
# )