#### Summarizing Kraken2 Output ####


#### add PSVs to gridfile ####
library(dplyr)
library(stringr)
library(readr)

# Read the data files
gridfile <- read.csv("/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/06. Data/02. Data (Fasta and gbff)/Data_simulation/grid_parameters_new.csv", sep = ",", stringsAsFactors = FALSE)
psv_data <- read_delim("/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/06. Data/02. Data (Fasta and gbff)/Final_output_TreeCluster_labels_avg_clade_t_0.0001_sorted.txt")

# Update the PSV_name in psv_data to only contain numbers
psv_data <- psv_data %>%
  rename(Filenames = SequenceName) %>%
  mutate(Filenames = str_replace(Filenames, "$", ".fna_post_extract.fna")) %>%
  rename(PSV_name = ProposedLabel) %>%
  select(-ClusterNumber, -Singleton_ClusterName) %>%
  mutate(PSV_number = str_extract(PSV_name, "\\d+"))  # Extract numbers from PSV_name

# Function to get PSV numbers for a given filename or set of filenames
get_psv_numbers <- function(filenames, psv_data) {
  filenames_split <- unlist(strsplit(filenames, "\\|"))
  psv_numbers <- sapply(filenames_split, function(f) {
    psv_row <- psv_data %>% filter(Filenames == f)
    if (nrow(psv_row) == 0) return(NA) else return(psv_row$PSV_number)
  })
  
  return(paste(psv_numbers, collapse = "|"))
}

# Apply the function to each row of the gridfile
gridfile <- gridfile %>%
  rowwise() %>%
  mutate(PSV = get_psv_numbers(Filenames, psv_data)) %>%
  ungroup()




#### Percentage of reads classified as correct PSV ####
# Define the folder containing all your output files
directory_path <- "/Users/maritbiesheuvel/Documents/output_kraken"

# List all files in the directory
file_list <- list.files(directory_path, pattern = "sampled_reads_.*\\.kreport2$", full.names = TRUE)

# Initialize an empty list to hold results
output_list <- list()

# Loop over each file in the file list
for (file in file_list) {
  # Read the data from the current file
  niter_data <- read.delim(file, header = FALSE)
  
  # Rename columns
  names(niter_data) <- c("perc_frag_covered_clade", "N_frag_covered_clade", "N_frag_assigned_clade", "rank", "NCBI_taxID", "PSVs")
  
  # Remove first 8 rows
  niter_data <- niter_data[-c(1:8), ]
  
  # Remove unnecessary columns
  niter_data <- subset(niter_data, select = -c(N_frag_covered_clade, N_frag_assigned_clade, rank, NCBI_taxID))
  
  # Extract the iteration number from the file name
  iteration <- as.integer(sub("sampled_reads_(\\d+)_mb_classification.kreport2", "\\1", basename(file)))
  
  # Add iteration column
  niter_data$niter <- iteration
  
  # Extract numbers using regular expression
  niter_data$PSV <- as.numeric(gsub("[^0-9]", "", niter_data$PSVs))
  niter_data$PSV[is.na(niter_data$PSV)] <- niter_data$PSVs[is.na(niter_data$PSV)]
  
  # Extract PSV for the current iteration from Grid file
  current_iteration_psvs <- strsplit(as.character(gridfile$PSV[gridfile$niter == iteration]), "\\|")[[1]]
  current_iteration_psvs <- as.numeric(current_iteration_psvs)
  
  # Initialize a data frame to hold the results for the current iteration
  iteration_results <- data.frame(niter = iteration)
  
  # Loop over each PSV and extract the corresponding row from niter_data
  for (i in seq_along(current_iteration_psvs)) {
    psv <- current_iteration_psvs[i]
    matching_row <- niter_data[niter_data$PSV == psv, ]
    
    # If matching row is found, add it to the results; else add a row with perc_frag_covered_clade = 0
    iteration_results[[paste0("perc_frag_covered_clade_", i)]] <- if (nrow(matching_row) > 0) {
      matching_row$perc_frag_covered_clade
    } else {
      0
    }
  }
  
  # Add the results for the current iteration to the output list
  output_list[[length(output_list) + 1]] <- iteration_results
}

# Sort the output list by niter
output_list <- output_list[order(sapply(output_list, function(x) x$niter))]

# Determine the maximum number of PSV columns across all data frames
max_columns <- max(sapply(output_list, function(df) ncol(df) - 1)) # Subtract 1 for the 'niter' column

# Function to standardize columns to the maximum number of columns
standardize_columns <- function(df, max_columns) {
  # Create column names for PSV_1 to PSV_max_columns
  required_columns <- paste0("perc_frag_covered_clade_", 1:max_columns)
  missing_columns <- setdiff(required_columns, names(df))
  df[missing_columns] <- NA # Fill missing columns with NA
  df <- df[order(names(df))] # Order columns
  return(df)
}

# Standardize all data frames to the maximum number of columns
output_list_standardized <- lapply(output_list, standardize_columns, max_columns = max_columns)

# Combine standardized data frames
output_data <- do.call(rbind, output_list_standardized)

# Rename columns to PSV_1 to PSV_9 for clarity
colnames(output_data) <- gsub("perc_frag_covered_clade_", "PSV_", colnames(output_data))

# Order output_data by niter from 1 to 16000
output_data <- output_data[order(output_data$niter), ]

# Add a column total_perc that sums up all perc_frag_covered_clade columns
output_data$total_perc <- rowSums(output_data[, grepl("PSV_", names(output_data))], na.rm = TRUE)

summary_perc <- summary(output_data$total_perc)
print(summary_perc)
sd_perc <- sd(output_data$total_perc)
print(sd_perc)

# add simulation parameters to output file 
# Ensure that niter is present in both dataframes and is of the same type
gridfile <- gridfile %>%
  select(niter, enrichment, pcv)

output_data <- output_data %>%
  right_join(gridfile, by = "niter")

output_data <- output_data %>%
  select(niter, enrichment, pcv, everything())


# Perc by enrichment 
summary_stats <- output_data %>%
  group_by(enrichment) %>%
  summarize(
    mean_value = mean(total_perc, na.rm = TRUE),
    median_value = median(total_perc, na.rm = TRUE),
    sd_value = sd(total_perc, na.rm = TRUE),
    q1_value = quantile(total_perc, 0.25, na.rm = TRUE),  # 1st quartile (25th percentile)
    q3_value = quantile(total_perc, 0.75, na.rm = TRUE),  # 3rd quartile (75th percentile)
    n = n()  # Count of observations in each group
  )
print(summary_stats)

# Perc by PSV 
summary_stats <- output_data %>%
  group_by(pcv) %>%
  summarize(
    mean_value = mean(total_perc, na.rm = TRUE),
    median_value = median(total_perc, na.rm = TRUE),
    sd_value = sd(total_perc, na.rm = TRUE),
    q1_value = quantile(total_perc, 0.25, na.rm = TRUE),  # 1st quartile (25th percentile)
    q3_value = quantile(total_perc, 0.75, na.rm = TRUE),  # 3rd quartile (75th percentile)
    n = n()  # Count of observations in each group
  )
print(summary_stats)



#PPV by enrichment and PSV
summary_stats <- output_data_topten_new %>%
  group_by(pcv, enrichment) %>%
  summarize(
    mean_value = mean(PPV, na.rm = TRUE),
    median_value = median(PPV, na.rm = TRUE),
    sd_value = sd(PPV, na.rm = TRUE),
    min_value = min(PPV, na.rm = TRUE),
    max_value = max(PPV, na.rm = TRUE),
    n = n()  # Count of observations in each group
  )
print(summary_stats)




# Define breaks and labels
breaks <- c(-Inf, 0, 0.05, 0.10, 0.15, 0.20, 0.40, 0.60, 0.80, 1.00, 1.50, 2.00, 2.50, 3.00, 4.00, 5.00, 10.00, 12.00, Inf)

labels <- c("0", ">0 - 0.05", "0.05 - 0.10", "0.10 - 0.15", "0.15 - 0.20", "0.20 - 0.40", "0.40 - 0.60", 
            "0.60 - 0.80", "0.80 - 1.00", "1.00 - 1.50", "1.50 - 2.00", "2.00 - 2.50", "2.50 - 3.00", 
            "3.00 - 4.00", "4.00 - 5.00", "5.00 - 10.00", "10.00 - 12.00", ">12.00")

# Apply cut() function
output_data$group <- cut(output_data$total_perc, breaks = breaks, labels = labels, include.lowest = TRUE, right = FALSE)

# Explicitly set zeros to the "0" category
output_data$group[output_data$total_perc == 0] <- "0"

# Create the table
table(output_data$group)





#### PSV meaningful present in kraken2 output file (>1%) ####

# Define the folder containing all your output files
directory_path <- "/Users/maritbiesheuvel/Documents/output_kraken"

# List all files in the directory
file_list <- list.files(directory_path, pattern = "sampled_reads_.*\\.kreport2$", full.names = TRUE)

# Initialize an empty list to store output data frames
output_list_new <- list()

# Loop over each file in the file list
for (file in file_list) {
  # Read the data from the current file
  niter_data <- read.delim(file, header = FALSE)
  
  # Rename columns
  names(niter_data)[names(niter_data) == "V1"] <- "perc_frag_covered_clade"
  names(niter_data)[names(niter_data) == "V2"] <- "N_frag_covered_clade"
  names(niter_data)[names(niter_data) == "V3"] <- "N_frag_assigned_clade"
  names(niter_data)[names(niter_data) == "V4"] <- "rank"
  names(niter_data)[names(niter_data) == "V5"] <- "NCBI_taxID"
  names(niter_data)[names(niter_data) == "V6"] <- "PSVs"
  
  # Remove first 8 rows
  niter_data <- niter_data[-c(1:8), ]
  
  # Remove unnecessary columns
  niter_data <- subset(niter_data, select = -c(N_frag_covered_clade, N_frag_assigned_clade, rank, NCBI_taxID))
  
  # Extract numbers using regular expression
  niter_data$PSV <- as.numeric(gsub("[^0-9]", "", niter_data$PSVs))
  niter_data$PSV[is.na(niter_data$PSV)] <- niter_data$PSVs[is.na(niter_data$PSV)]
  
  # Remove PSVs column
  niter_data <- subset(niter_data, select = -c(PSVs))
  
  # Filter rows where perc_frag_covered_clade >= 1
  niter_data_filtered <- niter_data[niter_data$perc_frag_covered_clade >= 1, ]
  
  # Count the number of PSVs after filtering
  total_psvs_filtered <- nrow(niter_data_filtered)
  
  # Extract iteration number from the file name
  iteration <- as.integer(sub("sampled_reads_(\\d+)_mb_classification.kreport2", "\\1", basename(file)))
  
  # Extract PSV for the current iteration from the grid file
  current_iteration_psvs <- strsplit(as.character(gridfile$PSV[gridfile$niter == iteration]), "\\|")[[1]]
  current_iteration_psvs <- as.numeric(current_iteration_psvs)
  
  # Check if any of the current_iteration_psvs are present in the filtered PSV column
  matching_rows <- which(niter_data_filtered$PSV %in% current_iteration_psvs)
  
  # Count the number of PSVs we looked for
  total_psvs <- length(current_iteration_psvs)
  
  # Count how many PSVs were found in the filtered data
  matches_found <- length(matching_rows)
  
  # Calculate the percentage of PSVs found in the filtered data
  percentage_matched <- (matches_found / total_psvs) * 100
  
  # Create output dataframe
  output_dataframe <- data.frame(
    niter = iteration,
    PSV = paste(current_iteration_psvs, collapse = "|"),
    Matches_found = matches_found,
    Total_PSVs = total_psvs,
    Percentage = percentage_matched,
    Total_PSVs_filtered = total_psvs_filtered  # Add column for total PSVs after filtering
  )
  
  # Add the output dataframe to the output list
  output_list_new[[iteration]] <- output_dataframe
}

# Combine all data frames in output_list_new into a single data frame
output_data_new <- do.call(rbind, output_list_new)

output_data_new$truepos <- output_data_new$Matches_found
output_data_new$falseneg <- output_data_new$Total_PSVs - output_data_new$Matches_found 
output_data_new$falsepos <- output_data_new$Total_PSVs_filtered - output_data_new$Matches_found
output_data_new$trueneg <- (351 - output_data_new$Total_PSVs) - output_data_new$falsepos  

# Sensitivity 
output_data_new$Se <- output_data_new$truepos / (output_data_new$truepos + output_data_new$falseneg) * 100
# Specificity 
output_data_new$Sp <- output_data_new$trueneg / (output_data_new$trueneg + output_data_new$falsepos) * 100
# PPV
output_data_new$PPV <- output_data_new$truepos / (output_data_new$truepos + output_data_new$falsepos) * 100
# NPV
output_data_new$NPV <- output_data_new$trueneg / (output_data_new$trueneg + output_data_new$falseneg) * 100
# FDR
output_data_new$FDR <- output_data_new$falsepos / (output_data_new$truepos + output_data_new$falsepos) * 100


# Summary - Overall 
params <- c("Se", "Sp", "PPV", "NPV", "FDR")
for (param in params) {
  summary_stats <- summary(output_data_new[[param]])
  std_dev <- sd(output_data_new[[param]], na.rm = TRUE)

  print(paste("Summary for", param, ":"))
  print(summary_stats)
  print(paste("Standard deviation for", param, ":"))
  print(std_dev)
}

# Summary - by enrichment 
gridfile <- gridfile %>%
  select(niter, enrichment, pcv)

output_data_new <- output_data_new %>%
  right_join(gridfile, by = "niter")

output_data_new <- output_data_new %>%
  select(niter, enrichment, pcv, everything())

#save as STATA file for later analysis
install.packages("haven")
library(haven)
write_dta(output_data_new, "/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/06. Data/06. Results/kraken2.dta")


library(dplyr)
library(purrr)

metrics <- c("PPV", "Se", "Sp", "NPV", "FDR")

# enrichment 
calculate_summary_stats <- function(data, metric) {
  data %>%
    group_by(enrichment) %>%
    summarize(
      mean_value = mean(.data[[metric]], na.rm = TRUE),
      median_value = median(.data[[metric]], na.rm = TRUE),
      sd_value = sd(.data[[metric]], na.rm = TRUE),
      Q1_value = quantile(.data[[metric]], probs = 0.25, na.rm = TRUE),
      Q3_value = quantile(.data[[metric]], probs = 0.75, na.rm = TRUE),
      n = n()  # Count of observations in each group
    )
}

# Loop over the metrics and calculate summary statistics for each
results <- map(metrics, ~ calculate_summary_stats(output_data_new, .x))

# Print results for each metric
names(results) <- metrics
results

# PSVs
# Define a function to calculate summary statistics for a given metric
calculate_summary_stats <- function(data, metric) {
  data %>%
    group_by(pcv) %>%
    summarize(
      mean_value = mean(.data[[metric]], na.rm = TRUE),
      median_value = median(.data[[metric]], na.rm = TRUE),
      sd_value = sd(.data[[metric]], na.rm = TRUE),
      Q1_value = quantile(.data[[metric]], probs = 0.25, na.rm = TRUE),
      Q3_value = quantile(.data[[metric]], probs = 0.75, na.rm = TRUE),
      n = n()  # Count of observations in each group
    )
}

# Loop over the metrics and calculate summary statistics for each
results <- map(metrics, ~ calculate_summary_stats(output_data_new, .x))

# Print results for each metric
names(results) <- metrics
results


## plot graph 
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

# Define the metrics to loop over
metrics <- c("Se", "Sp", "PPV", "NPV", "FDR")

# Define a function to calculate summary statistics for a given metric and grouping variable
# Define a function to calculate summary statistics for a given metric and grouping variable
calculate_summary_stats <- function(data, metric, group_by_var) {
  data %>%
    group_by(across(all_of(group_by_var))) %>%
    summarize(
      mean_value = mean(.data[[metric]], na.rm = TRUE),
      sd_value = sd(.data[[metric]], na.rm = TRUE),
      n = n(),  # Count of observations in each group
      .groups = 'drop'
    ) %>%
    mutate(
      se = sd_value / sqrt(n),  # Standard error
      ci_lower = mean_value - 1.96 * se,  # 95% CI lower bound
      ci_upper = mean_value + 1.96 * se   # 95% CI upper bound
    )
}

# Calculate summary stats for each metric and grouping variable
results_enrichment <- map(metrics, ~ calculate_summary_stats(output_data_new, .x, "enrichment"))
results_pcv <- map(metrics, ~ calculate_summary_stats(output_data_new, .x, "pcv"))

# Combine results into one data frame
results_enrichment <- bind_rows(results_enrichment, .id = "metric") %>%
  mutate(grouping_var = "enrichment")
results_pcv <- bind_rows(results_pcv, .id = "metric") %>%
  mutate(grouping_var = "pcv")

# Combine enrichment and pcv results
combined_results <- bind_rows(results_enrichment, results_pcv)



# Plot graph for ENRICHMENT
metrics_group1 <- c(1, 3)  
metrics_group2 <- c(2, 4, 5)  

metric_labels <- c("1" = "Se", "2" = "Sp", "3" = "PPV", "4" = "NPV", "5" = "FDR")
enrichment_levels <- c(0.3, 0.5, 0.7, 0.9)

# Se and PPV
filtered_group1 <- combined_results %>%
  filter(metric %in% metrics_group1, grouping_var == "enrichment") %>%
  mutate(
    metric = factor(metric, levels = metrics_group1, labels = metric_labels[as.character(metrics_group1)]),
    enrichment = factor(enrichment, levels = enrichment_levels) 
  )
# Sp, NPV and FDR
filtered_group2 <- combined_results %>%
  filter(metric %in% metrics_group2, grouping_var == "enrichment") %>%
  mutate(
    metric = factor(metric, levels = metrics_group2, labels = metric_labels[as.character(metrics_group2)]),
    enrichment = factor(enrichment, levels = enrichment_levels)  
  )

combined_filtered <- bind_rows(
  filtered_group1 %>% mutate(group = "Group 1"),
  filtered_group2 %>% mutate(group = "Group 2")
)

# Create the faceted plot
plot <- ggplot(combined_filtered, aes(x = enrichment, y = mean_value, color = metric, group = metric)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = metric), alpha = 0.2) +
  labs(
    x = "Enrichment proportion",
    y = "Mean (%) with 95% CI",
    title = "Kraken2",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),  
    strip.text = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ group, scales = "free_y")  # Separate plots with free y-axis scales

# Save the plot
ggsave("/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/Results/kraken2_enrichment.png", plot = plot, width = 12, height = 8, dpi = 300)
ggsave("/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/Results/kraken2_enrichment.pdf", plot = plot, width = 12, height = 8, dpi = 300)



# Plot graph for PSV's
metrics_group1 <- c(1, 3)  
metrics_group2 <- c(2, 4, 5)  

metric_labels <- c("1" = "Se", "2" = "Sp", "3" = "PPV", "4" = "NPV", "5" = "FDR")
psv_levels <- c(1, 3, 6, 9)

# Se and PPV
filtered_group1 <- combined_results %>%
  filter(metric %in% metrics_group1, grouping_var == "pcv") %>%
  mutate(
    metric = factor(metric, levels = metrics_group1, labels = metric_labels[as.character(metrics_group1)]),
    psv = factor(pcv, levels = psv_levels) 
  )
# Sp, NPV and FDR
filtered_group2 <- combined_results %>%
  filter(metric %in% metrics_group2, grouping_var == "pcv") %>%
  mutate(
    metric = factor(metric, levels = metrics_group2, labels = metric_labels[as.character(metrics_group2)]),
    psv = factor(pcv, levels = psv_levels)  
  )

combined_filtered <- bind_rows(
  filtered_group1 %>% mutate(group = "Group 1"),
  filtered_group2 %>% mutate(group = "Group 2")
)

# Create the faceted plot
plot <- ggplot(combined_filtered, aes(x = psv, y = mean_value, color = metric, group = metric)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = metric), alpha = 0.2) +
  labs(
    x = "Number of PSV inserted",
    y = "Mean (%) with 95% CI",
    title = "Kraken2",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),  
    strip.text = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ group, scales = "free_y")  # Separate plots with free y-axis scales

# Save the plot
ggsave("/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/Results/kraken2_PSV.png", plot = plot, width = 12, height = 8, dpi = 300)
ggsave("/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/Results/kraken2_PSV.pdf", plot = plot, width = 12, height = 8, dpi = 300)










                ######## OLD ##########

#### PSV in top 10 of Kraken2 output file ####

# Define the folder containing all your output files
directory_path <- "/Users/maritbiesheuvel/Documents/output_kraken"

# List all files in the directory
file_list <- list.files(directory_path, pattern = "sampled_reads_.*\\.kreport2$", full.names = TRUE)

# Initialize an empty list to store output data frames
output_list_topten_new <- list()

# Loop over each file in the file list
for (file in file_list) {
  # Read the data from the current file
  niter_data <- read.delim(file, header = FALSE)
  
  # Rename columns
  names(niter_data)[names(niter_data) == "V1"] <- "perc_frag_covered_clade"
  names(niter_data)[names(niter_data) == "V2"] <- "N_frag_covered_clade"
  names(niter_data)[names(niter_data) == "V3"] <- "N_frag_assigned_clade"
  names(niter_data)[names(niter_data) == "V4"] <- "rank"
  names(niter_data)[names(niter_data) == "V5"] <- "NCBI_taxID"
  names(niter_data)[names(niter_data) == "V6"] <- "PSVs"
  
  # Remove first 8 rows
  niter_data <- niter_data[-c(1:8), ]
  
  # Remove unnecessary columns
  niter_data <- subset(niter_data, select = -c(N_frag_covered_clade, N_frag_assigned_clade, rank, NCBI_taxID))
  
  # Extract numbers using regular expression
  niter_data$PSV <- as.numeric(gsub("[^0-9]", "", niter_data$PSVs))
  niter_data$PSV[is.na(niter_data$PSV)] <- niter_data$PSVs[is.na(niter_data$PSV)]
  
  # Remove PSVs column
  niter_data <- subset(niter_data, select = -c(PSVs))
  
  # Extract iteration number from the file name
  iteration <- as.integer(sub("sampled_reads_(\\d+)_mb_classification.kreport2", "\\1", basename(file)))
  
  # Extract PSV for the current iteration from Grid file
  current_iteration_psvs <- strsplit(as.character(gridfile$PSV[gridfile$niter == iteration]), "\\|")[[1]]
  current_iteration_psvs <- as.numeric(current_iteration_psvs)
  
  # Check if any of the current_iteration_psvs are present in the PSV column of niter_data
  matching_rows <- which(niter_data$PSV %in% current_iteration_psvs)
  
  # Check if any of the matching rows are in the top ten
  top_ten_matches <- sum(niter_data$PSV[1:10] %in% current_iteration_psvs)
  
  # Count the number of PSVs we looked for
  total_psvs <- length(current_iteration_psvs)
  
  # Calculate the percentage of PSVs found in the top ten
  percentage_top_ten <- (top_ten_matches / total_psvs) * 100
  
  # Create output dataframe
  output_dataframe <- data.frame(niter = iteration, PSV = paste(current_iteration_psvs, collapse = "|"), Top_ten = top_ten_matches, totPSV = total_psvs, Percentage = percentage_top_ten)
  
  # Add the output dataframe to the output list
  output_list_topten_new[[iteration]] <- output_dataframe
}

# Combine all data frames in output_list_topten into a single data frame
output_data_topten_new <- do.call(rbind, output_list_topten_new)


output_data_topten_new$falsepos <- output_data_topten_new$totPSV - output_data_topten_new$Top_ten
output_data_topten_new$truepos <- output_data_topten_new$Top_ten 
output_data_topten_new$trueneg <- 0

#### SENSITIVITY 
output_data_topten_new$Se <- (output_data_topten_new$Top_ten / output_data_topten_new$totPSV)*100

# If you want to summarize the sensitivity across all iterations (e.g., mean sensitivity)
mean_Se <- mean(output_data_topten_new$Se, na.rm = TRUE)

# If you want to get other statistics (e.g., median, standard deviation)
summary_Se <- summary(output_data_topten_new$Se)
std_dev_Se <- sd(output_data_topten_new$Se, na.rm = TRUE)

print(summary_Se)
print(std_dev_Se)

#### PPV
output_data_topten_new$PPV <- (output_data_topten_new$truepos / (output_data_topten_new$truepos + output_data_topten_new$falsepos))*100

# If you want to summarize the sensitivity across all iterations (e.g., mean sensitivity)
mean_PPV <- mean(output_data_topten_new$PPV, na.rm = TRUE)

# If you want to get other statistics (e.g., median, standard deviation)
summary_PPV <- summary(output_data_topten_new$PPV)
std_dev_PPV <- sd(output_data_topten_new$PPV, na.rm = TRUE)

print(summary_PPV)
print(std_dev_PPV)

#### FDR
output_data_topten_new$FDR <- (output_data_topten_new$falsepos / (output_data_topten_new$falsepos + output_data_topten_new$truepos))*100

# If you want to summarize the sensitivity across all iterations (e.g., mean sensitivity)
mean_FDR <- mean(output_data_topten_new$FDR, na.rm = TRUE)

# If you want to get other statistics (e.g., median, standard deviation)
summary_FDR <- summary(output_data_topten_new$FDR)
std_dev_FDR <- sd(output_data_topten_new$FDR, na.rm = TRUE)

print(summary_FDR)
print(std_dev_FDR)


# add simulation parameters to output file 
# Ensure that niter is present in both dataframes and is of the same type
gridfile <- gridfile %>%
  select(niter, enrichment, pcv)

output_data_topten_new <- output_data_topten_new %>%
  right_join(gridfile, by = "niter")

output_data_topten_new <- output_data_topten_new %>%
  select(niter, enrichment, pcv, everything())

# PPV by enrichment 
summary_stats <- output_data_topten_new %>%
  group_by(enrichment) %>%
  summarize(
    mean_value = mean(PPV, na.rm = TRUE),
    median_value = median(PPV, na.rm = TRUE),
    sd_value = sd(PPV, na.rm = TRUE),
    min_value = min(PPV, na.rm = TRUE),
    max_value = max(PPV, na.rm = TRUE),
    n = n()  # Count of observations in each group
  )
print(summary_stats)

# PPV by PSV 
summary_stats <- output_data_topten_new %>%
  group_by(pcv) %>%
  summarize(
    mean_value = mean(PPV, na.rm = TRUE),
    median_value = median(PPV, na.rm = TRUE),
    sd_value = sd(PPV, na.rm = TRUE),
    min_value = min(PPV, na.rm = TRUE),
    max_value = max(PPV, na.rm = TRUE),
    n = n()  # Count of observations in each group
  )
print(summary_stats)

#PPV by enrichment and PSV
summary_stats <- output_data_topten_new %>%
  group_by(pcv, enrichment) %>%
  summarize(
    mean_value = mean(PPV, na.rm = TRUE),
    median_value = median(PPV, na.rm = TRUE),
    sd_value = sd(PPV, na.rm = TRUE),
    min_value = min(PPV, na.rm = TRUE),
    max_value = max(PPV, na.rm = TRUE),
    n = n()  # Count of observations in each group
  )
print(summary_stats)

