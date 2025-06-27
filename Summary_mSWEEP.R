# load data

data4 <- read.csv("/Users/maritbiesheuvel/Documents/Mbovis/data4.csv", sep = ",", stringsAsFactors = FALSE)

output_data_new <- list()

# Loop over each unique iteration number in gridfile
for (i in unique(gridfile$niter)) {
  
  # Extract PSV names for the current iteration from the gridfile
  current_iteration_psvs <- strsplit(gridfile$PSV[gridfile$niter == i], "\\|")[[1]]
  
  # Add "PSV_" in front of each PSV number
  current_iteration_psvs <- paste0("PSV_", current_iteration_psvs)
  
  # Ensure the PSV names are in character format
  current_iteration_psvs <- as.character(current_iteration_psvs)
  
  # Define the column name for the current iteration
  psv_colname <- paste0("PSV_rel_abundance_", i)
  
  if (psv_colname %in% colnames(data4)) {
    
    # Extract relevant columns
    psv_data <- data4[, c("PSV_name", psv_colname)]
    
    # Filter PSVs with abundance >= 1%
    psv_data_filtered <- psv_data[psv_data[[psv_colname]] >= 0.01, ]

    # Check if the PSVs from the gridfile are in the filtered data
    filtered_matches <- sum(current_iteration_psvs %in% psv_data_filtered$PSV_name)
    
    # Count the number of PSVs we looked for
    total_psvs <- length(current_iteration_psvs)
    
    # Count the number of PSVs that were filtered (i.e., PSV abundance ≥ 1%)
    num_filtered_psvs <- nrow(psv_data_filtered)
    
    # Calculate the percentage of PSVs found in the filtered data
    percentage_filtered <- (filtered_matches / total_psvs) * 100
    
    # Create the output dataframe
    output_dataframe <- data.frame(
      niter = i,
      PSV = paste(current_iteration_psvs, collapse = "|"),
      Filtered_matches = filtered_matches,
      Total_PSV = total_psvs,
      Percentage_in_filtered = percentage_filtered,
      Num_filtered_PSVs = num_filtered_psvs  # New column for filtered PSV count
    )
    
    # Store the result in the list
    output_data_new[[as.character(i)]] <- output_dataframe
  }
}

# Combine all data frames in the list into a single data frame
final_output <- do.call(rbind, output_data_new)



final_output$truepos <- final_output$Filtered_matches
final_output$falseneg <- final_output$Total_PSV - final_output$Filtered_matches 
final_output$falsepos <- final_output$Num_filtered_PSVs - final_output$Filtered_matches
final_output$trueneg <- (351 - final_output$Total_PSV) - final_output$falsepos  

# Sensitivity 
final_output$Se <- final_output$truepos / (final_output$truepos + final_output$falseneg) * 100
# Specificity 
final_output$Sp <- final_output$trueneg / (final_output$trueneg + final_output$falsepos) * 100
# PPV
final_output$PPV <- final_output$truepos / (final_output$truepos + final_output$falsepos) * 100
# NPV
final_output$NPV <- final_output$trueneg / (final_output$trueneg + final_output$falseneg) * 100
# FDR
final_output$FDR <- final_output$falsepos / (final_output$truepos + final_output$falsepos) * 100


# Summary - Overall 
params <- c("Se", "Sp", "PPV", "NPV", "FDR")
for (param in params) {
  summary_stats <- summary(final_output[[param]])
  std_dev <- sd(final_output[[param]], na.rm = TRUE)
  
  print(paste("Summary for", param, ":"))
  print(summary_stats)
  print(paste("Standard deviation for", param, ":"))
  print(std_dev)
}

# Summary - by enrichment 
gridfile <- gridfile %>%
  select(niter, enrichment, pcv)

final_output <- final_output %>%
  right_join(gridfile, by = "niter")

final_output <- final_output %>%
  select(niter, enrichment, pcv, everything())


#save as STATA file for later analysis
#install.packages("haven")
library(haven)
write_dta(final_output, "/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/06. Data/06. Results/msweep.dta")


library(dplyr)
library(purrr)

metrics <- c("Se", "Sp", "PPV", "NPV", "FDR")

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
results <- map(metrics, ~ calculate_summary_stats(final_output, .x))

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
results <- map(metrics, ~ calculate_summary_stats(final_output, .x))

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
results_enrichment <- map(metrics, ~ calculate_summary_stats(final_output, .x, "enrichment"))
results_pcv <- map(metrics, ~ calculate_summary_stats(final_output, .x, "pcv"))

# Combine results into one data frame
results_enrichment <- bind_rows(results_enrichment, .id = "metric") %>%
  mutate(grouping_var = "enrichment")
results_pcv <- bind_rows(results_pcv, .id = "metric") %>%
  mutate(grouping_var = "pcv")

# Combine enrichment and pcv results
combined_results <- bind_rows(results_enrichment, results_pcv)



# Plot graph for ENRICHMENT
metrics_group1 <- c(1, 2, 4)  
metrics_group2 <- c(3)  
metrics_group3 <- c(5)  

metric_labels <- c("1" = "Se", "2" = "Sp", "3" = "PPV", "4" = "NPV", "5" = "FDR")
enrichment_levels <- c(0.3, 0.5, 0.7, 0.9)

# Se, Sp and NPV
filtered_group1 <- combined_results %>%
  filter(metric %in% metrics_group1, grouping_var == "enrichment") %>%
  mutate(
    metric = factor(metric, levels = metrics_group1, labels = metric_labels[as.character(metrics_group1)]),
    enrichment = factor(enrichment, levels = enrichment_levels) 
  )
# PPV
filtered_group2 <- combined_results %>%
  filter(metric %in% metrics_group2, grouping_var == "enrichment") %>%
  mutate(
    metric = factor(metric, levels = metrics_group2, labels = metric_labels[as.character(metrics_group2)]),
    enrichment = factor(enrichment, levels = enrichment_levels)  
  )
# FDR
filtered_group3 <- combined_results %>%
  filter(metric %in% metrics_group3, grouping_var == "enrichment") %>%
  mutate(
    metric = factor(metric, levels = metrics_group3, labels = metric_labels[as.character(metrics_group3)]),
    enrichment = factor(enrichment, levels = enrichment_levels)  
  )

combined_filtered <- bind_rows(
  filtered_group1 %>% mutate(group = "Group 1"),
  filtered_group2 %>% mutate(group = "Group 2"),
  filtered_group3 %>% mutate(group = "Group 3")
)

# Create the faceted plot
plot <- ggplot(combined_filtered, aes(x = enrichment, y = mean_value, color = metric, group = metric)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = metric), alpha = 0.2) +
  labs(
    x = "Enrichment proportion",
    y = "Mean (%) with 95% CI",
    title = "mSWEEP and Themisto",
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

print(plot)

# Save the plot
ggsave("/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/Results/mSWEEP_enrichment.png", plot = plot, width = 12, height = 8, dpi = 300)
ggsave("/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/Results/mSWEEP_enrichment.pdf", plot = plot, width = 12, height = 8, dpi = 300)



# Plot graph for PSV's
metrics_group1 <- c(1, 2, 4)  
metrics_group2 <- c(3)  
metrics_group3 <- c(5)  

metric_labels <- c("1" = "Se", "2" = "Sp", "3" = "PPV", "4" = "NPV", "5" = "FDR")
psv_levels <- c(1, 3, 6, 9)

# Se, Sp and NPV
filtered_group1 <- combined_results %>%
  filter(metric %in% metrics_group1, grouping_var == "pcv") %>%
  mutate(
    metric = factor(metric, levels = metrics_group1, labels = metric_labels[as.character(metrics_group1)]),
    psv = factor(pcv, levels = psv_levels) 
  )
# PPV
filtered_group2 <- combined_results %>%
  filter(metric %in% metrics_group2, grouping_var == "pcv") %>%
  mutate(
    metric = factor(metric, levels = metrics_group2, labels = metric_labels[as.character(metrics_group2)]),
    psv = factor(pcv, levels = psv_levels)  
  )
# FDR
filtered_group3 <- combined_results %>%
  filter(metric %in% metrics_group3, grouping_var == "pcv") %>%
  mutate(
    metric = factor(metric, levels = metrics_group3, labels = metric_labels[as.character(metrics_group3)]),
    psv = factor(pcv, levels = psv_levels)  
  )

combined_filtered <- bind_rows(
  filtered_group1 %>% mutate(group = "Group 1"),
  filtered_group2 %>% mutate(group = "Group 2"),
  filtered_group3 %>% mutate(group = "Group 3")
)

# Create the faceted plot
plot <- ggplot(combined_filtered, aes(x = psv, y = mean_value, color = metric, group = metric)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = metric), alpha = 0.2) +
  labs(
    x = "Number of PSV inserted",
    y = "Mean (%) with 95% CI",
    title = "mSWEEP and Themisto",
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

print(plot)

# Save the plot
ggsave("/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/Results/mSWEEP_PSV.png", plot = plot, width = 12, height = 8, dpi = 300)
ggsave("/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projectss/09. M. bovis sequencing/Results/mSWEEP_PSV.pdf", plot = plot, width = 12, height = 8, dpi = 300)

