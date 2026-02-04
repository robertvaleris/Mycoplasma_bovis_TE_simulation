# Enrichments             (4): 0.3, 0.5, 0.7 and 0.9
# N of PCVs               (4): 1, 3, 6 and 9 
# of iterations per combi    : 1,000
# Sequencing depth           : Randomly from poisson distribution 


start_enrich <- 0.3
iter = matrix(1:1000, nrow = 1000)
griddf <- expand.grid(enrichment = seq(from = start_enrich, by = 0.2, l = 4),
                      pcv = c(1, 3, 6, 9),
                      iter = 1000)
griddf


griddf2 <- as.data.frame(lapply(griddf, rep, griddf$iter))
griddf2$iter <- NULL

start_niter <- 1
niter = seq(from = start_niter, by = 1, l = 16000)

griddf2$niter <- niter

## add sequencing depth - random draw from Poisson distribution with mean = 20,000,000
set.seed(123)
seqdepth <- rpois(n = 16000, lambda = 20000000)
hist(seqdepth)

griddf2$seqdepth <- seqdepth


#final file
griddf2 <- griddf2[, c(3,1,2,4)]
griddf2$targetNreads <- round(griddf2$seqdepth * griddf2$enrichment)


####### add the names of genomes that need to be sampled to that list ####### 

# extract list of unique PSV names 
library(readxl)
PSV <- read_excel("Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projects /09. M. bovis sequencing/06. Data/02. Data (Fasta and gbff)/genomelist_memberships_final.xlsx")

# Set seed
set.seed(123)

# List of PSV names
psv_to_sample <- as.list(PSV$PSV)
psv_to_sample <- unique(psv_to_sample)

# Create a list to store sampled genomes for each sampled PSV
sampled_data <- list()

# Iterate through each row of griddf2
for (row_index in 1:nrow(griddf2)) {
  
  # Get the number of iterations for the current row
  niter <- griddf2$niter[row_index]
  
  # Extract the predetermined number of PCVs for the current iteration
  current_iteration_pcv <- griddf2$pcv[row_index]
  
  # Sample the predetermined number of PCVs for the current iteration
  sampled_pcv_list <- sample(psv_to_sample, size = current_iteration_pcv, replace = FALSE)
  
  # Initialize a list to store sampled filenames for each PCV
  sampled_filenames_list <- list()
  
  # Iterate over each sampled PCV
  for (selected_pcv in sampled_pcv_list) {
    
    # Get all genome filenames associated with the selected PCV
    selected_genomes <- PSV$filename[PSV$PSV == selected_pcv]
    
    # Sample 1 genome filename
    sampled_genome <- sample(selected_genomes, size = 1)
    
    # Append the sampled genome filename to the list
    sampled_filenames_list[[length(sampled_filenames_list) + 1]] <- sampled_genome
  }
  
  # Flatten the list to a vector
  sampled_filenames <- unlist(sampled_filenames_list)
  
  # Modify the filename for each sampled genome
  curated_filenames <- paste(sampled_filenames, "_post_extract.fna", sep = "")
  
  # Store the sampled data in a list for each iteration
  sampled_data[[length(sampled_data) + 1]] <- c(niter = niter, ID = list(sampled_pcv_list), Filenames = curated_filenames)
  
  # Add the sampled filenames to the corresponding row in griddf2
  griddf2$Filenames[row_index] <- list(curated_filenames)
}




## put it in a vector list 
# Create an empty list to store the vectors
list_of_vectors <- list()

# Iterate through each row of the dataframe
for (row_index in 1:nrow(griddf2)) {
  # Extract parameters for the current row
  enrichment <- griddf2$enrichment[row_index]
  pcv <- griddf2$pcv[row_index]
  seqdepth <- griddf2$seqdepth[row_index]
  targetNreads <- griddf2$targetNreads[row_index]
  niter <- griddf2$niter[row_index]
  filenames <- griddf2$Filenames[row_index]
  
  # Create a named vector
  current_vector <- c(enrichment = enrichment, pcv = pcv, seqdepth = seqdepth, targetNreads = targetNreads, niter = niter, filenames = filenames)
  
  # Add the vector to the list
  list_of_vectors[[row_index]] <- current_vector
}

# Display the list of vectors
#list_of_vectors

# Convert the list of vectors to a data frame
df <- do.call(rbind, list_of_vectors)

# Write the data frame to a CSV file
write.csv(df, "/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projects /09. M. bovis sequencing/06. Data/01. Scripts/04. Simulation/grid_parameters", row.names = FALSE)

#get all the names of the genomes of the list_of_vectors in one list
#create an empty list 
all_filenames <- list()
#Iterate through each vector in list_of_vectors
for (i in seq_along(list_of_vectors)) {
  #Extract filenames from the current vector
  filenames <- unlist(list_of_vectors[[i]]["filenames"])
  
  #Add filenames to list
  all_filenames <- c(all_filenames, filenames)
}

# Remove duplicates and create a unique list of filenames
unique_filenames <- unique(all_filenames)

install.packages("openxlsx")
library(openxlsx)
excel_filepath <- "/Users/maritbiesheuvel/Library/CloudStorage/OneDrive-UniversityofCalgary/02. PhD UCalgary/04. PhD-projects /09. M. bovis sequencing/06. Data/02. Data (Fasta and gbff)/Unique_genomes_for_simulation/unique_genomelist.xlsx"
wb <- createWorkbook()
addWorksheet(wb, "Filenames")

for (i in seq_along(unique_filenames)) {
  writeData(wb, sheet = "Filenames", x = unique_filenames[i], startRow = i, startCol = 1)
}
saveWorkbook(wb, file = excel_filepath)




