library(seqinr)

apply_sliding_window <- function(sublist_names, sequences, output_file) {
  window_size <- 150
  
  for (sublist_name in sublist_names) {
    sequence <- sequences[[sublist_name]]
    genome_length <- nchar(sequence)
    num_iterations <- genome_length - window_size + 1
    result_list <- vector("list", length = num_iterations)
    
    # Iterate over the loop
    for (i in 1:num_iterations) {
      # Calculate the start and end indices for the sliding window
      start_index <- i
      end_index <- start_index + window_size - 1
      
      # Extract the substring for the current sliding window
      current_substr <- substr(sequence, start_index, end_index)
      
      # Store the result in the list
      result_list[[i]] <- current_substr
    }
    
    # Concatenate the content of the sublist
    result_content <- paste(sapply(seq_along(result_list), function(n) {
      paste(">",sublist_name, "_read", n, "\n", result_list[[n]], sep = "")
    }), collapse = "\n")
    
    # Append the concatenated content to the output file
    cat(result_content, "\n\n", file = output_file, append = TRUE)
    
    # Post-process the file to remove blank lines
    lines <- readLines(output_file)
    filtered_lines <- lines[lines != ""]  # Exclude blank lines
    writeLines(filtered_lines, output_file)
    
  }
}

# Directory containing genome files
genome_directory <- "/home/marit.biesheuvel/allgenomic_fna"

# Output directory for result files 
output_directory <- "/home/marit.biesheuvel/mycoplasma_reads"

# List all files in the directory
genome_files <- list.files(path = genome_directory, pattern = "\\.fna$", full.names = TRUE)

# Loop over each genome file
for (genome_file in genome_files) {
  # Read the genome sequence from the file
  genome_sequence <- read.fasta(file = genome_file, as.string = TRUE, forceDNAtolower = FALSE, whole.header = F)  
  
  # Create the output folder (if it doesn't exist)
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }
  
  # Create an output file for the current genome
  genome_name <- tools::file_path_sans_ext(basename(genome_file))
  output_file <- file.path(output_directory, paste(genome_name, ".fna", sep = ""))
  
  # Open the output file for writing
  cat("", file = output_file)  # Create an empty file
  
  # Loop over each sublist in the current genome
  apply_sliding_window(sublist_names = names(genome_sequence),
                       sequences = genome_sequence,
                       output_file = output_file)
  
  # Print the output file path for reference
  cat("Result for", genome_name, "exported to", output_file, "\n")
}
