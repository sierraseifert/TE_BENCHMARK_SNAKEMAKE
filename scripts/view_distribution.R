# Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(yaml)

# Function to parse the input from the file
parse_input_from_file <- function(file_path) {
  # Read the file content
  lines <- readLines(file_path)
  
  # Define the bins for coverage
  bins <- c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%")
  
  # Parse the lines for each pipeline
  pipeline_data <- lapply(lines, function(line) {
    # Extract pipeline name (first word of the line)
    pipeline_name <- strsplit(line, ":")[[1]][1]
    
    # Extract numbers (separated by commas) and convert to numeric
    counts <- as.numeric(strsplit(strsplit(line, ":")[[1]][2], ",")[[1]])
    
    # Create a data frame for each pipeline
    df <- data.frame(
      Bin = bins,
      Count = counts,
      Pipeline = rep(pipeline_name, length(bins))
    )
    
    return(df)
  })
  
  # Combine all data frames into one
  df_combined <- do.call(rbind, pipeline_data)
  
  return(df_combined)
}

# Read the input from the 'coverage.txt' file
df_combined <- parse_input_from_file("output/coverage.txt")

# Create the plot
plot <- ggplot(df_combined, aes(x = Bin, y = Count, fill = Pipeline, group = Pipeline)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("RepeatModeler2" = "red", "Earl Grey" = "blue", "EDTA" = "green")) +
  labs(title = "Test Pipeline Coverage Breakdown",
       x = "Coverage Bin (%)", y = "Number of TEs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot to a file
png("output/plot_distribution.png", width = 3000, height = 2100, res = 300)
print(plot)  
dev.off()
