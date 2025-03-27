# this script filters an input csv by a specified target value and saves
# it as a new csv

###USER INPUT REQUIRED####################################################
# gather input
input_csv <- snakemake@input[['csv']] # path to input csv

output_csv <- snakemake@output[['csv']] # path to output csv

column <- "V5" # column in input_csv containing target value

target <- snakemake@params['te_type'] # target value in column

##########################################################################

# load libraries
library(dplyr)
library(tidyr)

# convert input csv to data frame
df <- read.csv(input_csv, header = TRUE)

# filter data frame
filtered_df <- df %>%
  filter(df[[column]] == target)

# save filtered data frame as output_csv
write.csv(filtered_df, file = output_csv, row.names = FALSE)
