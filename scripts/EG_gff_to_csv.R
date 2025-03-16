# this script converts a 9 column gff from Earl Grey to a 7 column csv 
# for ease of downstream analysis

# output csv schema (egcsv):
# V1: sequence name
# V2: "EarlGrey"
# V3: start
# V4: end
# V5: type
# V6: subtype
# V7: family # OR motif
# V8: if V5 = "LINE-dependent" or "Non-TE," additional info stored here

##########################################################################

# load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(Biostrings)
library(yaml)

##########################################################################

# gather input
config <- yaml.load_file("config.yaml")

input_gff <- config$EG_gff

output_csv <- sub("input/", "output/", input_gff)
output_csv <- sub("\\.gff$", ".csv", output_csv)

seq_name <- config$seq_name

##########################################################################
  
clean_eggff <- function(eggff, egcsv, name) {
  # convert gff to data frame
  egdf <- as.data.frame(read.csv(eggff, sep = "\t", header = 0))
  
  # set all values in V1 to name
  egdf$V1 <- name
  
  # set all values in V2 to "EarlGrey"
  egdf$V2 <- "EarlGrey"
  
  # save feature, start and end columns
  feature <- egdf$V3
  start <- egdf$V4
  end <- egdf$V5
  
  # shift start and end locations to V3 and V4
  egdf$V3 <- start
  egdf$V4 <- end
  egdf$V5 <- feature
  
  # clean V5 and V6 such that they are type and subtype, respectively
  egdf <- egdf %>%
    mutate(V5_split = str_split(V5, "/")) %>%
    mutate(
      new_col1 = map_chr(V5_split, ~ .x[1]),
      new_col2 = map_chr(V5_split, ~ ifelse(length(.x) > 1, .x[2], NA))
    ) %>%
    select(-V5_split)
  
  type <- egdf$new_col1
  subtype <- egdf$new_col2
  
  egdf$V5 <- type
  egdf$V6 <- subtype
  
  # clean V7 (note that simple repeats don't have family #s)
  egdf <- egdf %>%
    mutate(V7 = as.numeric(str_extract(V9, "(?<=FAMILY-)\\d+")))
  
  # delete unnecessary columns
  egdf <- egdf %>%
    select(-V8, -V9, -new_col1, -new_col2)
  
  # store in V8 what we're about to write over in V5
  egdf$V8 <- egdf$V5
  
  # replace SINE and Retroposon with LINE-dependent
  egdf$V5 <- gsub("^(SINE|Retroposon)$", "LINE-dependent", egdf$V5)
  
  # replace RC with DNA
  egdf$V5 <- gsub("^(RC)$", "DNA", egdf$V5)
  
  # replace Simple_repeat, Low_complexity, and Satellite with Non-TE
  egdf$V5 <- gsub("^(Simple_repeat|Low_complexity|Satellite)$", "Non-TE", egdf$V5)
  
  # raise error if there are remaining types not in our classification scheme
  invalid_values <- egdf$V5[!egdf$V5 %in% c("Unknown", "Non-TE", "DNA", "LTR", "LINE", "LINE-dependent")]
  if (length(invalid_values) > 0) {
    stop("Invalid values found in V5 column: ", paste(unique(invalid_values), collapse = ", "))
  }
  
  write.table(egdf, file = egcsv, row.names = FALSE, col.names = FALSE, sep = ",")
}

clean_eggff(input_gff, output_csv, seq_name)