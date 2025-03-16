# this script converts a 9 column gff from EDTA to a 7 column csv 
# for ease of downstream analysis

# output csv schema (edtacsv):
# V1: sequence name
# V2: "EDTA"
# V3: start
# V4: end
# V5: type
# V6: subtype
# V7: family #
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

input_gff <- config$EDTA_gff

output_csv <- sub("input/", "output/", input_gff)
output_csv <- sub("\\.gff$", ".csv", output_csv)

seq_name <- config$seq_name

##########################################################################

clean_edtagff <- function(edtagff, edtacsv, name) {
  gff_lines <- readLines(edtagff)[-c(1:6)]
  
  temp_file <- tempfile()
  writeLines(gff_lines, temp_file)
  
  # convert gff to data frame
  edtadf <- as.data.frame(read.csv(temp_file, sep = "\t", header = FALSE))
  
  unlink(temp_file)
  
  # set all values in V1 to name
  edtadf$V1 <- name
  
  # save feature, start and end columns
  feature <- edtadf$V3
  start <- edtadf$V4
  end <- edtadf$V5
  
  # shift start and end locations to V3 and V4
  edtadf$V3 <- start
  edtadf$V4 <- end
  edtadf$V5 <- feature
  
  # clean V5 and V6 such that they are type and subtype, respectively
  edtadf$V5 <- ifelse(
    grepl("/", edtadf$V9), 
    sub(".*Classification=([^/]+)\\/.*", "\\1", edtadf$V9), 
    sub(".*Classification=([^;]+);.*", "\\1", edtadf$V9) 
  )
  
  edtadf$V6 <- ifelse(
    grepl("/", edtadf$V9),  # Check if there is a '/'
    sub(".*Classification=[^/]+/([^;]+);.*", "\\1", edtadf$V9),
    NA 
  )
  
  # clean V7 such that it is the name of the TE (TEs that have the same name
  # belong to the same family?)
  edtadf$V7 <- sub(".*Name=([^;]+);.*", "\\1", edtadf$V9)
  
  # delete unnecessary columns
  edtadf <- edtadf %>%
    select(-V8, -V9)
  
  # store in V8 what we're about to write over in V5
  edtadf$V8 <- edtadf$V5
  
  # replace SINE with LINE-dependent... EDTA doesn't pick up on Retroposons?
  edtadf$V5 <- gsub("^(SINE)$", "LINE-dependent", edtadf$V5)
  
  # no Non-TE classification, just Unknown repeats
  # capitalize
  edtadf$V5 <- gsub("^(unknown)$", "Unknown", edtadf$V5)
  
  # replace MITE with DNA... is this valid?
  edtadf$V5 <- gsub("^(MITE)$", "DNA", edtadf$V5)
  
  # raise error if there are remaining types not in our classification scheme
  invalid_values <- edtadf$V5[!edtadf$V5 %in% c("Unknown", "Non-TE", "DNA", "LTR", "LINE", "LINE-dependent")]
  if (length(invalid_values) > 0) {
    stop("Invalid values found in V5 column: ", paste(unique(invalid_values), collapse = ", "))
  }
  
  write.table(edtadf, file = edtacsv, row.names = FALSE, col.names = FALSE, sep = ",")
}

clean_edtagff(input_gff, output_csv, seq_name)