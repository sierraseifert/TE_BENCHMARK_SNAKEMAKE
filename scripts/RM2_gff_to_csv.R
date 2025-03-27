# this script converts a 9 column gff from RepeatModeler2 to a 7 column csv 
# for ease of downstream analysis

# output csv schema (rmcsv):
# V1: sequence name
# V2: "RepeatModeler2"
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

input_gff <- config$RM2_gff

input_fasta <- config$RM2_fasta

output_csv <- sub("input/", "output/", input_gff)
output_csv <- sub("\\.gff$", ".csv", output_csv)
  
seq_name <- config$seq_name

##########################################################################

clean_rmgff <- function(rmgff, rmfasta, rmcsv, name) {
  # convert gff to data frame
  rmdf <- as.data.frame(read.csv(rmgff, sep = "\t", header = 0))

  # delete top 2 rows containing metadata
  rmdf <- rmdf[-c(1, 2), ]

  # set all values in V1 to name
  rmdf$V1 <- name

  # set all values in V2 to "RepeatModeler2"
  rmdf$V2 <- "RepeatModeler2"

  # keep V3 data in case (dispersed_repeat value, anything else?)
  keep <- rmdf$V3

  # shift start and end locations to V3 and V4
  rmdf$V3 <- rmdf$V4
  rmdf$V4 <- rmdf$V5

  # condense V9 values to family numbers to prepare for merge
  rmdf$V9 <- sub(".*Target Motif:([^ ]+).*", "\\1", rmdf$V9)

  # read fasta file to extract classification info by merging on family number
  seqs <- readDNAStringSet(rmfasta)
  titles <- names(seqs)
  split_titles <- strsplit(titles, "#")
  fastadf <- data.frame(
    number = sapply(split_titles, function(x) trimws(x[1])),
    type = sapply(split_titles, function(x) trimws(x[2])),
    stringsAsFactors = FALSE
  )

  # merge fastadf and rmdf on family number
  fastadf$type <- sub("\\s*\\(.*", "", fastadf$type)
  rmdf <- merge(rmdf, fastadf, by.x = "V9", by.y = "number", all.x = TRUE)

  # move type to V5
  rmdf$V5 <- rmdf$type

  # delete unnecessary columns and re order
  rmdf <- rmdf %>%
    select(-V9, -V7, -V8, -type, V9)
  colnames(rmdf)[which(colnames(rmdf) == "V9")] <- "V7"

  # if subtype exists, move it to V6
  rmdf$V6 <- NA
  rows_with_subtype <- grepl("/", rmdf$V5)
  rmdf$V6[rows_with_subtype] <- sub(".*?/", "", rmdf$V5[rows_with_subtype])
  rmdf$V5[rows_with_subtype] <- sub("/.*", "", rmdf$V5[rows_with_subtype])

  # store in V8 what we're about to write over in V5
  rmdf$V8 <- rmdf$V5

  # replace tRNA, rRNA, snRNA, scRNA, Simple_repeat, and Satellite with Non-TE
  rmdf$V5[is.na(rmdf$V5)] <- "Non-TE"
  rmdf$V5 <- gsub("^(tRNA|rRNA|snRNA|scRNA|Simple_repeat|Satellite)$", "Non-TE", rmdf$V5)

  # replace SINE and Retroposon with LINE-dependent
  rmdf$V5 <- gsub("^(SINE|Retroposon)$", "LINE-dependent", rmdf$V5)

  # replace RC with DNA
  rmdf$V5 <- gsub("^(RC)$", "DNA", rmdf$V5)

  # replace values containing "?" with Unknown
  rmdf$V5 <- ifelse(grepl("\\?", rmdf$V5), "Unknown", rmdf$V5)

  # raise error if there are remaining types not in our classification scheme
  invalid_values <- rmdf$V5[!rmdf$V5 %in% c("Non-TE", "Unknown", "DNA", "LTR", "LINE", "LINE-dependent")]
  if (length(invalid_values) > 0) {
    stop("STOP! There are values in the V5 column that are not in your
         classification scheme: ",
         paste(unique(invalid_values), collapse = ", "))
  }

  write.table(rmdf, file = rmcsv, row.names = FALSE, col.names = FALSE, sep = ",")
 }

clean_rmgff(input_gff, input_fasta, output_csv, seq_name)
