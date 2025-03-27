# this script converts a 9 column gff from Garlic to a 7 column csv 
# for ease of downstream analysis

# output csv schema (garliccsv):
# V1: sequence name
# V2: "Reference (Garlic)"
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

input_gff <- config$Garlic_gff

output_csv <- sub("input/", "output/", input_gff)
output_csv <- sub("\\.gff$", ".csv", output_csv)

seq_name <- config$seq_name

##########################################################################

clean_garlicgff <- function(garlicgff, garliccsv, name) {
  # convert gff to data frame
  garlicdf <- as.data.frame(read.csv(garlicgff, sep = "\t", header = 0))
  
  # set all values in V1 to name
  garlicdf$V1 <- name
  
  # set all values in V2 to Reference (Garlic)"
  garlicdf$V2 <- "Reference (Garlic)"
  
  # save feature, start and end columns
  feature <- garlicdf$V3
  start <- garlicdf$V4
  end <- garlicdf$V5
  
  # shift start and end locations to V3 and V4
  garlicdf$V3 <- start
  garlicdf$V4 <- end
  garlicdf$V5 <- feature
  
  # for all rows with no information on subtype, set V6 to NA
  rows_wo_subclass <- !grepl("/", garlicdf$V5)
  garlicdf$V7[rows_wo_subclass] <- str_extract(garlicdf$V5[rows_wo_subclass], "^[^:]+")
  garlicdf$V5[rows_wo_subclass] <- str_extract(garlicdf$V5[rows_wo_subclass], "(?<=:).*")
  garlicdf$V6[rows_wo_subclass] <- NA
  
  # clean the rest of the rows such that V5 is type, V6 is subtype, 
  # and V7 is family
  rows_w_colon <- grepl(":", garlicdf$V5)
  garlicdf <- garlicdf %>%
    mutate(
      family = str_extract(V5, "^[^:]+"),
      type = str_extract(V5, "(?<=:)[^/]+"),
      subtype = str_extract(V5, "(?<=/).*")
    ) %>%
    mutate(
      V7 = ifelse(grepl(":", V5), family, V7),
      V6 = ifelse(grepl(":", V5), subtype, V6),
      V5 = ifelse(grepl(":", V5), type, V5)
    )
  
  # delete unnecessary columns
  garlicdf <- garlicdf %>%
    select(-V8, -V9, -family, -type, -subtype)
  
  # store in V8 what we're about to write over in V5
  garlicdf$V8 <- garlicdf$V5
  
  # replace SINE and Retroposon with LINE-dependent
  garlicdf$V5 <- gsub("^(SINE|Retroposon)$", "LINE-dependent", garlicdf$V5)
  
  # replace tRNA, snRNA, rRNA, RNA, and scRNA with Non-TE
  garlicdf$V5 <- gsub("^(tRNA|snRNA|rRNA|RNA|scRNA)$", "Non-TE", garlicdf$V5)
  
  # replace RC with DNA
  garlicdf$V5 <- gsub("^(RC)$", "DNA", garlicdf$V5)
  
  # replace Other YPRIME with DNA YPRIME
  garlicdf$V5[garlicdf$V7 == "YPRIME"] <- "DNA"
  
  # replace Other NTS_DM with Non-TE NTS_DM
  garlicdf$V5[garlicdf$V7 == "NTS_DM"] <- "Non-TE"
  
  # replace values containing "?" with Unknown
  garlicdf$V5 <- ifelse(grepl("\\?", garlicdf$V5), "Unknown", garlicdf$V5)
  
  # raise error if there are remaining types not in our classification scheme
  invalid_values <- garlicdf$V5[!garlicdf$V5 %in% c("Unknown", "Non-TE", "DNA", "LTR", "LINE", "LINE-dependent")]
  if (length(invalid_values) > 0) {
    stop("Invalid values found in V5 column: ", paste(unique(invalid_values), collapse = ", "))
  }
  
  garlicdf <- garlicdf[garlicdf$V4 <= 100000000, ]
  
  write.table(garlicdf, file = garliccsv, row.names = FALSE, col.names = FALSE, sep = ",")

}

clean_garlicgff(input_gff, output_csv, seq_name)