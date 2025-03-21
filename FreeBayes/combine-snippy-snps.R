### Combine snippy filtered SNP files ###

library(dplyr)
library(tidyr)
library(purrr)

# Read the reference SNP file with SNP names and SNP position
ref_file <- read.delim("../cheetah_SNP-positions_B.txt", header = FALSE, sep = "\t", col.names = c("SNP", "Position"))

# List all TSV files in the directory
tsv_files <- list.files(pattern = "*.tsv")

# Function to read and rename columns
read_snp_file <- function(file) {
  df <- read.delim(file, header = TRUE, sep = "\t")
  sample_name <- tools::file_path_sans_ext(basename(file))
  names(df)[names(df) == "Evidence"] <- sample_name
  return(df)
}

# Read all SNP files into a list
snp_data_list <- map(tsv_files, read_snp_file)

# Merge all SNP files with the reference file
merged_data <- reduce(snp_data_list, full_join, by = c("SNP", "Position"))

# Ensure all SNPs in the reference file are present
final_data <- full_join(ref_file, merged_data, by = c("SNP", "Position"))

# Replace missing values with NA
final_data[is.na(final_data)] <- "NA"

# Write output using base R
write.table(final_data, "merged_snps.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
