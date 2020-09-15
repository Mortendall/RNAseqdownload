## code to prepare `DATASET` dataset goes here
library(here)

# GEO_data_file <- download_GEO_file("GSE126848")
# #selects 3 characteristica, "gender:ch1", "disease:ch1" and "description"(which contains numerical ID)
#
# metadata_file <- generate_metadata(c("gender:ch1", "disease:ch1", "description"))
# count_matrix <- count_matrix_assembly("txt")

#identified ID in count matrix as column "description" in metadata
ID_key_file <- ID_key_generator(metadata_file,"description", "disease:ch1")


