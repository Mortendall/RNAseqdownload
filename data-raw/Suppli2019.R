# # code to prepare `DATASET` dataset goes here
#
# GEO_data_file <- download_GEO_file("GSE126848")
# #selects 3 characteristica, "gender:ch1", "disease:ch1" and "description"(which contains numerical ID)
#
# metadata_file <- generate_metadata(c("gender:ch1", "disease:ch1", "description"))
# count_matrix <- count_matrix_assembly("txt")
#
# # identified ID in count matrix as column "description" in metadata
# ID_key_file <- ID_key_generator("description", "disease:ch1")

# Load in the count matrix
count_matrix <- count_matrix_load("count_matrix.csv")
metadata <- load_metadata("metadata.csv",count_matrix)

design_matrix <- Generate_design_matrix(metadata)

#use the design matrix to make a relevant contrast matrix
cont.matrix <- limma::makeContrasts(
  NAFLD_vs_healthy = NAFLD-healthy,
  obese_vs_NAFLD = obese-NAFLD,
  NASH_vs_NAFLD = NASH-NAFLD,
  obese_vs_healthy = obese-healthy,
  levels = design_matrix)


RNAseq_data <- RNAseq_processing(count_matrix, metadata, design_matrix,cont.matrix)

results<- lapply(RNAseq_data,annotated_dgeResults)

