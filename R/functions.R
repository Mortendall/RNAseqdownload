library(GEOquery)
library(tidyverse)
library(fs)
library(vroom)
library(Biobase)

#' Download GEO file
#'
#' @param GEOID An ID for the GEO database https://www.ncbi.nlm.nih.gov/geo/
#'
#' @return A table with all the available metadata. Can be used to select relevant parameters for the metadata file.

download_GEO_file <- function(GEOID){
  gds <- GEOquery::getGEO(GEOID ,GSEMatrix = T, getGPL = F)
  GEO_data <- Biobase::pData(gds[[1]])
  return(GEO_data)
}

#' Generate metadata
#'
#' @param GEO_data A data.frame containing the information about the experiment.
#'
#' @param characteristics_list  vector with the terms you wish to include in your metadata-file.
#' Can be things like group, ID, gender. Check the output of the function download_GEO_file.
#'
#' @return

generate_metadata <- function(GEO_data, characteristics_list){
  metadata <- GEO_data[characteristics_list]
  write.table(metadata, file = "metadata.csv", )
  return(metadata)
}

#' GEOdownloader
#'
#' @param GEOID An ID for the GEO database https://www.ncbi.nlm.nih.gov/geo/
#' @param metadata a metadata table containing relevant parameters for the experiment
#'
#' @return a count matrix containing raw count values for the experiment with group annotations.

GEODownloader <- function(GEOID, metadata){
  geo_download <- GEOquery::getGEOSuppFiles(GEOID, makeDirectory = F,baseDir = here("data-raw/"))

  geo_file <- fs::dir_ls(here("data-raw/"),
                     regexp = "txt",
                     recurse = TRUE)

  count_matrix <- vroom::vroom(geo_file)
  matrix_key <- count_matrix$key
  count_matrix <- count_matrix %>%
    dplyr::select(-1)
  rownames(count_matrix) <- matrix_key
 return(count_matrix)
}

  column_order = t(count_matrix[1,])
  column_order <- as.data.frame(column_order)
  column_order[1,]<- as.character(column_order[1,])

#' ID annotator
#'
#' @param metadata your metadata file
#' @param ID Sample ID
#' @param Group Group/condition for the sample
#'
#' @return returns ID mapped to disease condition. Can be used to organize the alignment
ID_key_generator <- function(metadata, ID, Group){
  column_ID <- metadata %>%
    select(c(ID, Group))
  return(column_ID)
}


