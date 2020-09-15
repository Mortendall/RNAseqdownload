library(GEOquery)
library(tidyverse)
library(fs)
library(vroom)
library(Biobase)
library(here)

#' Download GEO files
#'
#' @param GEOID An ID for the GEO database https://www.ncbi.nlm.nih.gov/geo/
#'
#' @return A table with all the available metadata. Can be used to select relevant parameters for the metadata file.


download_GEO_file <- function(GEOID){
  gds <- GEOquery::getGEO(GEOID ,GSEMatrix = T, getGPL = F)
  GEO_data <- Biobase::pData(gds[[1]])
  utils::write.table(GEO_data, file = here("data-raw/GEO_data.csv"), row.names = F)
  GEO_counts <- GEOquery::getGEOSuppFiles(GEOID, makeDirectory = F, baseDir = here("data-raw/"))
  return(GEO_data)
}

#' Generate metadata
#'
#' @param characteristics_list  vector with the terms you wish to include in your metadata-file.
#' Can be things like group, ID, gender. Check the output of the function download_GEO_file.
#'
#' @return

generate_metadata <- function(characteristics_list){
  GEO_data <- vroom::vroom(here("data-raw/GEO_data.csv"))
  metadata <- GEO_data[characteristics_list]
  utils::write.table(metadata, file = here("data-raw/metadata.csv"), row.names = F)
  return(metadata)
}

#' count_matrix_loader
#'
#' @param file_type your raw data file type (presently optimized for txt)
#'
#' @return a count matrix containing raw count values for the experiment with group annotations.

count_matrix_assembly <- function(file_type){
  geo_file <- fs::dir_ls(here("data-raw/"),
                     regexp = file_type,
                     recurse = TRUE)

  count_matrix <- vroom::vroom(geo_file)
  matrix_key <- count_matrix$key
  count_matrix <- count_matrix %>%
    dplyr::select(-1)
  rownames(count_matrix) <- matrix_key
  utils::write.table(count_matrix,
              file = here("data-raw/count_matrix.csv"),
              row.names = T,
              col.names = T)
 return(count_matrix)
}


#' ID annotator
#'
#' @param metadata your metadata file
#' @param ID Sample ID
#' @param Group Group/condition for the sample
#'
#' @return returns ID mapped to disease condition. Can be used to organize the alignment

ID_key_generator <- function(ID, Group){
  metadata <- vroom::vroom(here("data-raw/metadata.csv"),
                           delim = " ")
  column_ID <- metadata %>%
    select(c(ID, Group))
  return(column_ID)
}


