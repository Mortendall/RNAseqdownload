library(GEOquery)
library(tidyverse)
library(fs)
library(vroom)
library(Biobase)
library(here)
library(GEOquery)
library(tidyverse)
library(fs)
library(vroom)
library(Biobase)
library(here)
library(magrittr)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(edgeR)
library(openxlsx)


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
  utils::write.table(count_matrix,
                     file = here("data-raw/count_matrix.csv"),
                     row.names = F,
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


#' Count_matrix_loader
#'
#' @param file - name of the count matrix file
#'
#' @return a count matrix tibble with gene IDs as row names

count_matrix_load <- function(file) {
  data_file <- fs::dir_ls(here::here("data-raw/"),
                          regexp = file,
                          recurse = T)
  count_matrix <- vroom::vroom(data_file)
  colnames(count_matrix) <- str_remove(colnames(count_matrix), "^0+")
  count_matrix_key <- count_matrix$key
  count_matrix <- count_matrix %>%
    dplyr::select(-key)
  rownames(count_matrix) = count_matrix_key
  return(count_matrix)
}

#' Load metadata and sort them according to count matrix input
#'
#' @param file_name the name of teh metadata file (default "metadata.csv")
#' @param count_matrix the name of the count matrix object
#'
#' @return metadata file sorted according to count matrix order

load_metadata <- function(file_name, count_matrix) {
  data_file <- fs::dir_ls(here::here("data-raw/"),
                          regexp = file_name,
                          recurse = T)
  metadata <-
    vroom::vroom(
      data_file,
      col_types = cols(
        `gender:ch1` = col_character(),
        `disease:ch1` = col_character(),
        description = col_character()
      )
    )
  data_order <- colnames(count_matrix)
  data_order <- as_tibble(data_order)
  colnames(data_order) = "description"
  metadata <- full_join(data_order, metadata, "description")
  return(metadata)
}

#' Generate design matrix
#'
#' @param metadata a metadata object generated through the load_metadata function
#'
#' @return a design matrix file

Generate_design_matrix <- function(metadata){
  group <- as.matrix(metadata[3])
  design <- model.matrix( ~ 0 + group, metadata)
  colnames(design) <-
    str_remove_all(colnames(design), "\\(|\\)|group|:")
  return(design)
}

#' RNAseq analysis
#'
#' @param count_matrix count matrix generated from count_matrix_load
#' @param design design matrix made through the Generate_design_matrix
#' @param metadata metadata generated from load_metadata function
#' @param contrast_matrix generated through the limma::makeContrast function
#'
#' @return

RNAseq_processing <- function(count_matrix, metadata, design, contrast_matrix) {
  group <- as.matrix(metadata[3])
  RNAseq <- DGEList(counts = count_matrix, group = group)
  keep <- filterByExpr(RNAseq)
  RNAseq <- RNAseq[keep, , keep.lib.sizes = F]
  RNAseq <- calcNormFactors(RNAseq)
  RNAseq <- estimateDisp(RNAseq,design, robust = T)
  et <- exactTest(RNAseq)
  efit <- glmQLFit(RNAseq, design, robust = T)
  dgeResults <- apply(cont.matrix, 2, . %>%
                        glmQLFTest(glmfit = efit, contrast = .) %>%
                        topTags(n = Inf, p.value = 1) %>%
                        extract2("table") %>%
                        as.data.table(keep.rownames = TRUE))
  return(dgeResults)
}

#' Annotate results
#'
#' @param Results_file an object containing ensembl IDs in the first column
#'
#' @return A list annotated with gene symbols

annotated_dgeResults <- function(Results_file) {
  dat <- data.table::copy(Results_file)
  data.table::setnames(dat, names(dat)[1], "ENSEMBL")
  ens2symbol <-
    clusterProfiler::bitr(dat$ENSEMBL,
         fromType = 'ENSEMBL',
         toType = 'SYMBOL',
         OrgDb = "org.Hs.eg.db")
  ens2symbol %<>% data.table %>% data.table::setkey(ENSEMBL)
  dat <- dplyr::full_join(dat, ens2symbol)
  return(dat)
}

#' Gene ontology enrichment analysis of genes generated from a results file
#'
#' @param result_list list of data.tables generated from edgeR. Must be data.table and contain a SYMBOL annotation column
#'
#' @return a list containing enrichresults for each element in the results file list

goAnalysis <- function(result_list){
  bg_list <- result_list[[1]][,"SYMBOL"]
  bg_list = clusterProfiler::bitr(
    bg_list$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db",
    drop = T
  )
  goResult_list <- vector(mode = "list", length = length(result_list))
  for(i in 1:length(result_list)){
    sig_list<- result_list[[i]] %>%
      filter(FDR<0.05)

    eg <- clusterProfiler::bitr(
      sig_list$SYMBOL,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = "org.Hs.eg.db",
      drop = T
    )
    goResults <- enrichGO(gene = eg$ENTREZID,
                          universe = bg$ENTREZID,
                          OrgDb = org.Hs.eg.db,
                          ont = "BP")
    goResult_list[[i]]<- goResults
  }
  for (i in 1:length(goResult_list)){
    names(goResult_list)[i]<-names(result_list)[i]
  }
  return(goResult_list)

}

#' File exporter - exports GOresults as an excel sheet, and prints dotplot and cnet plots
#'
#' @param goList a list object containing one or more enrichResults
#'
#' @return

printGOterms <- function(goList){
  goSheets<- vector(mode = "list", length = length(goList))
  for (i in 1:length(goSheets)){
    goSheets[[i]] <- goList[[i]]@result
    names(goSheets)[i]<-names(goList)[i]
  }
  write.xlsx(goSheets, file = here("data/NASH_NAFLD_GOterms.xlsx"), asTable = TRUE)
  dir.create(here("data/figures"), showWarnings = F)
  for (i in 1:length(goList)){
    dotplot <- dotplot(goList[[i]], title = names(goList)[i],size = 1)
    ggsave(dotplot, filename = paste(here("data/figures"),"/dotplot_",names(goList[i]),".png", sep = ""),width = 12, height = 12, units = "cm", scale = 2.5)
    goList_anno <- setReadable(goList[[i]], OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    cnetplot <- cnetplot(goList_anno, title = names(goList)[i], size = 1)
    ggsave(cnetplot, filename = paste(here("data/figures"),"/cnetplot_",names(goList[i]),".png", sep = ""),scale = 2.5)
  }
}

