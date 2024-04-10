
#' Subset SCAE object
#'
#' @description
#' Function used for subsetting the different layers stored in a
#' SingleCellAlleleExperiment object. Valid subset values are:
#' `subset=c("nonimmune", "alleles", "immune_genes", "functional_groups")`.
#'
#' @param scae SCAE object
#' @param subset Character string specifying a data layer. Valid values are
#' `subset=c("nonimmune", "alleles", "immune_genes", "functional_groups")`.
#' @importFrom methods validObject is
#' @return SCAE object
#' @examples
#' example_data_5k <- scaeData::scaeDataGet(dataset="pbmc_5k")
#' lookup_name <- "pbmc_5k_lookup_table.csv"
#' lookup <- read.csv(system.file("extdata", lookup_name, package="scaeData"))
#'
#' scae <- read_allele_counts(example_data_5k$dir,
#'                           sample_names="example_data_wta",
#'                           filter_mode="no",
#'                           lookup_file=lookup,
#'                           barcode_file=example_data_5k$barcodes,
#'                           gene_file=example_data_5k$features,
#'                           matrix_file=example_data_5k$matrix,
#'                           filter_threshold=0,
#'                           verbose=TRUE)
#'
#' scae
#'
#' scae_nonimmune_subset <- scae_subset(scae, subset="nonimmune")
#' scae_nonimmune_subset
#'
#' scae_alleles_subset <- scae_subset(scae, subset="alleles")
#' scae_alleles_subset
#'
#' scae_immune_genes_subset <- scae_subset(scae, subset="immune_genes")
#' scae_immune_genes_subset
#'
#' scae_functional_groups_subset <- scae_subset(scae, subset="functional_groups")
#' scae_functional_groups_subset
#'
#' @export
scae_subset <- function(scae, subset=c("nonimmune", "alleles", "immune_genes", "functional_groups")){

  if (is.character(subset) & is(scae, 'SingleCellAlleleExperiment')){
    scae_sub<- switch(subset,
                      "nonimmune"=get_nigenes(scae),
                      "alleles"=scae_subset_alleles(scae),
                      "immune_genes"=get_agenes(scae),
                      "functional_groups"=scae_subset_functional(scae),
                      message("Invalid layer specified, Choose from `nonimmune`, `alleles`, `immune_genes`, `functional_groups`"))
    methods::validObject(scae_sub)
    return(scae_sub)
  }else{
    stop("Input must be a character string,  Choose from `nonimmune`, `alleles`, `immune_genes`, `functional_groups`")
  }
}


##--------------------Internal functions used in scae_subset()----------------##

#' Get allele rows
#'
#' @description
#' Getter function returning subsampled SCAE object with all rows containing
#' raw allele information. These rows are identified by "I" in
#' rowData(scae)$NI_I and "A" in rowData(scae)$Quant_type.
#'
#' @param scae A \code{\link{SingleCellAlleleExperiment}} object.
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#' @return A SingleCellAlleleExperiment object.
scae_subset_alleles <- function(scae) {
  rd_qt <- rowData(scae)$Quant_type
  rd_ni_i <- rowData(scae)$NI_I
  subset_rows <- stats::complete.cases(rd_ni_i, rd_qt)
  ## alleles of the genes with extended quantification
  alleles <- scae[subset_rows & rd_ni_i == "I" & startsWith(rd_qt,"A"), ]
  return(alleles)
}

#' Get immune gene rows
#'
#' @description
#' Getter function returning subsampled SCAE object with all rows containing
#' immune gene information. These rows are identified by "I" in
#' rowData(scae)$NI_I and "G" in rowData(scae)$Quant_type.
#'
#' @param scae A \code{\link{SingleCellAlleleExperiment}} object.
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#' @return A SingleCellAlleleExperiment object.
get_agenes <- function(scae) {
  rd_qt <- rowData(scae)$Quant_type
  rd_ni_i <- rowData(scae)$NI_I
  subset_rows <- stats::complete.cases(rd_ni_i, rd_qt)
  ## genes with extended quantification
  agenes <- scae[subset_rows & rd_ni_i == "I" & rd_qt == "G", ]
  return(agenes)
}

#' Get non-immune rows
#'
#' @description
#' Getter function returning subsampled SCAE object with all rows containing
#' non-immune gene information. These rows are identified by "NI" in
#' rowData(scae)$NI_I and "G" in rowData(scae)$Quant_type.
#'
#' @param scae A \code{\link{SingleCellAlleleExperiment}} object.
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#' @return A SingleCellAlleleExperiment object.
get_nigenes <- function(scae) {
  rd_qt <- rowData(scae)$Quant_type
  rd_ni_i <- rowData(scae)$NI_I
  subset_rows <- stats::complete.cases(rd_ni_i, rd_qt)
  ## classical genes
  agenes <- scae[subset_rows & rd_ni_i == "NI" & rd_qt == "G", ]
  return(agenes)
}

#' Get functional class rows
#'
#' @description
#' Getter function returning subsampled SCAE object with all rows containing
#' functional class information. These rows are identified by "I" in
#' rowData(scae)$NI_I and "F" in rowData(scae)$Quant_type.
#'
#' @param scae A \code{\link{SingleCellAlleleExperiment}} object.
#' @importFrom stats complete.cases
#' @importFrom SingleCellExperiment rowData
#' @return A SingleCellAlleleExperiment object.
scae_subset_functional <- function(scae) {
  rd_qt <- rowData(scae)$Quant_type
  rd_ni_i <- rowData(scae)$NI_I
  subset_rows <- stats::complete.cases(rd_ni_i, rd_qt)
  # functional groups of the genes with extended quantification
  func <- scae[subset_rows & rd_ni_i == "I" & rd_qt == "F", ]
  return(func)
}
