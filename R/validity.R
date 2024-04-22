validSingleCellAlleleExperiment <- function(scae_object) {
  msg <- NULL

  if (!"NI_I" %in% colnames(rowData(scae_object))){
    msg <- c(msg, "No column 'NI_I' found in the rowData slot.")
  }

  if (!"Quant_type" %in% colnames(rowData(scae_object))){
    msg <- c(msg, "No column 'Quant_type' found in the rowData slot")
  }

  if (sum(!rowData(scae_object)$NI_I %in% c("NI", "I"))>0){
    msg <- c(msg, "Extended rowData column 'NI_I' contains invalid values. Only \"NI\" and \"I\" are valid.")
  }

  if (sum(!rowData(scae_object)$Quant_type %in% c("A", "G", "F"))>0){
    msg <- c(msg, "Extended rowData column 'Quant_type' contains invalid values. Only \"A\", \"G\" and \"F\" are valid.")
  }

  if (is.null(msg)) {
    TRUE
  } else msg
}

#' @importFrom S4Vectors setValidity2
S4Vectors::setValidity2("SingleCellAlleleExperiment",
                        validSingleCellAlleleExperiment)


#' Check package installation for optional functionalities
#'
#' @param log A logical parameter to decide if a logcounts assay should be computed
#' based on library factors computed with `scuttle::computeLibraryFactors()`.
#' @param gene_symbols A logical parameter to decide whether to compute additional
#' gene gene symbols in case the raw data only contains ENSEMBL gene identifiers.
#' @return Error messages if cases are met
check_valid_optional_package <- function(log, gene_symbols){
  if (log){
    if (!requireNamespace("scuttle", quietly=TRUE)) {
      stop("Package 'scuttle' needed when using 'log=TRUE'.
         Install: BiocManager::install(\"scuttle\")")
    }
  }
  if (gene_symbols){
    if (!requireNamespace("org.Hs.eg.db", quietly=TRUE)){
      stop("Package \"org.Hs.eg.db\" needed when using 'gene_symbols=TRUE'.
              Install: BiocManager::install(\"org.Hs.eg.db\")")
    }
    if (!requireNamespace("AnnotationDbi", quietly=TRUE)){
      stop("Package \"AnnotationDbi\" needed when using 'gene_symbols=TRUE'.
              Install: BiocManager::install(\"org.Hs.eg.db\")")
    }
  }
}

