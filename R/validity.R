validSingleCellAlleleExperiment <- function(scae_object) {
  msg <- NULL

  if (!"NI_I" %in% colnames(rowData(scae_object))){
    msg <- c(msg, "No column 'NI_I' found in the rowData slot.")
  }

  if (!"Quant_type" %in% colnames(rowData(scae_object))){
    msg <- c(msg, "No column 'Quant_type' found in the rowData slot")
  }

  if (sum(!rowData(scae_object)$NI_I %in% c("NI", "I"))>0){
    msg <- c(msg, "Extended rowData column 'NI_I' contains invalid values")
  }

  if (sum(!rowData(scae_object)$Quant_type %in% c("A", "G", "F"))>0){
    msg <- c(msg, "Extended rowData column 'Quant_type' contains invalid values")
  }

  if (is.null(msg)) {
    TRUE
  } else msg

}

#' @importFrom S4Vectors setValidity2
S4Vectors::setValidity2("SingleCellAlleleExperiment",
                        validSingleCellAlleleExperiment)
