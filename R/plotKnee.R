
#------------------------------Knee plot---------------------------------------#

#' Knee plot
#'
#' @description
#' Creates a knee plot ranking the barcodes according to their total UMI count. The plot is
#' used to determine a threshold for filtering barcodes in the preprocessing step.
#'
#' @param matrix A sparse \code{\link{Matrix}} object containing the quantification data.
#' @param genes A data.frame object containing gene identifiers.
#' @param barcodes A data.frame object containing barcode identifiers.
#'
#' @importFrom S4Vectors metadata
#'
#' @return A list including a data.frame with barcode rank information, the corresponding knee and inflection point.
get_knee_info <- function(matrix, genes, barcodes){

  barcodes <- barcodes
  features <- genes
  matrix <- matrix

  #advanced knee plot
  br_out <- DropletUtils::barcodeRanks(matrix)
  knee_point <- S4Vectors::metadata(br_out)$knee
  inflection_point <- S4Vectors::metadata(br_out)$inflection

  knee_list <- c(knee_df=br_out, knee_point=knee_point, inflection_point=inflection_point)
  return(knee_list)
}
