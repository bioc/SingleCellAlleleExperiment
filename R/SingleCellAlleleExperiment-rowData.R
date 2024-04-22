
#' rowData setter for the SingleCellAlleleExperiment class
#'
#' @description
#' Setter function for the rowData slot for the \code{\link{SingleCellAlleleExperiment}} class.
#'
#' @details
#' If you set `rowData(scae)<- NULL` the mandatory columns "NI_I" and "Quant_type"
#' will be kept silently, setting all other columns to NULL.
#'
#' If you want to change the content of the mandatory "NI_I" and "Quant_type" columns
#' check the valid values:
#'
#' - NI_I: c("NI" and "I") are valid values.
#' - Quant_type: c("A", "G" "F") are valid values.
#'
#' @seealso [SingleCellAlleleExperiment]
#'
#' @param x A \code{\link{SingleCellAlleleExperiment}} object
#' @param value Value of valid type and content (see validty.R)
#' @importFrom SummarizedExperiment "rowData<-"
#' @importFrom methods validObject
#' @return A \code{\link{SingleCellAlleleExperiment}} object
#' @export
setReplaceMethod("rowData",
                 "SingleCellAlleleExperiment",
                 function(x, value) {
    old <- rowData(x)
    sce <- as(x, "SingleCellExperiment")
    rowData(sce) <- value

    no_ni_i_col <- !"NI_I" %in% colnames(rowData(sce))
    no_quant_type_col <- !"Quant_type" %in% colnames(rowData(sce))

    if (no_ni_i_col && no_quant_type_col){
      stop(sprintf("\n%s is an invalid value.
                   This would delete the mandatory columns \"NI_I\" and \"Quant_type\".
                   Consider a different value or specify the columns.
                   Valid values for \"NI_I\": c(\"I\", \"NI\").
                   Valid values for \"Quant_type\": c(c(\"A\", \"G\", \"F\")",
                   value))
    }else if (no_ni_i_col) {
      stop("\nCan not set \"NI_I\" to NULL. Mandatory column needs to be kept.")
    }else if (no_quant_type_col) {
      stop("Can not set \"Quant_type\" to NULL. Mandatory column needs to be kept.")
    }

    scae <- as(sce, "SingleCellAlleleExperiment")
    methods::validObject(scae)
    scae
 })

#' rowData-NULL-setter for the SingleCellAlleleExperiment class
#'
#' @description
#' Setter function for the rowData slot for the \code{\link{SingleCellAlleleExperiment}} class.
#'
#' @param x A \code{\link{SingleCellAlleleExperiment}} object
#' @param value NULL
#' @importFrom SummarizedExperiment "rowData<-"
#' @importFrom methods validObject
#' @return A \code{\link{SingleCellAlleleExperiment}} object
#' @export
setReplaceMethod("rowData",
                 c("SingleCellAlleleExperiment", "NULL"),
                 function(x, value) {
  if (is.null(value)){
      retain_cols <- c("NI_I", "Quant_type")
      old <- rowData(x)
      value <- old[, retain_cols]
      rowData(x) <- value
  }
   methods::validObject(x)
 return(x)
})
