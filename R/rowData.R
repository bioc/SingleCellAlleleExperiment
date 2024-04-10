
#' Setting rowData
#'
#' @description
#' Function used to change and set new values in the `rowData` slot.
#' Incorporates the validity checks on the output SCAE object to catch
#' invalid changes.
#'
#' @param x A \code{\link{SingleCellAlleleExperiment}} object
#' @param ... additional parameters passed to rowData<-
#' @param value Value of valid type and content (see validty.R)
#' @importFrom SummarizedExperiment "rowData<-"
#' @importFrom methods validObject callNextMethod
#' @return A \code{\link{SingleCellAlleleExperiment}} object
#' @export
setReplaceMethod("rowData", "SingleCellAlleleExperiment", function(x, ..., value) {
  scae <- callNextMethod()
  methods::validObject(scae)
})
