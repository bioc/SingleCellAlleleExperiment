#' @rdname SingleCellAlleleExperiment
#' 
#' @exportClass SingleCellAlleleExperiment SingleCellAlleleExperiment
#' 
#' @importFrom SingleCellExperiment SingleCellExperiment
setClass("SingleCellAlleleExperiment",
         contains = "SingleCellExperiment")
