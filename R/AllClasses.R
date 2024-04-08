#' @rdname SingleCellAlleleExperiment
#'
#' @exportClass SingleCellAlleleExperiment SingleCellAlleleExperiment
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
.scae <- setClass("SingleCellAlleleExperiment",
                  contains = "SingleCellExperiment")
