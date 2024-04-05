validSingleCellAlleleExperiment <- function(object) {
  msg <- NULL
  
  ## here: grow the message as we find things that don't look like as they should
  
  if (is.null(msg)) {
    TRUE
  } else msg
  
  
}

#' @importFrom S4Vectors setValidity2
S4Vectors::setValidity2("SingleCellAlleleExperiment", 
                        validSingleCellAlleleExperiment)

