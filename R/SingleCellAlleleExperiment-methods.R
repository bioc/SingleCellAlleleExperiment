#' @name SingleCellAlleleExperiment-misc
#'
#' @title Miscellaneous SingleCellAlleleExperiment methods
#'
#' @description
#' Miscellaneous methods for the \code{\link{SingleCellAlleleExperiment}} class 
#' and its descendants that do not fit into any other documentation category 
#' such as, for example, show methods.
#'
#' @param object a \code{\link{SingleCellAlleleExperiment}} object
#'
#' @return Returns NULL
NULL


# SingleCellExperiment show method ----------------------------------------

#' @importFrom S4Vectors coolcat
#' @importFrom methods callNextMethod
#' @importMethodsFrom SingleCellExperiment show
.scae_show <- function(object) {
  n_tot_feats <- nrow(object)
  n_a_feats <- nrow(get_agenes(object))
  n_ni_feats <- nrow(get_nigenes(object))
  n_fun_feats <- nrow(scae_subset_functional(object))
  
  callNextMethod()
  cat(
   "----------\n", 
   "Including a total of ", n_tot_feats, " features:\n",
   sep = ""
  )
  coolcat("Immune genes (%d): %s\n", 
          rownames(get_agenes(object)))
  coolcat("Allele-level information (%d): %s\n", 
          rownames(scae_subset_alleles(object)))
  coolcat("Functional level information (%d): %s\n", 
          rownames(scae_subset_functional(object)))
  coolcat("Non-immune genes (%d): %s\n", 
          rownames(get_nigenes(object)))
}

#' @rdname SingleCellAlleleExperiment-misc
#' 
#' @export
setMethod("show",
          signature = signature(object = "SingleCellAlleleExperiment"),
          definition = .scae_show)

