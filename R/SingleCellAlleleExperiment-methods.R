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
  n_a_feats <- nrow(scae_subset(object, "alleles"))
  n_ag_feats <- nrow(get_agenes(object))
  n_fun_feats <- nrow(scae_subset_functional(object))
  n_immune_total <- n_a_feats + n_ag_feats + n_fun_feats

  callNextMethod()
  cat(
   "---------------\n",
   "Including a total of ", n_immune_total, " immune related features:\n",
   sep = ""
  )
  coolcat("Allele-level information (%d): %s\n",
          rownames(scae_subset_alleles(object)))
  coolcat("Immune genes (%d): %s\n",
          rownames(get_agenes(object)))
  coolcat("Functional level information (%d): %s\n",
          rownames(scae_subset_functional(object)))
}

#' @rdname SingleCellAlleleExperiment-misc
#'
#' @export
setMethod("show",
          signature = signature(object = "SingleCellAlleleExperiment"),
          definition = .scae_show)
