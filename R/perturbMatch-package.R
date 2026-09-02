#' perturbMatch: match transcriptomics signatures to single-gene perturbations
#'
#' @description
#' perturbMatch scores a query transcriptomics signature against a reference
#' database of single-gene perturbation signatures, ranking the perturbations
#' whose expression response most resembles the query and so nominating the
#' regulators most likely to drive the observed change.
#'
#' @section Getting started:
#' \describe{
#'   \item{\code{\link{getReferenceDatabase}}}{load a reference database, from
#'     \pkg{ExperimentHub} or from the bundled demonstration subset.}
#'   \item{\code{\link{prepareQuery}}}{turn a differential expression result
#'     into the ranked query \code{\link{calcEnrichment}} expects.}
#'   \item{\code{\link{calcEnrichment}}}{score the query against every
#'     reference signature, returning a \link[=PerturbMatch-class]{PerturbMatch}
#'     object.}
#'   \item{\code{\link{plotSimilarity}}, \code{\link{regulatorGSEA}}}{summarise
#'     and interpret the ranked regulators.}
#' }
#'
#' See \code{vignette("perturbMatch", package = "perturbMatch")} for a worked
#' example.
#'
#' @author Jamie Soul \email{jamie.soul@@liverpool.ac.uk}
#'
#' @references
#' Soul J, Young DA (2026). Automated generation of a gene perturbation
#' transcriptomic atlas using large language models. \emph{bioRxiv}.
#' \doi{10.64898/2026.08.08.743502}
#'
#' @keywords internal
"_PACKAGE"
