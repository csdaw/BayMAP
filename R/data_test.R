#' A small PAR-CLIP data set as an example
#'
#' @description A small subsample on chromosome 1 of a publically available data
#'   set of \cite{Kishore et al. (2011)}.
#' @details
#' \describe{
#'   \item{\code{chr}}{a factor with the chromosomes as levels \code{1}, \dots,  \code{22}, \code{X}, \code{Y}}
#'   \item{\code{pos}}{a positive integer with the genomic position}
#'   \item{\code{mutation}}{a factor with substitution types as levels, e.g., \code{TC}}
#'   \item{\code{count}}{a positive integer with the number of observed substitutions}
#'   \item{\code{coverage}}{a positive integer with the number of observed reads for this position}
#' }
#' 
#' @author Eva-Maria Huessler, \email{eva-maria.huessler@uni-duesseldorf.de}
#'
#' @references \cite{Kishore et al. (2011)}, (SRA accession numbers: SRR189784)
#'
#' @keywords datasets
"data_test"