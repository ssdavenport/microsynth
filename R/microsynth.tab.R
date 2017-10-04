#' Summarize microsynth Results
#'
#' This convenience function displays output that compares treatment
#' characteristics
#' to synthetic control and the population, and shows estimated results in a
#' similar format as they appear when saved to .csv or .xlsx.
#'
#' @param
#' microsynth.res the list resulting from a call to \code{\link{microsynth}}
#'
#' @return
#'
#' @export
microsynth.tab <- function(microsynth.res) {
    out <- list(microsynth.res$w$Summary, microsynth.res$Results)
    names(out) <- c("Matching Summary", "Results")
    out
}
