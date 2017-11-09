#' S3 class microsynth.
#' @exportClass microsynth
#' @export print.microsynth
#' @export summary.microsynth
makemicrosynth <- function(x) {
    class(x) <- c("microsynth")
    x
}

#' This convenience function displays output that compares treatment
#' characteristics
#' to synthetic control and the population, and shows estimated results in a
#' similar format as they appear when saved to .csv or .xlsx.
#'
print.microsynth <- function(microsynth_output) {
  out <- list(microsynth_output$w$Summary, microsynth_output$Results)
  names(out) <- c("Matching Summary", "Results")
  cat("Weight Balance Table: \n")
  print(out[[1]])
  cat("\nResults: \n")
  print(out[[2]])
}

#' This convenience function displays output that compares treatment
#' characteristics
#' to synthetic control and the population, and shows estimated results in a
#' similar format as they appear when saved to .csv or .xlsx.
#'
summary.microsynth <- function(microsynth_output) {
  out <- list(microsynth_output$w$Summary, microsynth_output$Results)
  names(out) <- c("Matching Summary", "Results")
  cat("Weight Balance Table: \n")
  print(out[[1]])
  cat("\nResults: \n")
  print(out[[2]])
}
