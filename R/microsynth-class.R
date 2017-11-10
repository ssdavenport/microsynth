#' Summarizing microsynth Fits and Results
#'
#' Summary method for class "microsynth".
#'
#' @return The functions \code{print.microsynth} and
#' \code{summary.microsynth} display information about the microsynth
#' fit and estimation results, if available.
#'
#' The output includes two parts: 1) a matching summary that compares
#' characteristics of the treatment to the synthetic control and the
#' population; and 2) estimated results, in a similar format as they
#' appear when saved to .csv or .xlsx., once for each specified
#' post-intervention evaluation time.
#' @export
summary.microsynth <- function(microsynth_output) {
  out <- list(microsynth_output$w$Summary, microsynth_output$Results)
  names(out) <- c("Matching Summary", "Results")
  cat("Weight Balance Table: \n")
  print(out[[1]])
  cat("\nResults: \n")
  print(out[[2]])
}

#' Displaying microsynth Fits and Results
#'
#' Print method for class "microsynth".
#'
#' @return The functions \code{print.microsynth} and
#' \code{summary.microsynth} display information about the microsynth
#' fit and estimation results, if available.
#'
#' The output includes two parts: 1) a matching summary that compares
#' characteristics of the treatment to the synthetic control and the
#' population; and 2) estimated results, in a similar format as they
#' appear when saved to .csv or .xlsx., once for each specified
#' post-intervention evaluation time.
#' @export
print.microsynth <- function(microsynth_output) {
  out <- list(microsynth_output$w$Summary, microsynth_output$Results)
  names(out) <- c("Matching Summary", "Results")
  cat("Weight Balance Table: \n")
  print(out[[1]])
  cat("\nResults: \n")
  print(out[[2]])
}
