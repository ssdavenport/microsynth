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
#'
#' @param object A \code{microsynth} object produced by \code{microsynth()}
#' @param ... further arguments passed to or from other methods.
#'
#' @export
summary.microsynth <- function(object, ...) {

  cat("Weight Balance Table: \n")
  print(object$w$Summary, ...)

  # Be careful about case when $Results does not exist.
  if(!is.null(object$Results)) {
    cat("\nResults: \n")
    print.res(object$Results)
  }
}

#' Displaying microsynth Fits and Results
#'
#' Print method for class "microsynth".
#'
#' @param x A \code{microsynth} object produced by \code{microsynth()}
#' @param ... further arguments passed to or from other methods.
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
print.microsynth <- function(x, ...) {

  cat("Weight Balance Table: \n")
  print(x$w$Summary, ...)

  # Be careful about case when $Results does not exist.
  if(!is.null(x$Results)) {
    cat("\nResults: \n")
    print.res(x$Results)
  }
}

# This is a helper function for print/summary
print.res <- function(results.tables, ...) {
  # First, extract the end.post value (otherwise it will interfere with colnames)
  # Allow for an iterated process if there are 2 end times (and 2 tables)


  end.post <- names(results.tables)
  names(results.tables) <- NULL
  for(j in 1:length(end.post)) { # Iterate for each table (if end.post is a vector length >1)

    # Format the table, column by column
    res <- data.frame(results.tables[[j]])
    for(i in 1:NCOL(res)) {
      nam <- colnames(res)[i]
      if(nam=="Trt"){tmp <- round(res[,i],2)}
      else if(nam=="Con")
      {
        tmp <- formatC(round(res[,i],2), digits = 2, format = "f", flag = "0")
      }
      else if(grepl("pVal",nam))
      {
        tmp <- formatC(round(res[,i],4), digits = 4, format = "f", flag = "0")
      }
      else
      {
        tmp <- formatC(100*res[,i], digits = 1, format = "f", flag = "0")
        tmp <- paste(tmp,"%",sep="")
      }
      tmp[grepl("NA",tmp)]=NA
      tmp[is.na(tmp)]="--"
      res[,i] <- tmp
    }
    # Produce output (for each table)
    cat(paste0("\nend.post = ", end.post[j], "\n"))
    print(res, ...)
  }
}
