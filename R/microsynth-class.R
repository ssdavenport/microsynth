#' Summarizing microsynth Fits and Results
#'
#' Summary method for class 'microsynth'.
#'
#' @return The functions \code{print.microsynth} and
#' \code{summary.microsynth} displays information about the microsynth
#' fit and estimation results, if available.
#'
#' The output includes two parts: 1) a matching summary that compares
#' characteristics of the treatment to the synthetic control and the
#' population; and 2) estimated results, in a similar format as they
#' appear when saved to .csv or .xlsx., once for each specified
#' post-intervention evaluation time.
#'
#' @param object A \code{microsynth} object produced by \code{microsynth()}
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @examples
#'
#' # Use seattledmi, block-level panel data, to evaluate a crime intervention.
#'
#' # Declare time-variant (outcome) and time-invariant variables for matching
#' cov.var <- c('TotalPop', 'BLACK', 'HISPANIC', 'Males_1521',
#'        'HOUSEHOLDS', 'FAMILYHOUS', 'FEMALE_HOU', 'RENTER_HOU', 'VACANT_HOU')
#'
#' match.out <- c('i_felony', 'i_misdemea', 'i_drugs', 'any_crime')
#' set.seed(99199) # for reproducibility
#'
#' # Perform matching and estimation, without permutations or jackknife
#' # runtime: < 1 min
#'
#' \donttest{
#' sea1 <- microsynth(seattledmi,
#'                   idvar='ID', timevar='time', intvar='Intervention',
#'                   start.pre=1, end.pre=12, end.post=16,
#'                   match.out=match.out, match.covar=cov.var,
#'                   result.var=match.out, omnibus.var=match.out,
#'                   test='lower',
#'                   n.cores = min(parallel::detectCores(), 2))
#'
#' # View results
#' summary(sea1)
#' }
#'
#' @export
summary.microsynth <- function(object, ...) {

    cat("Weight Balance Table: \n\n")
    print(object$w$Summary, ...)

    # Be careful about case when $Results does not exist.
    if (!is.null(object$Results)) {
        cat("\nResults: \n")
        print_res(object$Results)
    }
}



#' Displaying microsynth Fits and Results
#'
#' Print method for class 'microsynth'.
#'
#' @param x A \code{microsynth} object produced by \code{microsynth()}
#' @param ... further arguments passed to or from other methods.
#'
#' @return The functions \code{print.microsynth} and
#' \code{summary.microsynth} display information about the microsynth
#' fit and estimation results, if available.
#'
#' The output includes two parts: 1) a display of key input parameters;
#' and 2) estimated results, in a similar format as they
#' appear when saved to .csv or .xlsx., once for each specified
#' post-intervention evaluation time.
#'
#' @examples
#'
#' # Use seattledmi, block-level panel data, to evaluate a crime intervention.
#'
#' # Declare time-variant (outcome) and time-invariant variables for matching
#' cov.var <- c('TotalPop', 'BLACK', 'HISPANIC', 'Males_1521',
#'        'HOUSEHOLDS', 'FAMILYHOUS', 'FEMALE_HOU', 'RENTER_HOU', 'VACANT_HOU')
#'
#' match.out <- c('i_felony', 'i_misdemea', 'i_drugs', 'any_crime')
#' set.seed(99199) # for reproducibility
#'
#' # Perform matching and estimation, without permutations or jackknife
#' # runtime: < 1 min
#'
#' \donttest{
#' sea1 <- microsynth(seattledmi,
#'                   idvar='ID', timevar='time', intvar='Intervention',
#'                   start.pre=1, end.pre=12, end.post=16,
#'                   match.out=match.out, match.covar=cov.var,
#'                   result.var=match.out, omnibus.var=match.out,
#'                   test='lower',
#'                   n.cores = min(parallel::detectCores(), 2))
#'
#' # View results
#' print(sea1)
#' }
#'
#' @export
print.microsynth <- function(x, ...) {

    cat(paste0("\tmicrosynth object\n\n"))

    cat("Scope:\n")
    cat(paste0("\tUnits:\t\t\tTotal: ", x$info$nUnits, "\tTreated: ", x$info$nTreatment, "\tUntreated: ", x$info$nControl, "\n"))
    cat(paste0("\tStudy Period(s):\tPre-period: ", x$info$start.pre, " - ", x$info$end.pre, "\tPost-period: ", x$info$end.pre + 1, " - ", x$info$end.post, "\n"))
    cat(paste0("\tConstraints:\t\tExact Match: ", x$info$num.constr[1], "\t\tMinimized Distance: ", x$info$num.constr[2]))
    cat(paste0("\nTime-variant outcomes:", "\n\tExact Match: ", paste0(x$info$match, collapse = ", "), " (", length(x$info$match), ")", "\n\tMinimized Distance: ", paste0(x$info$match.min, collapse = ", "), "(", length(x$info$match.min),
        ")", "\n"))
    cat(paste0("Time-invariant covariates:", "\n\tExact Match: ", paste0(x$info$covar, collapse = ", "), " (", length(x$info$covar), ")", "\n\tMinimized Distance: ", paste0(x$info$covar.min, collapse = ", "), "(",
        length(x$info$covar.min), ")"))

    if (!is.null(x$Results)) {
        cat("\n\nResults:")
        print_res(x$Results)
    }
}

# This is a helper function for print/summary
print_res <- function(results.tables, ...) {
    # First, extract the end.post value (otherwise it will interfere with colnames) Allow for an iterated process if there are 2 end times (and 2 tables)

    end.post <- names(results.tables)
    names(results.tables) <- NULL
    for (j in 1:length(end.post)) {
        # Iterate for each table (if end.post is a vector length >1)

        # Format the table, column by column
        res <- data.frame(results.tables[[j]])
        for (i in 1:NCOL(res)) {
            nam <- colnames(res)[i]
            if (nam == "Trt") {
                tmp <- round(res[, i], 2)
            } else if (nam == "Con") {
                tmp <- formatC(round(res[, i], 2), digits = 2, format = "f", flag = "0")
            } else if (grepl("pVal", nam)) {
                tmp <- formatC(round(res[, i], 4), digits = 4, format = "f", flag = "0")
            } else {
                tmp <- formatC(100 * res[, i], digits = 1, format = "f", flag = "0")
                tmp <- paste(tmp, "%", sep = "")
            }
            tmp[grepl("NA", tmp)] = NA
            tmp[is.na(tmp)] = "--"
            res[, i] <- tmp
        }
        # Produce output (for each table)
        cat(paste0("\nend.post = ", end.post[j], "\n"))
        print(res, ...)
    }
}
