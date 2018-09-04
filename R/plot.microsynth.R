#' @title
#' Plotting for microsynth objects.
#'
#' @description
#' Using a \code{microsynth} object as an input, this function gives time
#' series plots of selected outcomes.
#'
#' @details
#' Plots are given over both pre- and
#' intervention time periods and shown in terms of raw outcome values or
#' treatment/control differences.  Time series of permutation groups may be
#' overlaid to help illustrate statistical uncertainty.
#'
#' Only required input is a parameter \code{ms} which is a microsythn object.
#'
#' @param ms A microsynth object
#'
#' @param plot.var A vector of variable names giving the outcome variables that
#'   are shown in plots.  If \code{plot.var = NULL}, all outcome variables that
#'   are included in \code{ms} are plotted.  Only variables contained in
#'   the input \code{result.var} as used in the creation of \code{ms} can
#'   be plotted using \code{plot.microsynth()}.
#'
#' @param file A character string giving the name of file that will be
#'   created in the home directory containing plots.
#'  The name should have a \code{.pdf} extension.
#'
#' @param start.pre An integer indicating the time point that corresponds to the
#'   earliest time period that will be plotted.
#'   When \code{start.pre = NULL}, it is reset to the
#'   minimum time appearing in \code{ms}.
#'
#' @param end.pre An integer that gives the final time point of the
#'   pre-intervention period.  That is, \code{end.pre} is the last time at
#'   which treatment and synthetic control will were matched to one another.
#'   All time points
#'   following \code{end.pre} are considered to be post-intervention and the
#'   behavior of outcomes will be compared between the treatment and synthetic
#'   control groups across those time periods.
#'   If \code{end.pre = NULL} the end of the pre-intevention period will be
#'   determined from the object \code{ms}.
#'
#' @param end.post An integer that gives final time point that will be plottd.
#'   When \code{end.post = NULL} (the default), it is reset
#'   to the maximum time that appears in \code{ms}.
#'
#' @param sep If \code{sep = TRUE}, separate plots will be generated for each
#'   outcome.  Applicable only if plots are saved to file (
#'   \code{plot.file} is \code{non-NULL}). To change display of plots produced
#'   as output, use \code{\link[graphics]{par}}.
#'
#' @param legend.spot The location of the legend in the plots.
#'
#' @param plot.first The number of permutation groups to plot.
#'
#' @param height The height of the graphics region (in inches)
#'   when a pdf is created.
#'
#' @param width The width of the graphics region (in inches)
#'   when a pdf is created.
#'
#' @examples
#' set.seed(99199) # for reproducibility
#' # create a microsynth object
#' sea2 <- microsynth(seattledmi, idvar='ID', timevar='time',
#'         intvar='Intervention', start.pre=1, end.pre=12, end.post=c(16),
#'         match.out=match.out, match.covar=cov.var, result.var=match.out,
#'         omnibus.var=match.out, test='lower', perm=250,
#'         jack=TRUE, result.file='ExResults2.xlsx')
#'
#' # Create plots within the GUI.
#'
#' plot.microsynth(sea2)
#'
#' # Create plots and output to a single file (plots.pdf).
#'
#' plot.microsynth(sea2, file = "plots.pdf", sep = FALSE)
#'
#' # Create plots for only one outcome and output to separate files.
#'
#' plot.microsynth(sea2, plot.var = "any_crime", file = "plots.pdf", sep = FALSE)
#'
#' @export

plot.microsynth <- function (ms,
                             plot.var = NULL, start.pre = NULL, end.pre = NULL, end.post = NULL,
                             file=NULL, sep = TRUE, plot.first = NULL, legend.spot = "bottomleft",
                             height = NULL, width = NULL) {

plotdat.t <- ms$Plot.Stats$Treatment
plotdat.c <- ms$Plot.Stats$Control
plotdat.a <- ms$Plot.Stats$All
plotdat.d <- ms$Plot.Stats$Difference
scale.by <- ms$Plot.Stats$scale.by

if(length(height) == 0) {
  if (!sep) {
    height <- 11
  } else {
    height <- 5
  }
}

if(length(width) == 0) {
  if (!sep) {
    width <- 8.5
  } else {
    width <- 5
  }
}

if(length(plot.first) == 0) {
  plot.first <- dim(plotdat.d)[2] - 1
}

time.names <- colnames(plotdat.t)

if (length(start.pre) > 0) {
  start.pre <- match(as.character(start.pre), time.names)
} else {
  start.pre <- 1
}
if (length(end.pre) > 0) {
    end.pre <- match(as.character(end.pre), time.names)
} else {
  end.pre <- ms$Plot.Stats$end.pre
}
if (length(end.post) > 0) {
  end.post <- match(as.character(end.post), time.names)
} else {
  end.post <- length(time.names)
}

all.vars <- rownames(plotdat.t)

if (length(plot.var) == 0) {
  plot.var <- all.vars
}

plot.var <- remove.vars(plot.var, all.vars, "plot.var")

#xnams <- as.numeric(colnames(test1))
xnams <- 1:length(time.names)
tuse <- xnams <= end.post & xnams >= start.pre
use <- xnams <= end.pre
nuse <- xnams >= end.pre
if (end.post > end.pre) {
  fuse <- !use & xnams <= end.post
}
else if (end.post == end.pre) {
  fuse <- xnams >= end.pre & xnams <= end.post
}
else {
  stop("end.post is less than end.pre")
}

if (length(file) == 0) {
  graphics::par(mfrow = c(1, 2), ask = FALSE)
}
else {
  if (substr(file, nchar(file) - 3, nchar(file)) ==
      ".pdf") {
    file <- substr(file, 1, nchar(file) - 4)
  }
  if (!sep) {
    file <- paste(file, ".pdf", sep = "")
    grDevices::pdf(file = file, width = width, height = height)
    graphics::par(mfrow = c(3, 2), ask = FALSE)
  }
}

  for (j in 1:length(plot.var)) {
    ylab1 <- main1 <- main <- plot.var[j]
    ylab2 <- "Treatment - Control"
    #use.mu <- which(mu[, plot.var[j]] > 0)
    use.mu <- 1:dim(plotdat.d)[2]
    for (i in 1:dim(plotdat.d)[2]) {
      tmp <- plotdat.d[plot.var[j], i, ]
      if (i == 1) {
        if (sep & length(file) > 0) {
          grDevices::pdf(file = paste(file, "_", main1,
                                      "_TC.pdf", sep = ""), width = width, height = width)
        }
        tmp1 <- plotdat.t[plot.var[j], ]
        tmp2 <- plotdat.c[plot.var[j], ]
        tmp3 <- scale.by * plotdat.a[plot.var[j], ]
        lty <- c(1, 2, 4)
        col <- c(2, 1, 3)
        lwd <- c(2, 2, 2)
        ylim <- c(min(tmp1, tmp2), max(tmp1, tmp2))
        ylim <- c(min(ylim, tmp3), max(ylim, tmp3))
        ylim[2] <- 1.2 * ylim[2]
        xxnams1 <- as.numeric(as.character(time.names[xnams[tuse]]))
        xxnams2 <- as.numeric(as.character(time.names[xnams[tuse & use]]))
        xxnams3 <- as.numeric(as.character(time.names[xnams[tuse & nuse]]))
        iend.pre <- as.numeric(as.character(time.names[end.pre]))
        if (sum(is.na(xxnams1)) > 0) {
          xxnams1 <- xnams[tuse]
          xxnams2 <- xnams[tuse & use]
          xxnams3 <- xnams[tuse & nuse]
          iend.pre <- end.pre
        }
        xlim <- c(min(xxnams1), max(xxnams1))
        graphics::plot(xxnams1, tmp1[tuse], type = "l",
                       lty = lty[1], col = col[1], lwd = lwd[1],
                       xlim = xlim, xlab = "", ylab = ylab1, main = main,
                       ylim = ylim)
        graphics::abline(v = iend.pre, lty = 2, col = 2)
        graphics::lines(xxnams1, tmp2[tuse], type = "l",
                        lty = lty[2], col = col[2], lwd = lwd[2])
        graphics::lines(xxnams1, tmp3[tuse], type = "l",
                        lty = lty[3], col = col[3], lwd = lwd[3])
        leg <- c("Treatment", "Synthetic Control",
                 "All blocks (scaled)")
        graphics::legend(legend.spot, legend = leg,
                         col = col, lty = lty, lwd = lwd, cex = 0.8,
                         bty = "n")
        if (sep & length(file) > 0) {
          grDevices::dev.off()
        }
        if (is.element(i, use.mu)) {
          if (sep & length(file) > 0) {
            grDevices::pdf(file = paste(file, "_",
                                        main1, "_Diff.pdf", sep = ""), width = 5,
                           height = 5)
          }
          ylim1 <- c(min(tmp), max(tmp))
          bigtmp <- plotdat.d[plot.var[j], use.mu, ]
          ylim2 <- 2 * c(stats::quantile(bigtmp, 0.05,
                                         na.rm = TRUE), stats::quantile(bigtmp,
                                                                        0.95, na.rm = TRUE))
          ylim <- c(min(ylim1[1], ylim2[1], na.rm = TRUE),
                    max(ylim1[2], ylim2[2], na.rm = TRUE))
          xlim <- c(min(xxnams1), max(xxnams1))
          graphics::plot(xxnams2, tmp[tuse &
                                        use], type = "l", ylim = ylim, xlim = xlim,
                         col = 2, lty = 2, main = main, xlab = "",
                         ylab = ylab2)
        }
        else {
          if (!sep & length(file) > 0) {
            graphics::plot(1, 1, main = plot.var[j])
          }
        }
      }
      else {
        if (i <= min((dim(plotdat.d)[2]), plot.first)) {
          if (is.element(i, use.mu) & is.element(1,
                                                 use.mu)) {
            graphics::lines(xxnams1, tmp[tuse],
                            col = "azure3")
          }
        }
      }
      if (i == dim(plotdat.d)[2] & is.element(1, use.mu)) {
        tmp <- plotdat.d[plot.var[j], 1, ]
        graphics::lines(xxnams2, tmp[tuse &
                                       use], lty = 1, col = 1, lwd = 2)
        graphics::lines(xxnams3, tmp[tuse &
                                       nuse], col = 2, lwd = 2)
        graphics::abline(v = iend.pre, col = 2, lwd = 1,
                         lty = 2)
        graphics::abline(h = 0, lwd = 1, lty = 2)
        if (sep & length(file) > 0) {
          grDevices::dev.off()
        }
      }
    }
  }
  if (!sep & length(file) > 0) {
    grDevices::dev.off()
  }
}
