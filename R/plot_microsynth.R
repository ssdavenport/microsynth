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
#' Only required input is a parameter \code{ms} which is a microsynth object.
#'
#' @param ms A microsynth object
#'
#' @param plot.var A vector of variable names giving the outcome variables that
#'   are shown in plots.  If \code{plot.var = NULL}, all outcome variables that
#'   are included in \code{ms} are plotted.  Only variables contained in
#'   the input \code{result.var} as used in the creation of \code{ms} can
#'   be plotted using \code{plot()}.
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
#' @param at A vector that gives the location of user-specified x-axis labels.
#'   \code{at} should be a (numeric) subset of the named time points contained
#'   in \code{ms} (e.g., \code{colnames(ms$Plot.Stats$Treatment)}).
#'
#' @param labels A vector of the same length as \code{at} that gives the names
#'   of the labels that will be marked at the times indicated by \code{at} in
#'   the plots.
#'
#' @param all A scalar character string giving the unit name for cases.
#'   If \code{NULL}, a third curve showing the overall outcome levels is
#'   not plotted.
#'
#' @param main.tc A scalar (or a vector of the same length as \code{plot.var})
#'   character string giving the title to be used for the first plots
#'   (that show treatment and control).  Defaults to \code{plot.var}.
#'
#' @param main.diff A scalar (or a vector of the same length as \code{plot.var})
#'   character string giving the title to be used for the second plots
#'   (that show differences between treatment and control).
#'   Defaults to \code{plot.var}.
#'
#' @param xlab.tc A scalar (or a vector of the same length as \code{plot.var})
#'   character string giving the x-axis labels to be used for the first plots
#'   (that show treatment and control).  Defaults to \code{''}.
#'
#' @param xlab.diff A scalar (or a vector of the same length as \code{plot.var})
#'   character string giving the x-axis labels to be used for the second plots
#'   (that show differences between treatment and control).
#'   Defaults to \code{''}.
#'
#' @param ylab.tc A scalar (or a vector of the same length as \code{plot.var})
#'   character string giving the y-axis labels to be used for the first plots
#'   (that show treatment and control).  Defaults to \code{plot.var}.
#'
#' @param ylab.diff A scalar (or a vector of the same length as \code{plot.var})
#'   character string giving the y-axis labels to be used for the second plots
#'   (that show differences between treatment and control).
#'   Defaults to \code{'Treatment - Control'}.
#'
#' @examples
#'
#' # Declare time-variant (outcome) and time-invariant variables for matching
#' cov.var <- c('TotalPop', 'BLACK', 'HISPANIC', 'Males_1521',
#'        'HOUSEHOLDS', 'FAMILYHOUS', 'FEMALE_HOU', 'RENTER_HOU', 'VACANT_HOU')
#'
#' match.out <- c('i_felony', 'i_misdemea', 'i_drugs', 'any_crime')
#'
#' set.seed(99199) # for reproducibility
#'
#' \donttest{
#'
#' # Perform matching and estimation, without permutations or jackknife
#' # runtime: <1 min
#' sea1 <- microsynth(seattledmi,
#'                   idvar="ID", timevar="time", intvar="Intervention",
#'                   start.pre=1, end.pre=12, end.post=16,
#'                   match.out=match.out, match.covar=cov.var,
#'                   result.var=match.out, omnibus.var=match.out,
#'                   test="lower",
#'                   n.cores = min(parallel::detectCores(), 2))
#'
#' # Plot with default settings in the GUI.
#' plot_microsynth(sea1)
#'
#' # Make plots, display, and save to a single file (plots.pdf).
#' plot_microsynth(sea1, file = file.path(tempdir(), 'plots.pdf'), sep = FALSE)
#'
#' # Make plots for only one outcome, display, and save to a single file.
#' plot_microsynth(sea1, plot.var = "any_crime",
#'      file = file.path(tempdir(), 'plots.pdf'), sep = FALSE)
#' }
#'
#' @export

plot_microsynth <- function (ms,
                             plot.var = NULL, start.pre = NULL, end.pre = NULL, end.post = NULL,
                             file=NULL, sep = TRUE, plot.first = NULL, legend.spot = "bottomleft",
                             height = NULL, width = NULL, at = NULL, labels = NULL,
                             all = "cases", main.tc = NULL, main.diff = NULL, xlab.tc = NULL,
                             xlab.diff = NULL, ylab.tc = NULL, ylab.diff = NULL) {

  if(!is.element("Plot.Stats",names(ms))) {
    stop("object ms does not contain output regarding results (e.g., only weights were generated).")
  }

  plotdat.t <- ms$Plot.Stats$Treatment
  plotdat.c <- ms$Plot.Stats$Control
  if(length(all) > 0) {
    plotdat.a <- ms$Plot.Stats$All
  } else {
    plotdat.a <- ms$Plot.Stats$Treatment
  }
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
    if (is.element(as.character(start.pre), time.names)) {
      start.pre <- match(as.character(start.pre), time.names)
    } else {
      start.pre <- which(as.numeric(time.names) >= as.numeric(start.pre))[1]
    }
  } else {
    start.pre <- 1
  }
  if (length(end.pre) > 0) {
    if (is.element(as.character(end.pre), time.names)) {
      end.pre <- match(as.character(end.pre), time.names)
    } else {
      end.pre <- which(as.numeric(time.names) <= as.numeric(end.pre))
      end.pre <- end.pre[length(end.pre)]
    }
  } else {
    end.pre <- ms$Plot.Stats$end.pre
    end.pre <- match(as.character(end.pre), time.names)
  }
  if (length(end.post) > 0) {
    if (is.element(as.character(end.post), time.names)) {
      end.post <- match(as.character(end.post), time.names)
    } else {
      end.post <- which(as.numeric(time.names) <= as.numeric(end.post))
      end.post <- end.post[length(end.post)]
    }
  } else {
    end.post <- length(time.names)
  }

  if (start.pre == end.post) {
    stop("Do not plot:  start.pre == end.post")
  }

  all.vars <- rownames(plotdat.t)

  if (length(plot.var) == 0) {
    plot.var <- all.vars
  }

  plot.var <- remove.vars(plot.var, all.vars, "plot.var")

  if(length(xlab.tc) == 0) {
    xlab.tc <- rep("", length(plot.var))
  } else if(length(xlab.tc) == 1) {
    xlab.tc <- rep(xlab.tc, length(plot.var))
  }

  if(length(xlab.diff) == 0) {
    xlab.diff <- rep("", length(plot.var))
  } else if(length(xlab.diff) == 1) {
    xlab.diff <- rep(xlab.diff, length(plot.var))
  }

  if(length(ylab.tc) == 0) {
    ylab.tc <- plot.var
  } else if(length(ylab.tc) == 1) {
    ylab.tc <- rep(ylab.tc, length(plot.var))
  }

  if(length(ylab.diff) == 0) {
    ylab.diff <- rep("Treatment - Control", length(plot.var))
  } else if(length(ylab.diff) == 1) {
    ylab.diff <- rep(ylab.diff, length(plot.var))
  }

  if(length(main.tc) == 0) {
    main.tc <- plot.var
  } else if(length(main.tc) == 1) {
    main.tc <- rep(main.tc, length(plot.var))
  }

  if(length(main.diff) == 0) {
    main.diff <- plot.var
  } else if(length(main.diff) == 1) {
    main.diff <- rep(main.diff, length(plot.var))
  }

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
    file.type <- "pdf"
    if (substr(file, nchar(file) - 3, nchar(file)) ==
        ".pdf") {
      file <- substr(file, 1, nchar(file) - 4)
    } else if (substr(file, nchar(file) - 3, nchar(file)) ==
               ".png") {
      file <- substr(file, 1, nchar(file) - 4)
      file.type <- "png"
    }
    if (!sep) {
      if (file.type == "pdf") {
        file <- paste(file, ".pdf", sep = "")
        grDevices::pdf(file = file, width = width, height = height)
      } else if (file.type == "png") {
        file <- paste(file, ".png", sep = "")
        grDevices::png(file = file, width = width, height = height,
                       units = "in", res = 500)
      }
      graphics::par(mfrow = c(3, 2), ask = FALSE)
    }
  }

  for (j in 1:length(plot.var)) {
    #use.mu <- which(mu[, plot.var[j]] > 0)
    use.mu <- 1:dim(plotdat.d)[2]
    for (i in 1:dim(plotdat.d)[2]) {
      tmp <- plotdat.d[plot.var[j], i, ]
      if (i == 1) {
        if (sep & length(file) > 0) {
          if (file.type == "pdf") {
            grDevices::pdf(file = paste(file, "_", plot.var[j],
                                        "_TC.pdf", sep = ""), width = width, height = height)
          } else if (file.type == "png") {
            grDevices::png(file = paste(file, "_", plot.var[j],
                                        "_TC.png", sep = ""), width = width, height = height,
                           units = "in", res = 500)
          }
        }
        tmp1 <- plotdat.t[plot.var[j], ]
        tmp2 <- plotdat.c[plot.var[j], ]
        tmp3 <- scale.by * plotdat.a[plot.var[j], ]
        lty <- c(1, 2, 4)
        col <- c(2, 1, 3)
        lwd <- c(2, 2, 2)
        ylim <- c(min(tmp1[tuse], tmp2[tuse]), max(tmp1[tuse], tmp2[tuse]))
        ylim <- c(min(ylim, tmp3[tuse]), max(ylim, tmp3[tuse]))
        ylim[2] <- 1.2 * ylim[2]
        xxnams1 <- as.numeric(as.character(time.names[xnams]))
        iend.pre <- as.numeric(as.character(time.names[end.pre]))
        if (sum(is.na(xxnams1)) > 0) {
          xxnams1 <- xnams
          iend.pre <- end.pre
        }
        if(length(at) == 0) {
          #at <- xxnams1
          xaxt <- "s"
        } else {
          xaxt <- "n"
        }
        #if(length(labels) == 0) {
        #  labels <- at
        #}
        xxnams2 <- xxnams1[tuse & use]
        xxnams3 <- xxnams1[tuse & nuse]
        xxnams1 <- xxnams1[tuse]
        xlim <- c(min(xxnams1), max(xxnams1))
        graphics::plot(xxnams1, tmp1[tuse], type = "l",
                       lty = lty[1], col = col[1], lwd = lwd[1],
                       xlim = xlim, xlab = xlab.tc[j], ylab = ylab.tc[j], main = main.tc[j],
                       ylim = ylim, xaxt = xaxt)
        if(length(at) > 0) {
          at1 <- intersect(at, xxnams1)
          labels1 <- labels[is.element(at, xxnams1)]
          graphics::axis(1, at=at1,labels=labels1)
        }
        graphics::abline(v = iend.pre, lty = 2, col = 2)
        graphics::lines(xxnams1, tmp2[tuse], type = "l",
                        lty = lty[2], col = col[2], lwd = lwd[2])
        if(length(all) > 0) {
          graphics::lines(xxnams1, tmp3[tuse], type = "l",
                          lty = lty[3], col = col[3], lwd = lwd[3])
          leg3 <- paste("All ",all," (scaled)", sep = "")
          leg <- c("Treatment", "Synthetic Control", leg3)
        } else {
          leg <- c("Treatment", "Synthetic Control")
          col <- col[1:2]
          lty <- lty[1:2]
          lwd <- lwd[1:2]
        }
        graphics::legend(legend.spot, legend = leg,
                         col = col, lty = lty, lwd = lwd, cex = 0.8,
                         bty = "n")
        if (sep & length(file) > 0) {
          grDevices::dev.off()
        }
        if (is.element(i, use.mu)) {
          if (sep & length(file) > 0) {
            if (file.type == "pdf") {
              grDevices::pdf(file = paste(file, "_",
                                          plot.var[j], "_Diff.pdf", sep = ""), width = width,
                             height = height)
            } else if (file.type == "png") {
              grDevices::png(file = paste(file, "_",
                                          plot.var[j], "_Diff.png", sep = ""), width = width,
                             height = height, units = "in", res = 500)
            }
          }
          ylim1 <- c(min(tmp[tuse]), max(tmp[tuse]))
          bigtmp <- plotdat.d[plot.var[j], use.mu, ]
          bigtmp <- matrix(bigtmp, length(use.mu), dim(plotdat.d)[3])
          ylim2a <- 2 * min(apply(bigtmp[ ,tuse, drop = FALSE], 2, stats::quantile, probs = .05, na.rm = TRUE))
          ylim2b <- 2 * max(apply(bigtmp[ ,tuse, drop = FALSE], 2, stats::quantile, probs = .95, na.rm = TRUE))
          ylim <- c(min(ylim1[1], ylim2a, na.rm = TRUE),
                    max(ylim1[2], ylim2b, na.rm = TRUE))
          xlim <- c(min(xxnams1), max(xxnams1))
          graphics::plot(xxnams2, tmp[tuse &
                                        use], type = "l", ylim = ylim, xlim = xlim,
                         col = 2, lty = 2, main = main.diff[j], xlab = xlab.diff[j],
                         ylab = ylab.diff[j], xaxt = xaxt)
          if(length(at) > 0) {
            at1 <- intersect(at, xxnams1)
            labels1 <- labels[is.element(at, xxnams1)]
            graphics::axis(1, at=at1,labels=labels1)
          }
        }
        else {
          if (!sep & length(file) > 0) {
            graphics::plot(1, 1, main = main.diff[j])
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




