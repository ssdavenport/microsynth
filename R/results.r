
# Main function used to produce basic summary statistics; called by microsynth()
get.stats <- function(bigdat, w, inter, keep.groups, result.var = dimnames(bigdat)[[2]], 
    end.pre, period = 1, end.post = 80, file = NULL, sep = TRUE, start.pre = 25, legend.spot = "bottomleft", 
    omnibus.var = result.var, cut.mse = 1, scale.var = "Intercept", twosided = FALSE, time.names = NULL) {
    if (length(time.names) == 0) {
        time.names <- as.character(1:dim(bigdat)[3])
    }
    
    use.omnibus <- length(omnibus.var) > 0
    all.var <- union(result.var, omnibus.var)
    plot.it <- all.var
    stat5 <- stat4 <- stat2 <- stat1 <- mu <- matrix(NA, NCOL(w), length(result.var) + sum(use.omnibus))
    rownames(stat5) <- rownames(stat4) <- rownames(stat2) <- rownames(stat1) <- rownames(mu) <- colnames(w)
    if (use.omnibus) {
        colnames(stat5) <- colnames(stat4) <- colnames(stat2) <- colnames(stat1) <- colnames(mu) <- c(result.var, 
            "Omnibus")
    } else {
        colnames(stat5) <- colnames(stat4) <- colnames(stat2) <- colnames(stat1) <- colnames(mu) <- c(result.var)
    }
    bigdat1 <- make.quarter3(bigdat, period = period, end.pre = end.pre)
    keep <- keep.groups
    keep[1] <- TRUE
    synth <- list()
    inter <- inter[, !grepl("Jack", colnames(w)), drop = FALSE]
    keep <- keep[!grepl("Jack", colnames(w))]
    w <- w[, !grepl("Jack", colnames(w)), drop = FALSE]
    for (i in 1:NCOL(w)) {
        Inter <- inter[, i]
        use.con <- !is.na(Inter) & Inter == FALSE
        use.tre <- !is.na(Inter) & Inter == TRUE
        w.tmp <- w[use.con, i]
        condat1 <- bigdat1[use.con, all.var, , drop = FALSE]
        intdat1 <- bigdat1[use.tre, all.var, , drop = FALSE]
        test2 <- apply(w.tmp * condat1, c(2, 3), sum)
        test1 <- apply(intdat1, c(2, 3), sum)
        if (i == 1) {
            if (scale.var == "Intercept") {
                scale.by <- sum(use.tre)/sum(use.tre | use.con)
            } else {
                scale.by <- sum(bigdat1[use.tre, scale.var, 1])/sum(bigdat1[use.tre | use.con, 
                  scale.var, 1])
            }
            alldat1 <- bigdat1[use.tre | use.con, all.var, , drop = FALSE]
            test3 <- apply(alldat1, c(2, 3), sum)
        }
        if (i == 1) {
            xnams <- as.numeric(colnames(test1))
            tuse <- xnams <= end.post & xnams >= start.pre
            use <- xnams <= end.pre
            nuse <- xnams >= end.pre
            if (end.post > end.pre) {
                fuse <- !use & xnams <= end.post
            } else if (end.post == end.pre) {
                fuse <- xnams >= end.pre & xnams <= end.post
            } else {
                stop("end.post is less than end.pre")
            }
            if (length(plot.it) > 0) {
                no.jack <- which(!grepl("Jack", colnames(w)))
                keep1 <- keep[no.jack]
                plotdat.d <- array(NA, c(length(plot.it), length(no.jack), length(xnams)))
                dimnames(plotdat.d) <- list(plot.it, colnames(w)[no.jack], time.names[xnams])
                plotdat.a <- plotdat.t <- plotdat.c <- array(NA, c(length(plot.it), length(xnams)))
                dimnames(plotdat.a) <- dimnames(plotdat.t) <- dimnames(plotdat.c) <- list(plot.it, 
                  time.names[xnams])
            } else {
                plotdat.t <- plotdat.c <- plotdat.a <- plotdat.d <- NULL
            }
        }
        if (use.omnibus) {
            use.cols <- 1:(NCOL(stat1) - 1)
        } else {
            use.cols <- 1:NCOL(stat1)
        }
        mu[i, use.cols] <- sum(fuse) * rowMeans(test1[result.var, tuse, drop = FALSE])
        stat4[i, use.cols] <- rowSums(test1[result.var, fuse, drop = FALSE])
        stat5[i, use.cols] <- rowSums(test2[result.var, fuse, drop = FALSE])
        stat1[i, use.cols] <- stat4[i, use.cols, drop = FALSE] - stat5[i, use.cols, drop = FALSE]
        stat2[i, use.cols] <- stat1[i, use.cols, drop = FALSE]/rowSums(test2[result.var, 
            fuse, drop = FALSE])
        if (length(plot.it) > 0) {
            if (i == 1) {
                plotdat.t[plot.it, ] <- test1[plot.it, ]
                plotdat.c[plot.it, ] <- test2[plot.it, ]
                plotdat.a[plot.it, ] <- test3[plot.it, ]
                i1 <- i
            }
            if (is.element(i, no.jack)) {
                plotdat.d[plot.it, i1, ] <- test1[plot.it, ] - test2[plot.it, ]
                i1 <- i1 + 1
            }
        }
        if (use.omnibus) {
            if (!twosided) {
                stat1[i, NCOL(stat1)] <- sum(stat1[i, omnibus.var])
                stat2[i, NCOL(stat2)] <- sum(stat1[i, omnibus.var])/sum(mu[i, omnibus.var])
            } else {
                stat1[i, NCOL(stat1)] <- sum((stat1[i, omnibus.var])^2)
                stat2[i, NCOL(stat2)] <- sum((stat1[i, omnibus.var]/mu[i, omnibus.var])^2)
            }
        }
    }
    stats <- list(stat1[keep, , drop = FALSE], stat2[keep, , drop = FALSE], stat4[keep, 
        , drop = FALSE], stat5[keep, , drop = FALSE], list(Treatment = plotdat.t, Control = plotdat.c, 
        All = plotdat.a, Difference = plotdat.d, end.pre = end.pre, scale.by = scale.by))
    return(stats)
}


# Main function used to produce complex (survey) statistics; called by microsynth()
get.stats1 <- function(bigdat, w, inter, keep.groups, all.var, end.pre, period = 1, end.post = 80, 
    omnibus.var = NULL, G = 25, twosided = FALSE, printFlag = TRUE, n.cores = 1) {
    use.omnibus <- length(omnibus.var) > 0
    dof <- NA
    
    jack <- sum(grepl("Jack", colnames(w)))
    use.jack <- as.numeric(jack > 0)
    keep.groups <- keep.groups[!grepl("Jack", colnames(w))]
    w.jack <- w[, grepl("Jack", colnames(w))]
    w <- w[, !grepl("Jack", colnames(w)), drop = FALSE]
    inter.jack <- inter[, grepl("Jack", colnames(inter))]
    inter <- inter[, !grepl("Jack", colnames(inter)), drop = FALSE]
    
    if (length(G) == 0) {
        G <- min(table(inter[, 1]))
    }
    
    keep <- keep.groups
    
    if (sum(!keep) > 0) {
        w <- w[, keep]
        inter <- inter[, keep]
    }
    
    boot <- sum(grepl("Perm", colnames(w)))
    boot.nams <- colnames(w)[grepl("Perm", colnames(w))]
    
    if (use.jack == 1) {
        all.nams <- c(colnames(w)[1], "Jackknife", boot.nams)
    } else {
        all.nams <- colnames(w)
    }
    
    tot.col <- 1 + use.jack + boot
    boot.lower <- 2 + use.jack
    boot.upper <- 1 + use.jack + boot
    
    delta.out <- stat2 <- stat1 <- matrix(NA, tot.col, length(all.var) + sum(use.omnibus))
    rownames(delta.out) <- rownames(stat2) <- rownames(stat1) <- all.nams
    if (use.omnibus) {
        colnames(delta.out) <- colnames(stat2) <- colnames(stat1) <- c(all.var, "Omnibus")
    } else {
        colnames(delta.out) <- colnames(stat2) <- colnames(stat1) <- c(all.var)
    }
    
    for.max <- tot.col
    if (n.cores > 1 & tot.col > 1 + use.jack) {
        for.max <- 1 + use.jack
    }
    
    for (i in 1:for.max) {
        G.tmp <- G
        if (i == 1) {
            i.tmp <- i
        } else {
            i.tmp <- i - use.jack
        }
        is.jack <- use.jack == 1 & i == 2
        tmp <- proc.time()
        Inter <- inter[, i.tmp]
        use <- !is.na(Inter)
        w.tmp <- w[use, i.tmp]
        is.tre <- !is.na(Inter) & Inter == TRUE
        is.tre <- is.tre[use]
        
        if (i == 1) {
            test <- make.quarter2(bigdat[use, , , drop = FALSE], tre = is.tre, w = w.tmp, 
                period = period, end.pre = end.pre)
            if (end.post > end.pre) {
                use.test <- test[, 1] > end.pre & test[, 1] <= end.post
            } else if (end.post == end.pre) {
                use.test <- test[, 1] >= end.pre & test[, 1] <= end.post
            } else {
                stop("end.post is less than end.pre")
            }
            time <- test[use.test, 1, drop = FALSE]
            time.tmp <- length(table(time))
            time.tmp1 <- time.tmp
            treat <- test[use.test, 2, drop = FALSE]
            w.tmp <- test[use.test, 3, drop = FALSE]
            treat1 <- treat
            w.tmp1 <- w.tmp
            test <- test[use.test, -(1:3), drop = FALSE]
        } else if (is.jack) {
            w.tmp <- w.tmp1
            treat <- treat1
            remain <- end.pre%%period
            newcol <- (dim(bigdat)[3] - remain)%/%period
            time.tmp <- time.tmp1
            w.jack.tmp <- do.call("rbind", rep(list(w.jack), newcol))
            w.jack.tmp <- w.jack.tmp[use.test, , drop = FALSE]
            w.jack.tmp[is.na(w.jack.tmp)] <- 0
            G.tmp <- NCOL(w.jack.tmp)
        } else {
            remain <- end.pre%%period
            newcol <- (dim(bigdat)[3] - remain)%/%period
            time.tmp <- time.tmp1
            w.tmp <- rep(w.tmp, newcol)
            w.tmp <- w.tmp[use.test]
            treat <- rep(as.numeric(is.tre), newcol)
            treat <- treat[use.test]
        }
        
        for (j in 1:length(all.var)) {
            test.tmp <- data.frame(y = test[, all.var[j]], treat = as.numeric(treat), time = factor(time))
            time.tmp <- 1
            if (!is.jack) {
                design <- survey::svydesign(ids = ~0, data = test.tmp, weights = w.tmp)
            } else {
                sup.out <- suppressWarnings({
                  design <- survey::svrepdesign(data = test.tmp, repweights = w.jack.tmp, 
                    weights = w.tmp, combined.weights = TRUE, type = "JK1", mse = TRUE)
                })
            }
            usevars <- colnames(test.tmp)
            if (length(levels(factor(time))) == 1) {
                usevars <- setdiff(usevars, "time")
            }
            for (k in 2:length(usevars)) {
                if (k == 2) {
                  form = paste("y~", usevars[k], sep = "")
                } else {
                  form <- paste(form, "+", usevars[k], sep = "")
                }
            }
            if (!length(levels(factor(time))) == 1) {
                form <- paste(form, "-1", sep = "")
            } else {
                time.tmp <- 2
            }
            if (j == 1) {
                form1 <- substr(form, 2, 10000)
            }
            form <- stats::formula(form)
            sup.out <- suppressWarnings({
                mod <- survey::svyglm(form, design, family = stats::gaussian())
            })
            coefs <- summary(mod)$coefficients
            if (i == 1) {
                if (j == 1) {
                  out.coefs <- matrix(NA, length(all.var), NCOL(coefs))
                  rownames(out.coefs) <- all.var
                  colnames(out.coefs) <- colnames(coefs)
                }
                out.coefs[j, ] <- coefs[time.tmp, ]
            }
            if (sum(coefs[, "Std. Error"]) == Inf | sum(coefs[, "Std. Error"]) == -Inf) {
                is.inf <- TRUE
            } else {
                is.inf <- FALSE
            }
            # stat1[i, j] <- coefs[time.tmp, 't value']
            stat1[i, j] <- coefs[time.tmp, "Estimate"]/coefs[time.tmp, "Std. Error"]
            yC <- mean(coefs[-time.tmp, "Estimate"])
            yT <- yC + coefs[time.tmp, "Estimate"]
            if (yT/yC > 0) {
                stat2[i, j] <- log(yT/yC)
            } else {
                stat2[i, j] <- NA
            }
            delta.out[i, j] <- my.delta(mu = coefs[, "Estimate"], Sigma = stats::vcov(mod))
        }
        if (use.omnibus) {
            if (i == 1) {
                dum.tmp <- stat1[i, omnibus.var]
                keep <- !is.na(dum.tmp)
                if (sum(!keep) > 0) {
                  if (printFlag) {
                    message("\nThe following variables yield survey statistics with value NA. \nThese will be removed from the omnibus statistic: \n", 
                      appendLF = FALSE)
                  }
                  if (printFlag) {
                    message(paste0(utils::capture.output(omnibus.var[!keep]), collapse = "\n"))
                  }
                  if (printFlag) {
                    message("\n", appendLF = FALSE)
                  }
                }
                omnibus.var <- omnibus.var[keep]
            }
            test.tmp <- data.frame(test[, omnibus.var, drop = FALSE], treat = as.numeric(treat), 
                time = factor(time))
            form2 <- paste(omnibus.var, sep = "", collapse = ",")
            form2 <- paste("cbind(", form2, ")", sep = "")
            form1 <- paste(form2, form1, sep = "")
            form1 <- stats::formula(form1)
            if (i == 1) {
                reps <- assign.groups(as.numeric(treat), G = G.tmp)
            }
            coefs <- matrix(NA, length(omnibus.var), G.tmp)
            rownames(coefs) <- omnibus.var
            colnames(coefs) <- 1:G.tmp
            Sigma <- matrix(0, length(omnibus.var), length(omnibus.var))
            rownames(Sigma) <- colnames(Sigma) <- omnibus.var
            if (!is.jack) {
                for (g in 1:G.tmp) {
                  allmod <- stats::lm(form1, data = test.tmp, weights = w.tmp, subset = which(reps != 
                    g))
                  coefs[, g] <- as.matrix(stats::coef(allmod))["treat", ]
                }
            } else {
                for (g in 1:G.tmp) {
                  test.tmp.tmp <- test.tmp[!is.na(inter.jack[, g]), ]
                  w.jack.tmp.tmp <- w.jack.tmp[!is.na(inter.jack[, g]), g]
                  allmod <- stats::lm(form1, data = test.tmp.tmp, weights = w.jack.tmp.tmp)
                  coefs[, g] <- as.matrix(stats::coef(allmod))["treat", ]
                }
            }
            thetas <- rowMeans(coefs)
            dof <- length(thetas)
            coefs <- coefs - thetas
            thetas <- as.matrix(thetas)
            for (g in 1:G.tmp) {
                Sigma <- Sigma + tcrossprod(coefs[, g, drop = FALSE])
            }
            Sigma <- (G.tmp - 1)/G.tmp * Sigma
            if (i == 1) {
                rm.var <- find.sing(Sigma)
                keep.var <- 1:NCOL(Sigma)
                if (length(rm.var) > 0) {
                  keep.var <- keep.var[-rm.var]
                }
            }
            thetas <- thetas[keep.var, , drop = FALSE]
            Sigma <- Sigma[keep.var, keep.var, drop = FALSE]
            if (!twosided) {
                a <- as.matrix(diag(Sigma)^-0.5)
                stat1[i, NCOL(stat1)] <- crossprod(a, thetas)/sqrt(t(a) %*% Sigma %*% a)
            } else {
                stat1[i, NCOL(stat1)] <- crossprod(thetas, solve(Sigma)) %*% thetas
            }
        }
        tmp <- proc.time() - tmp
        if (i == 1) {
            if (printFlag) {
                message("Completed survey statistics for main weights: Time = ", round(tmp[3], 
                  2), "\n", sep = "", appendLF = FALSE)
            }
            if (jack == 0 & boot > 0) {
                if (printFlag) {
                  message("Calculating survey statistics for permutation groups...\n", appendLF = FALSE)
                }
                tmp.boot <- proc.time()
            }
            if (is.inf) {
                if (printFlag) {
                  message("WARNING: Infinite standard errors yielded by main weights.\n", 
                    appendLF = FALSE)
                }
            }
        } else if (is.jack) {
            if (printFlag) {
                message("Completed survey statistics for jackknife: Time = ", round(tmp[3], 
                  2), "\n", sep = "", appendLF = FALSE)
            }
            if (boot > 0) {
                if (printFlag) {
                  message("Calculating survey statistics for permutation groups...\n", appendLF = FALSE)
                }
                tmp.boot <- proc.time()
            }
            if (is.inf) {
                if (printFlag) {
                  message("WARNING: Infinite standard errors yielded by jackknife weights.\n", 
                    appendLF = FALSE)
                }
            }
        } else {
            if (i == boot.lower) {
                if (printFlag) {
                  message("Completed survey statistics for permutation group:\n", i - use.jack - 
                    1, sep = "", appendLF = FALSE)
                }
            } else if ((i - use.jack - 1)%%20 != 1 & i != boot.upper) {
                if (printFlag) {
                  message(", ", i - use.jack - 1, sep = "", appendLF = FALSE)
                }
            } else if ((i - use.jack - 1)%%20 == 1 & i != boot.upper) {
                if (printFlag) {
                  message(", \n", i - use.jack - 1, sep = "", appendLF = FALSE)
                }
            } else if ((i - use.jack - 1)%%20 != 1 & i == boot.upper) {
                if (printFlag) {
                  message(", ", i - use.jack - 1, "\n", sep = "", appendLF = FALSE)
                }
            } else {
                if (printFlag) {
                  message(", \n", i - use.jack - 1, "\n", sep = "", appendLF = FALSE)
                }
            }
            if (i == boot.upper) {
                if (i == boot.lower) {
                  if (printFlag) {
                    message("\n", appendLF = FALSE)
                  }
                }
                tmp.boot <- proc.time() - tmp.boot
                if (printFlag) {
                  message("Completed survey statistics for permutation groups: Time = ", 
                    round(tmp.boot[3], 2), "\n", sep = "", appendLF = FALSE)
                }
            }
            if (is.inf) {
                if (printFlag) {
                  message("WARNING: Infinite standard errors yielded by permutation group ", 
                    i, ".\n", sep = "", appendLF = FALSE)
                }
            }
        }
    }
    
    if (tot.col > for.max) {
        if (printFlag) {
            message("Parallelizing with n.cores = ", n.cores, "...\n", sep = "", appendLF = FALSE)
        }
        requireNamespace("parallel", quietly = TRUE)
        cl <- parallel::makeCluster(n.cores)
        
        list.out <- parallel::parLapply(cl = cl, X = (for.max + 1):tot.col, get.stats1.sub, 
            G, use.jack, boot.upper, boot.lower, inter, w, end.post, end.pre, period, bigdat, 
            all.var, use.omnibus, omnibus.var, twosided, test, use.test, time, time.tmp1, 
            reps, keep.var, twosided, printFlag = FALSE, tmp.boot, all.nams1 = colnames(delta.out))
        
        parallel::stopCluster(cl)
        
        tmp.boot <- proc.time() - tmp.boot
        if (printFlag) {
            message("Completed survey statistics for permutation groups: Time = ", round(tmp.boot[3], 
                2), "\n", sep = "", appendLF = FALSE)
        }
        
        for (index in (for.max + 1):tot.col) {
            delta.out[index, ] <- list.out[[index - for.max]][[1]]
            stat1[index, ] <- list.out[[index - for.max]][[2]]
            stat2[index, ] <- list.out[[index - for.max]][[3]]
        }
    }
    
    stat2[stat2 == -Inf | stat2 == Inf] <- NA
    return(list(stat1, stat2, delta.out, dof, out.coefs))
}



# Sub-function of get.stats1()
get.stats1.sub <- function(X, G, use.jack, boot.upper, boot.lower, inter, w, end.post, end.pre, 
    period, bigdat, all.var, use.omnibus, omnibus.var, two.sided, test, use.test, time, 
    time.tmp1, reps, keep.var, twosided, printFlag, tmp.boot, all.nams1) {
    
    i <- X
    G.tmp <- G
    if (i == 1) {
        i.tmp <- i
    } else {
        i.tmp <- i - use.jack
    }
    is.jack <- use.jack == 1 & i == 2
    tmp <- proc.time()
    Inter <- inter[, i.tmp]
    use <- !is.na(Inter)
    w.tmp <- w[use, i.tmp]
    is.tre <- !is.na(Inter) & Inter == TRUE
    is.tre <- is.tre[use]
    
    if (i == 1) {
        test <- make.quarter2(bigdat[use, , , drop = FALSE], tre = is.tre, w = w.tmp, period = period, 
            end.pre = end.pre)
        if (end.post > end.pre) {
            use.test <- test[, 1] > end.pre & test[, 1] <= end.post
        } else if (end.post == end.pre) {
            use.test <- test[, 1] >= end.pre & test[, 1] <= end.post
        } else {
            stop("end.post is less than end.pre")
        }
        time <- test[use.test, 1, drop = FALSE]
        time.tmp <- length(table(time))
        time.tmp1 <- time.tmp
        treat <- test[use.test, 2, drop = FALSE]
        w.tmp <- test[use.test, 3, drop = FALSE]
        treat1 <- treat
        w.tmp1 <- w.tmp
        test <- test[use.test, -(1:3), drop = FALSE]
    } else if (is.jack) {
        w.tmp <- w.tmp1
        treat <- treat1
        remain <- end.pre%%period
        newcol <- (dim(bigdat)[3] - remain)%/%period
        time.tmp <- time.tmp1
        w.jack.tmp <- do.call("rbind", rep(list(w.jack), newcol))
        w.jack.tmp <- w.jack.tmp[use.test, , drop = FALSE]
        w.jack.tmp[is.na(w.jack.tmp)] <- 0
        G.tmp <- NCOL(w.jack.tmp)
    } else {
        remain <- end.pre%%period
        newcol <- (dim(bigdat)[3] - remain)%/%period
        time.tmp <- time.tmp1
        w.tmp <- rep(w.tmp, newcol)
        w.tmp <- w.tmp[use.test]
        treat <- rep(as.numeric(is.tre), newcol)
        treat <- treat[use.test]
    }
    
    delta.out.out <- stat2.out <- stat1.out <- rep(NA, length(all.var) + sum(use.omnibus))
    names(delta.out.out) <- names(stat2.out) <- names(stat1.out) <- all.nams1
    for (j in 1:length(all.var)) {
        test.tmp <- data.frame(y = test[, all.var[j]], treat = as.numeric(treat), time = factor(time))
        time.tmp <- 1
        if (!is.jack) {
            design <- survey::svydesign(ids = ~0, data = test.tmp, weights = w.tmp)
        } else {
            sup.out <- suppressWarnings({
                design <- survey::svrepdesign(data = test.tmp, repweights = w.jack.tmp, 
                  weights = w.tmp, combined.weights = TRUE, type = "JK1", mse = TRUE)
            })
        }
        usevars <- colnames(test.tmp)
        if (length(levels(factor(time))) == 1) {
            usevars <- setdiff(usevars, "time")
        }
        for (k in 2:length(usevars)) {
            if (k == 2) {
                form = paste("y~", usevars[k], sep = "")
            } else {
                form <- paste(form, "+", usevars[k], sep = "")
            }
        }
        if (!length(levels(factor(time))) == 1) {
            form <- paste(form, "-1", sep = "")
        } else {
            time.tmp <- 2
        }
        if (j == 1) {
            form1 <- substr(form, 2, 10000)
        }
        form <- stats::formula(form)
        sup.out <- suppressWarnings({
            mod <- survey::svyglm(form, design, family = stats::gaussian())
        })
        coefs <- summary(mod)$coefficients
        if (i == 1) {
            if (j == 1) {
                out.coefs <- matrix(NA, length(all.var), NCOL(coefs))
                rownames(out.coefs) <- all.var
                colnames(out.coefs) <- colnames(coefs)
            }
            out.coefs[j, ] <- coefs[time.tmp, ]
        }
        if (sum(coefs[, "Std. Error"]) == Inf | sum(coefs[, "Std. Error"]) == -Inf) {
            is.inf <- TRUE
        } else {
            is.inf <- FALSE
        }
        # stat1[i, j] <- coefs[time.tmp, 't value']
        stat1.out[j] <- coefs[time.tmp, "Estimate"]/coefs[time.tmp, "Std. Error"]
        yC <- mean(coefs[-time.tmp, "Estimate"])
        yT <- yC + coefs[time.tmp, "Estimate"]
        if (yT/yC > 0) {
            stat2.out[j] <- log(yT/yC)
        } else {
            stat2.out[j] <- NA
        }
        delta.out.out[j] <- my.delta(mu = coefs[, "Estimate"], Sigma = stats::vcov(mod))
    }
    if (use.omnibus) {
        if (i == 1) {
            dum.tmp <- stat1[i, omnibus.var]
            keep <- !is.na(dum.tmp)
            if (sum(!keep) > 0) {
                if (printFlag) {
                  message("\nThe following variables yield survey statistics with value NA. \nThese will be removed from the omnibus statistic: \n", 
                    appendLF = FALSE)
                }
                if (printFlag) {
                  message(paste0(utils::capture.output(omnibus.var[!keep]), collapse = "\n"))
                }
                if (printFlag) {
                  message("\n", appendLF = FALSE)
                }
            }
            omnibus.var <- omnibus.var[keep]
        }
        test.tmp <- data.frame(test[, omnibus.var, drop = FALSE], treat = as.numeric(treat), 
            time = factor(time))
        form2 <- paste(omnibus.var, sep = "", collapse = ",")
        form2 <- paste("cbind(", form2, ")", sep = "")
        form1 <- paste(form2, form1, sep = "")
        form1 <- stats::formula(form1)
        if (i == 1) {
            reps <- assign.groups(as.numeric(treat), G = G.tmp)
        }
        coefs <- matrix(NA, length(omnibus.var), G.tmp)
        rownames(coefs) <- omnibus.var
        colnames(coefs) <- 1:G.tmp
        Sigma <- matrix(0, length(omnibus.var), length(omnibus.var))
        rownames(Sigma) <- colnames(Sigma) <- omnibus.var
        if (!is.jack) {
            for (g in 1:G.tmp) {
                allmod <- stats::lm(form1, data = test.tmp, weights = w.tmp, subset = which(reps != 
                  g))
                coefs[, g] <- as.matrix(stats::coef(allmod))["treat", ]
            }
        } else {
            for (g in 1:G.tmp) {
                test.tmp.tmp <- test.tmp[!is.na(inter.jack[, g]), ]
                w.jack.tmp.tmp <- w.jack.tmp[!is.na(inter.jack[, g]), g]
                allmod <- stats::lm(form1, data = test.tmp.tmp, weights = w.jack.tmp.tmp)
                coefs[, g] <- as.matrix(stats::coef(allmod))["treat", ]
            }
        }
        thetas <- rowMeans(coefs)
        dof <- length(thetas)
        coefs <- coefs - thetas
        thetas <- as.matrix(thetas)
        for (g in 1:G.tmp) {
            Sigma <- Sigma + tcrossprod(coefs[, g, drop = FALSE])
        }
        Sigma <- (G.tmp - 1)/G.tmp * Sigma
        if (i == 1) {
            rm.var <- find.sing(Sigma)
            keep.var <- 1:NCOL(Sigma)
            if (length(rm.var) > 0) {
                keep.var <- keep.var[-rm.var]
            }
        }
        thetas <- thetas[keep.var, , drop = FALSE]
        Sigma <- Sigma[keep.var, keep.var, drop = FALSE]
        if (!twosided) {
            a <- as.matrix(diag(Sigma)^-0.5)
            stat1.out[length(stat1.out)] <- crossprod(a, thetas)/sqrt(t(a) %*% Sigma %*% 
                a)
        } else {
            stat1.out[length(stat1.out)] <- crossprod(thetas, solve(Sigma)) %*% thetas
        }
    }
    tmp <- proc.time() - tmp
    if (i == 1) {
        if (printFlag) {
            message("Completed survey statistics for main weights: Time = ", round(tmp[3], 
                2), "\n", sep = "", appendLF = FALSE)
        }
        if (jack == 0 & boot > 0) {
            if (printFlag) {
                message("Calculating survey statistics for permutation groups...\n", appendLF = FALSE)
            }
            tmp.boot <- proc.time()
        }
        if (is.inf) {
            if (printFlag) {
                message("WARNING: Infinite standard errors yielded by main weights.\n", 
                  appendLF = FALSE)
            }
        }
    } else if (is.jack) {
        if (printFlag) {
            message("Completed survey statistics for jackknife: Time = ", round(tmp[3], 
                2), "\n", sep = "", appendLF = FALSE)
        }
        if (boot > 0) {
            if (printFlag) {
                message("Calculating survey statistics for permutation groups...\n", appendLF = FALSE)
            }
            tmp.boot <- proc.time()
        }
        if (is.inf) {
            if (printFlag) {
                message("WARNING: Infinite standard errors yielded by jackknife weights.\n", 
                  appendLF = FALSE)
            }
        }
    } else {
        if (i == boot.lower) {
            if (printFlag) {
                message("Completed survey statistics for permutation group:\n", i - use.jack - 
                  1, sep = "", appendLF = FALSE)
            }
        } else if ((i - use.jack - 1)%%20 != 1 & i != boot.upper) {
            if (printFlag) {
                message(", ", i - use.jack - 1, sep = "", appendLF = FALSE)
            }
        } else if ((i - use.jack - 1)%%20 == 1 & i != boot.upper) {
            if (printFlag) {
                message(", \n", i - use.jack - 1, sep = "", appendLF = FALSE)
            }
        } else if ((i - use.jack - 1)%%20 != 1 & i == boot.upper) {
            if (printFlag) {
                message(", ", i - use.jack - 1, "\n", sep = "", appendLF = FALSE)
            }
        } else {
            if (printFlag) {
                message(", \n", i - use.jack - 1, "\n", sep = "", appendLF = FALSE)
            }
        }
        if (i == boot.upper) {
            if (i == boot.lower) {
                if (printFlag) {
                  message("\n", appendLF = FALSE)
                }
            }
            tmp.boot <- proc.time() - tmp.boot
            if (printFlag) {
                message("Completed survey statistics for permutation groups: Time = ", round(tmp.boot[3], 
                  2), "\n", sep = "", appendLF = FALSE)
            }
        }
        if (is.inf) {
            if (printFlag) {
                message("WARNING: Infinite standard errors yielded by permutation group ", 
                  i, ".\n", sep = "", appendLF = FALSE)
            }
        }
    }
    return(list(delta.out = delta.out.out, stat1 = stat1.out, stat2 = stat2.out))
}



# Sub-function of get.stats(); reshape data
make.quarter3 <- function(dat, period = 1, end.pre) {
    n <- dim(dat)[1]
    p <- dim(dat)[2]
    q <- dim(dat)[3]
    add.back <- min(as.numeric(dimnames(dat)[[3]]))
    
    remain <- end.pre%%period
    newcol <- (q - remain)%/%period
    times <- period * (1:newcol) + remain
    times <- times + add.back - 1
    
    out <- array(NA, c(n, p, newcol))
    dimnames(out) <- list(dimnames(dat)[[1]], dimnames(dat)[[2]], times)
    start <- remain + 1
    for (i in 1:newcol) {
        stop <- start + period - 1
        out[, , i] <- apply(dat[, , start:stop, drop = FALSE], c(1, 2), sum)
        start <- stop + 1
    }
    return(out)
}



# Sub-function of get.stats1(), get.stats1.sub(); reshape data
make.quarter2 <- function(dat, tre, w, period = 1, end.pre) {
    n <- dim(dat)[1]
    p <- dim(dat)[2]
    q <- dim(dat)[3]
    add.back <- min(as.numeric(dimnames(dat)[[3]]))
    remain <- end.pre%%period
    newcol <- (q - remain)%/%period
    times <- period * (1:newcol) + remain
    times <- times + add.back - 1
    
    out <- matrix(NA, n * length(times), p + 3)
    colnames(out) <- c("Time", "Treatment", "w", dimnames(dat)[[2]])
    start <- remain + 1
    for (i in 1:length(times)) {
        stop <- start + period - 1
        tmp <- dat[, , start:stop, drop = FALSE]
        here <- ((i - 1) * n + 1):(i * n)
        out[here, 1] <- times[i]
        out[here, 2] <- as.numeric(tre)
        out[here, 3] <- w
        out[here, 4:NCOL(out)] <- apply(tmp, c(1, 2), sum)
        start <- stop + 1
    }
    
    return(out)
}



# Sub-function of get.stats1(), get.stats1.sub(); calculate confidence intervals
my.delta <- function(mu, Sigma) {
    inds1 <- which(rownames(Sigma) == "treat")
    inds2 <- which(rownames(Sigma) != "treat")
    inds <- c(inds1, inds2)
    mu <- mu[inds]
    Sigma <- Sigma[inds, inds]
    mu1 <- mu[1]
    mu2 <- mu[2:length(mu)]
    m <- length(mu2)
    MU2 <- mean(mu2)
    D <- rep(NA, length(mu))
    D[1] <- 1/(mu1 + MU2)
    D[2:length(D)] <- (D[1] - (1/MU2))/m
    D <- as.matrix(D)
    out <- crossprod(D, Sigma) %*% D
    return(c(out))
}
