
get.w <- function (bigdat, covar.var, covar.var1 = NULL, dum, dum1 = NULL,
                   boot = 0, jack = 0, Int, int.val = 1, trim = NULL, maxit = 500,
                   cal.epsilon = 1e-04, end.pre, bounds = c(-Inf, Inf), calfun = "raking",
                   qpmeth = "LowRankQP", check.feas = FALSE, use.backup = TRUE,
                   scale.var = "Intercept", cut.mse = 1, time.names = NULL,
                   keep.int = FALSE, printFlag = TRUE, n.cores = 1)
{
  # Primary function used for calculating weights
  n <- dim(bigdat)[1]
  n.int <- sum(Int == int.val)
  back.state <- back.state1 <- back.state2 <- back.state3 <- ""
  fin.boots <- FALSE
  boots <- rep.G <- tmp.jack <- NULL
  if (boot > 0) {
    n.choose <- choose(n, n.int)
    if (boot > n.choose - 1) {
      boot <- n.choose - 1
      if(printFlag){message("Resetting perm = ", boot, "\n", sep = "",
                            appendLF = FALSE)}
    }
    if (n.choose - 1 <= max(1e+06, boot)) {
      boots <- utils::combn(1:n, n.int)
      check.combn <- function(x, y) {
        return(sum(!is.element(x, y)))
      }
      is.trt.area <- which(apply(boots, 2, check.combn,
                                 x = which(Int == int.val)) == 0)
      boots <- boots[, base::sample((1:NCOL(boots))[-is.trt.area],
                                    boot), drop = FALSE]
      fin.boots <- TRUE
    }
  }
  use.model <- 1
  if (check.feas & use.backup) {
    tmp <- proc.time()
    if(printFlag){message("Checking feasibility of first model...", "\n",
                          sep = "", appendLF = FALSE)}
    is.sol <- is.feasible(bigdat, covar.var, dum, Int = Int,
                          int.val = int.val, end.pre = end.pre, eps = 1e-04)
    tmp <- proc.time() - tmp
    if (!is.sol) {
      use.model <- 2
      if(printFlag){message("First model is infeasible: Time = ", round(tmp[3],
                                                                        2), "\n", sep = "", appendLF = FALSE)}
      dum.tmp <- merge.dums(dum, dum1)
      tmp <- proc.time()
      if(printFlag){message("Checking feasibility of second model...",
                            "\n", sep = "", appendLF = FALSE)}
      is.sol <- is.feasible(bigdat, covar.var, dum.tmp[[1]],
                            Int = Int, int.val = int.val, end.pre = end.pre,
                            eps = 1e-04)
      tmp <- proc.time() - tmp
      if (!is.sol) {
        use.model <- 3
        if(printFlag){message("Second model is infeasible: Time = ",
                              round(tmp[3], 2), "\n", sep = "", appendLF = FALSE)}
        if(printFlag){message("Will use third model.", "\n", sep = "",
                              appendLF = FALSE)}
      }
      else {
        if(printFlag){message("Second model is feasible: Time = ",
                              round(tmp[3], 2), "\n", sep = "", appendLF = FALSE)}
      }
    }
    else {
      if(printFlag){message("First model is feasible: Time = ", round(tmp[3],
                                                                      2), "\n", sep = "", appendLF = FALSE)}
    }
  }
  newdat <- get.newdat(bigdat, dum = dum, dum1 = dum1, covar.var = covar.var,
                       covar.var1 = covar.var1, end.pre = end.pre, time.names = time.names)
  newdat1 <- newdat[[2]]
  newdat <- newdat[[1]]
  duma <- merge.dums(dum, dum1)
  newdata <- get.newdat(bigdat, dum = duma[[1]], dum1 = duma[[2]], covar.var = covar.var,
                        covar.var1 = covar.var1, end.pre = end.pre, time.names = time.names)
  newdat1a <- newdata[[2]]
  newdata <- newdata[[1]]
  dumb <- merge.dums(duma[[1]], duma[[2]])
  covar.var1b <- union(covar.var, covar.var1)
  covar.varb <- NULL
  newdatb <- get.newdat(bigdat, dum = dumb[[1]], dum1 = dumb[[2]], covar.var = covar.varb,
                        covar.var1 = covar.var1b, end.pre = end.pre, time.names = time.names)
  newdat1b <- newdatb[[2]]
  newdatb <- newdatb[[1]]
  colnam <- "Main"
  if (is.logical(jack)) {
    if (jack) {
      jack <- min(table(Int == int.val))
    }
    else {
      jack <- 0
    }
  }
  if (jack > min(table(Int == int.val))) {
    jack <- min(table(Int == int.val))
    if(printFlag){message("Resetting jack = ", jack, "\n", sep = "", appendLF = FALSE)}
  }
  if (jack > 0) {
    rep.G <- assign.groups(Int == int.val, G = jack)
    colnam <- c(colnam, paste("Jack", 1:jack, sep = ""))
  }
  if (boot > 0) {
    colnam <- c(colnam, paste("Perm", 1:boot, sep = ""))
  }
  wghts <- Inter <- matrix(NA, n, boot + jack + 1)
  mse <- matrix(NA, 6, boot + jack + 1)
  mod <- rep(NA, boot + jack + 1)
  rownames(wghts) <- rownames(Inter) <- dimnames(bigdat)[[1]]
  rownames(mse) <- c("First model: Primary variables", "First model: Second variables",
                     "Secondary model: Primary variables", "Second model: Secondary variables",
                     "Third model: Primary variables", "Third model: Secondary variables")
  names(mod) <- colnames(wghts) <- colnames(Inter) <- colnames(mse) <- colnam
  inds <- 1:(boot + jack + 1)
  jack.lower <- 2
  jack.upper <- 1 + jack
  boot.lower <- 2 + jack
  boot.upper <- 1 + jack + boot
  
  for.max <- boot + jack + 1
  if(n.cores > 1 & boot + jack  > 0) {
    for.max <- 1
  }
  
  for (i in 1:for.max) {
    use.model.i <- use.model
    tmp <- proc.time()
    if (i == 1) {
      samp <- Int == int.val
      use <- rep(TRUE, n)
    }
    else if (grepl("Jack", colnam[i])) {
      samp <- Int == int.val
      g <- as.numeric(gsub("Jack", "", colnam[i]))
      use <- rep.G != g
      rep.meth <- "jackknife"
    }
    else if (grepl("Perm", colnam[i])) {
      g <- as.numeric(gsub("Perm", "", colnam[i]))
      if (!fin.boots) {
        samp <- base::sample(1:n, sum(Int == int.val), replace = FALSE,
                             prob = NULL)
        samp <- is.element(1:n, samp)
      }
      else {
        samp <- is.element(1:n, boots[, g])
      }
      use <- rep(TRUE, n)
      rep.meth <- "permutation"
    }
    Inter[samp & use, i] <- TRUE
    Inter[!samp & use, i] <- FALSE
    if (use.model.i == 1) {
      mod[i] <- "First"
      ws <- get.w.sub(newdat = newdat, newdat1 = newdat1,
                      end.pre = end.pre, samp = samp, use = use,
                      n = NROW(newdat), maxit = maxit, calfun = calfun,
                      bounds = bounds, epsilon = cal.epsilon, trim = trim,
                      qpmeth = qpmeth, scale.var = scale.var, printFlag = printFlag)
      if (ws$mse > cut.mse | is.na(ws$mse)) {
        if (use.backup) {
          cat.back <- ".  Using second model."
          cat.back1 <- ".  Used second model."
        }
        else {
          cat.back <- cat.back1 <- ".  Consider setting use.backup = TRUE."
        }
        if (i > 1) {
          back.state1 <- paste("First model was infeasible for ",
                               rep.meth, " group ", g, cat.back1, "\n",
                               sep = "")
        }
        else {
          back.state1 <- paste("First model is infeasible for primary weights",
                               cat.back, "\n", sep = "")
        }
        if (use.backup) {
          use.model.i <- 2
        }
      }
    }
    if (use.model.i == 2 & use.backup) {
      mod[i] <- "Second"
      ws <- get.w.sub(newdat = newdata, newdat1 = newdat1a,
                      end.pre = end.pre, samp = samp, use = use,
                      n = NROW(newdat), maxit = maxit, calfun = calfun,
                      bounds = bounds, epsilon = cal.epsilon, trim = trim,
                      qpmeth = qpmeth, scale.var = scale.var, printFlag = printFlag)
      if (ws$mse > cut.mse | is.na(ws$mse)) {
        if (i > 1) {
          back.state2 <- paste("Second model was infeasible for ",
                               rep.meth, " group ", g, ".  Used third model.",
                               "\n", sep = "")
          back.state3 <- paste("First and second model were infeasible for ",
                               rep.meth, " group ", g, ".  Used third model.",
                               "\n", sep = "")
        }
        else {
          back.state2 <- paste("Second model is infeasible for primary weights.  Will use third model.",
                               "\n", sep = "")
          back.state3 <- paste("First and second model are infeasible for primary weights.  Will use third model.",
                               "\n", sep = "")
        }
        use.model.i <- 3
      }
    }
    if (use.model.i == 3 & use.backup) {
      mod[i] <- "Third"
      ws <- get.w.sub(newdat = newdatb, newdat1 = newdat1b,
                      end.pre = end.pre, samp = samp, use = use,
                      n = NROW(newdat), maxit = maxit, calfun = calfun,
                      bounds = bounds, epsilon = cal.epsilon, trim = trim,
                      qpmeth = qpmeth, scale.var = scale.var, printFlag = printFlag)
    }
    if (back.state1 != "" & back.state2 != "") {
      back.state.final <- back.state3
    }
    else if (back.state1 != "" & back.state2 == "") {
      back.state.final <- back.state1
    }
    else if (back.state1 == "" & back.state2 != "") {
      back.state.final <- back.state2
    }
    else {
      back.state.final <- ""
    }
    back.state <- paste(back.state, back.state.final, sep = "")
    back.state1 <- back.state2 <- back.state3 <- ""
    mses <- get.mse(newdat, newdat1, samp, use, ws$wghts,
                    ws$wghts.init, ws$scale.by)
    msesa <- get.mse(newdata, newdat1a, samp, use, ws$wghts,
                     ws$wghts.init, ws$scale.by)
    msesb <- get.mse(newdatb, newdat1b, samp, use, ws$wghts,
                     ws$wghts.init, ws$scale.by)
    mse[1, i] <- mses$mse
    mse[2, i] <- mses$mse1
    mse[3, i] <- msesa$mse
    mse[4, i] <- msesa$mse1
    mse[5, i] <- msesb$mse
    mse[6, i] <- msesb$mse1
    wghts[, i] <- ws$wghts
    tmp <- proc.time() - tmp
    if (i == 1) {
      if (back.state != "")
        if(printFlag){message(back.state, appendLF = FALSE)}
      back.state <- ""
      use.model <- use.model.i
      if(printFlag){message("Created main weights for synthetic control: Time = ",
                            round(tmp[3], 2), "\n\n", sep = "", appendLF = FALSE)}
      if(printFlag){message("Matching summary for main weights:\n", appendLF = FALSE)}
      if (use.model.i == 1) {
        num.exact <- NCOL(newdat)
        if(length(newdat1) == 0) {
          num.prox <- 0
        } else {
          num.prox <- NCOL(newdat1)
        }
        printstuff <- mses$printstuff
        if(printFlag){message(paste0(utils::capture.output(round(printstuff,
                                                                 4)), collapse = "\n"), appendLF = FALSE)}
        if(printFlag){message("\n", appendLF = FALSE)}
      }
      else if (use.model.i == 2) {
        num.exact <- NCOL(newdata)
        if(length(newdat1a) == 0) {
          num.prox <- 0
        } else {
          num.prox <- NCOL(newdat1a)
        }
        printstuff <- msesa$printstuff
        if(printFlag){message(paste0(utils::capture.output(round(printstuff,
                                                                 4)), collapse = "\n"), appendLF = FALSE)}
        if(printFlag){message("\n", appendLF = FALSE)}
      }
      else if (use.model.i == 3) {
        num.exact <- NCOL(newdatb)
        if(length(newdat1b) == 0) {
          num.prox <- 0
        } else {
          num.prox <- NCOL(newdat1b)
        }
        printstuff <- msesb$printstuff
        if(printFlag){message(paste0(utils::capture.output(round(printstuff,
                                                                 4)), collapse = "\n"), appendLF = FALSE)}
        if(printFlag){message("\n", appendLF = FALSE)}
      }
      if (jack > 0) {
        if(printFlag){message("Calculating weights for jackknife replication groups...\n",
                              appendLF = FALSE)}
        tmp.jack <- proc.time()
      }
      else if (boot > 0) {
        if(printFlag){message("Calculating weights for permutation groups...\n",
                              appendLF = FALSE)}
        tmp.boot <- proc.time()
      }
    }
    else if (i >= jack.lower & i <= jack.upper) {
      if (i == jack.lower) {
        if(printFlag){message("Completed weights for jackknife group:\n",
                              i - 1, sep = "", appendLF = FALSE)}
      }
      else if ((i - 1)%%20 != 1 & i != jack.upper) {
        if(printFlag){message(", ", i - 1, sep = "", appendLF = FALSE)}
      }
      else if ((i - 1)%%20 == 1 & i != jack.upper) {
        if(printFlag){message(", \n", i - 1, sep = "", appendLF = FALSE)}
      }
      else if ((i - 1)%%20 != 1 & i == jack.upper) {
        if(printFlag){message(", ", i - 1, "\n", sep = "", appendLF = FALSE)}
      }
      else {
        if(printFlag){message(", \n", i - 1, "\n", sep = "", appendLF = FALSE)}
      }
      if (i == jack.upper) {
        if (i == jack.lower) {
        }
        if (back.state != "")
          if(printFlag){message(back.state, appendLF = FALSE)}
        back.state <- ""
        tmp.jack <- proc.time() - tmp.jack
        if(printFlag){message("Completed weights for all jackknife replication groups: Time = ",
                              round(tmp.jack[3], 2), "\n\n", sep = "", appendLF = FALSE)}
        if (boot > 0) {
          if(printFlag){message("Calculating weights for permutation groups...\n",
                                appendLF = FALSE)}
          tmp.boot <- proc.time()
        }
      }
    }
    else if (i >= boot.lower & i <= boot.upper) {
      if (i == boot.lower) {
        if(printFlag){message("Completed weights for permutation group:\n",
                              i - jack - 1, sep = "", appendLF = FALSE)}
      }
      else if ((i - jack - 1)%%20 != 1 & i != boot.upper) {
        if(printFlag){message(", ", i - jack - 1, sep = "", appendLF = FALSE)}
      }
      else if ((i - jack - 1)%%20 == 1 & i != boot.upper) {
        if(printFlag){message(", \n", i - jack - 1, sep = "", appendLF = FALSE)}
      }
      else if ((i - jack - 1)%%20 != 1 & i == boot.upper) {
        if(printFlag){message(", ", i - jack - 1, "\n", sep = "", appendLF = FALSE)}
      }
      else {
        if(printFlag){message(", \n", i - jack - 1, "\n", sep = "",
                              appendLF = FALSE)}
      }
      if (i == boot.upper) {
        if (i == boot.lower) {
          if(printFlag){message("\n", appendLF = FALSE)}
        }
        tmp.boot <- proc.time() - tmp.boot
        if(printFlag){message("Completed weights for all permutation groups: Time = ",
                              round(tmp.boot[3], 2), "\n\n", sep = "", appendLF = FALSE)}
      }
    }
  }
  
  if(boot + jack + 1 > for.max) {
    if(printFlag){message("Parallelizing with n.cores = ",n.cores,"...\n", sep = "", appendLF = FALSE)}
    requireNamespace("parallel", quietly = TRUE)
    cl <- parallel::makeCluster(n.cores)
    
    list.out <- parallel::parLapply(cl = cl, X = (for.max + 1):(boot + jack + 1), get.w.sub.par, 
                                    use.model, Int, int.val, colnam, rep.G, n, fin.boots, 
                                    boots, newdat, newdat1, newdata, newdat1a, newdatb, newdat1b, 
                                    end.pre, maxit, calfun, bounds, cal.epsilon, trim, 
                                    qpmeth, scale.var, printFlag = FALSE, cut.mse, use.backup, jack, 
                                    back.state, back.state1, back.state2, back.state3, tmp.jack, 
                                    boot, boot.lower, boot.upper, jack.lower, jack.upper, rownams = rownames(wghts), rownams1 = rownames(mse))
    
    parallel::stopCluster(cl)
    
    if(length(tmp.jack) == 0){tmp.jack <- tmp.boot}
    tmp.boot <- proc.time() - tmp.jack
    if(printFlag){message("Completed calculation of all replication weights: Time = ",
                          round(tmp.boot[3], 2), "\n\n", sep = "", appendLF = FALSE)}
    
    for(index in (for.max + 1):(boot + jack + 1)) {
      wghts[, index] <- list.out[[index - for.max]][[1]]
      Inter[, index] <- list.out[[index - for.max]][[2]]
      mse[, index] <- list.out[[index - for.max]][[3]]
      mod[index] <- list.out[[index - for.max]][[4]]
    }
  }
  
  keep.groups <- mse[2 * (use.model - 1) + 1,] < cut.mse & !is.na(mse[2 * (use.model - 1) + 1,])
  keep.groups[!grepl("Perm", colnames(mse))] <- TRUE
  
  if(printFlag & boot > 0){message("Removing ", sum(!keep.groups), 
                                   " permutation groups with pre-intervention MSE > cut.mse.\n\n",
                                   sep = "", appendLF = FALSE)}
  
  out <- list(Weights = wghts, Intervention = Inter, MSE = mse,
              Model = mod, Summary = printstuff, keep.groups = keep.groups,
              num.constr = c(num.exact = num.exact, num.prox = num.prox))
  return(out)
}



# Sub-function of get.w()
merge.dums <- function(dum, dum1) {
  nam1 <- names(dum)
  nam2 <- names(dum1)
  
  if (length(nam1) > 0) {
    nam <- union(nam1, nam2)
    
    dum.out <- list()
    dum.out[[length(nam)]] <- NA
    dum.new <- list()
    dum.new[[length(nam1)]] <- NA
    names(dum.new) <- nam1
    
    dum[[length(nam) + 1]] <- NA
    dum1[[length(nam) + 1]] <- NA
    
    max.l.dum <- 0
    
    for (i in 1:length(nam)) {
      here1 <- which(nam1 == nam[i])
      here2 <- which(nam2 == nam[i])
      
      if (length(here1) == 0 & length(here2) == 0) {
        dum.out[[i]] <- NULL
      } else if (length(here1) > 0 & length(here2) == 0) {
        max.l.dum <- max(max.l.dum, length(dum[[here1]]))
        dum.new[[here1]] <- sum(dum[[here1]])
        dum.out[[i]] <- dum[[here1]]
      } else if (length(here1) == 0 & length(here2) > 0) {
        dum.out[[i]] <- dum1[[here2]]
      } else {
        max.l.dum <- max(max.l.dum, length(dum[[here1]]))
        dum.new[[here1]] <- sum(dum[[here1]])
        
        dum1c <- cumsum(dum[[here1]])
        dum2c <- cumsum(dum1[[here2]])
        
        dumc <- union(dum1c, dum2c)
        dumc <- dumc[order(dumc, decreasing = FALSE)]
        
        if (length(dumc) > 1) {
          dumc <- c(dumc[1], dumc[2:length(dumc)] - dumc[1:(length(dumc) - 1)])
        }
        
        dum.out[[i]] <- dumc
      }
    }
    
    if (max.l.dum <= 1) {
      dum.new <- list()
    }
    
    names(dum.out) <- nam
  } else {
    dum.new <- dum
    dum.out <- dum1
  }
  return(list(dum = dum.new, dum1 = dum.out))
}



# Sub-function of get.w(); check if a constraint model has a feasible solution
is.feasible <- function(bigdat, covar.var, dum, Int, int.val = 1, end.pre, eps = 0.001) {
  n <- dim(bigdat)[1]
  
  newdat <- get.newdat(bigdat, dum = dum, covar.var = covar.var, end.pre = end.pre)[[1]]
  
  intdat <- newdat[Int == int.val, ]
  condat <- newdat[Int != int.val, ]
  targets <- colSums(intdat)
  is.sol <- check.feasible2(t(condat), targets, eps = eps)
  return(is.sol)
}


# Sub-function of is.feasible(); check if a constraint model has a feasible solution
check.feasible2 <- function(A, b, eps = 1e-07, M = 10000, meth = "LowRankQP") {
  if (NROW(A) <= NCOL(A)) {
    rem <- find.sing(tcrossprod(A))
  } else {
    rem <- NULL
  }
  leave <- setdiff(1:NROW(A), rem)
  A <- A[leave, ]
  b <- b[leave]
  A3 <- cbind(A, diag(NROW(A)), (-1) * diag(NROW(A)))
  a <- c(rep(0, NCOL(A)), rep(1, 2 * NROW(A)))
  
  if (meth == "LowRankQP") {
    requireNamespace("LowRankQP", quietly = TRUE)
    Vmat <- matrix(0, NCOL(A3), 1)
    uvec <- rep(M, NCOL(A3))
    sup.out <- utils::capture.output(all.root <- LowRankQP::LowRankQP(Vmat = Vmat, dvec = a, Amat = A3, bvec = b, uvec = uvec,
                                                                      method = "SMW", verbose = FALSE, niter = maxit))
    sol <- all.root$alpha
  }
  if (meth == "simplex") {
    requireNamespace("boot", quietly = TRUE)
    simp <- boot::simplex(a = a, A3 = A3, b3 = b)
    sol <- simp$soln
  }
  # w <- sol[1:NCOL(A)]
  sol <- sol[-(1:NCOL(A))]
  sol <- sqrt(mean(sol * sol))
  if (sol >= eps) {
    out <- FALSE
  } else {
    out <- TRUE
  }
  return(out)
}


# Sub-function of is.feasible(); reshape data inaccordance with a stated model structure 
# Produces two data frames: 
## newdat (data used for exact constraints) and
## newdat1 (data used for proximate constraints)
get.newdat <- function(bigdat, dum = NULL, dum1 = NULL, covar.var = NULL, covar.var1 = NULL, end.pre, time.names = NULL) {
  n <- dim(bigdat)[1]
  
  if (length(time.names) == 0) {
    time.names <- as.character(1:dim(bigdat)[3])
  }
  
  result.var <- names(dum)
  newdat <- list()
  
  if (length(result.var) > 0) {
    for (j in 1:length(result.var)) {
      dum.tmp <- dum[[j]]
      lows <- end.pre - cumsum(dum.tmp) + 1
      highs <- c(end.pre, lows[-length(lows)] - 1)
      lows.neg <- which(lows < 1)
      if(length(lows.neg) > 0) {
        if(highs[lows.neg[1]] >= 1) {
          lows[lows.neg[1]] <- 1
          lows.neg <- lows.neg[-1]
        }
        if(length(lows.neg) > 0) { 
          lows <- lows[-lows.neg]
          highs <- highs[-lows.neg]
        }
      }
      newdat[[j]] <- matrix(NA, dim(bigdat)[1], length(lows))
      colnames(newdat[[j]]) <- rep("", length(lows))
      for(i in 1:length(lows)) {
        newdat[[j]][,i] <- rowSums(as.matrix(bigdat[, result.var[j], lows[i]:highs[i]]))
        if (lows[i] == highs[i]) {
          colnames(newdat[[j]])[i] <- paste(result.var[j], ".", time.names[lows[i]], sep = "")
        } else {
          colnames(newdat[[j]])[i] <- paste(result.var[j], ".", time.names[lows[i]], ".", time.names[highs[i]], sep = "")
        }
      }
    }
  }
  
  covar.dat <- NULL
  if (length(covar.var) > 0) {
    covar.dat <- bigdat[, covar.var, 1]
    covar.dat <- as.matrix(covar.dat)
    colnames(covar.dat) <- covar.var
  }
  
  if (length(newdat) == 0) {
    newdat <- NULL
  } else {
    newdat <- data.frame(newdat, check.names = FALSE)
  }
  
  if (length(covar.dat) > 0 & length(newdat) > 0) {
    newdat <- data.frame(covar.dat, newdat, check.names = FALSE)
  } else if (length(covar.dat) > 0 & length(newdat) == 0) {
    newdat <- data.frame(covar.dat, check.names = FALSE)
  }
  
  if (length(newdat) > 0) {
    newdat <- data.frame(Intercept = rep(1, n), newdat, check.names = FALSE)
  } else {
    newdat <- data.frame(Intercept = rep(1, n))
  }
  
  newdat1 <- list()
  if (length(dum1) > 0) {
    result.var1 <- names(dum1)
    for (j in 1:length(result.var1)) {
      dum.tmp <- dum1[[j]]
      lows <- end.pre - cumsum(dum.tmp) + 1
      highs <- c(end.pre, lows[-length(lows)] - 1)
      lows.neg <- which(lows < 1)
      if(length(lows.neg) > 0) {
        if(highs[lows.neg[1]] >= 1) {
          lows[lows.neg[1]] <- 1
          lows.neg <- lows.neg[-1]
        }
        if(length(lows.neg) > 0) { 
          lows <- lows[-lows.neg]
          highs <- highs[-lows.neg]
        }
      }
      newdat1[[j]] <- matrix(NA, dim(bigdat)[1], length(lows))
      colnames(newdat1[[j]]) <- rep("", length(lows))
      for(i in 1:length(lows)) {
        newdat1[[j]][,i] <- rowSums(as.matrix(bigdat[, result.var1[j], lows[i]:highs[i]]))
        if (lows[i] == highs[i]) {
          colnames(newdat1[[j]])[i] <- paste(result.var1[j], ".", time.names[lows[i]], sep = "")
        } else {
          colnames(newdat1[[j]])[i] <- paste(result.var1[j], ".", time.names[lows[i]], ".", time.names[highs[i]], sep = "")
        }
      }
    }
  }
  
  covar.dat1 <- NULL
  if (length(covar.var1) > 0) {
    covar.dat1 <- bigdat[, covar.var1, 1]
    covar.dat1 <- as.matrix(covar.dat1)
    colnames(covar.dat1) <- covar.var1
  }
  
  if (length(newdat1) == 0) {
    newdat1 <- NULL
  } else {
    newdat1 <- data.frame(newdat1, check.names = FALSE)
  }
  
  if (length(covar.dat1) > 0 & length(newdat1) > 0) {
    newdat1 <- data.frame(covar.dat1, newdat1, check.names = FALSE)
  } else if (length(covar.dat1) > 0 & length(newdat1) == 0) {
    newdat1 <- data.frame(covar.dat1, check.names = FALSE)
  }
  
  if (length(newdat1) > 0) {
    newdat1 <- data.frame(newdat1, check.names = FALSE)
  } else {
    newdat1 <- NULL
  }
  
  return(list(newdat = newdat, newdat1 = newdat1))
  
}



# Secondary function used for calculating weights
get.w.sub <- function(newdat = NULL, newdat1 = NULL, bigdat = NULL, dum = NULL, dum1 = NULL, covar.var = NULL, covar.var1 = NULL,
                      end.pre, samp, use, n = NROW(newdat), maxit = 500, calfun = "raking", bounds = c(-Inf, Inf), epsilon = 1e-04, trim = NULL, qpmeth = "LowRankQP",
                      scale.var = "Intercept", printFlag = TRUE) {
  tmp <- proc.time()
  
  if (length(newdat) == 0) {
    newdat <- get.newdat(bigdat = bigdat, dum = dum, dum1 = dum1, covar.var = covar.var, covar.var1 = covar.var1, end.pre = end.pre)
    newdat1 <- newdat[[2]]
    newdat <- newdat[[1]]
  }
  
  mult <- sum(samp)/sum(use & samp)
  # mult <- 1
  
  if (length(newdat1) > 0) {
    newdat1 <- as.matrix(newdat1)
  }
  
  intdat <- newdat[samp & use, , drop = FALSE]
  scale.by <- sum(intdat[, scale.var])
  condat <- newdat[!samp & use, , drop = FALSE]
  
  alldat <- newdat[, , drop = FALSE]
  scale.by <- scale.by/sum(alldat[, scale.var])
  
  targets <- colSums(intdat)
  targets <- mult * targets
  # targets <- colSums(newdat[samp, , drop = FALSE])
  init <- rep(mult * NROW(intdat)/NROW(condat), NROW(condat))
  
  if (length(newdat1) > 0) {
    intdat1 <- newdat1[samp & use, , drop = FALSE]
    condat1 <- newdat1[!samp & use, , drop = FALSE]
    alldat1 <- newdat1[, , drop = FALSE]
    targets1 <- colSums(intdat1)
    targets1 <- mult * targets1
    # targets1 <- colSums(newdat1[samp, , drop = FALSE])
  }
  
  usevars <- colnames(condat)
  
  if (NCOL(as.matrix(condat)) <= NROW(as.matrix(condat))) {
    rem <- find.sing(crossprod(as.matrix(condat)))
  } else {
    rem <- NULL
  }
  keep <- !is.element(1:NCOL(condat), rem)
  
  form <- paste("~", paste(usevars[keep], collapse = "+", sep = ""), "-1", sep = "")
  form <- stats::formula(form)
  
  ws.init <- ws <- init
  
  if (NCOL(newdat) > 1) {
    tryCatch({
      caldesign <- survey::svydesign(ids = ~0, data = condat[, keep], weights = init)
      cali2 <- survey::calibrate(design = caldesign, maxit = maxit, formula = form, population = targets[keep], data = condat[,
                                                                                                                              keep], calfun = calfun, bounds = bounds, force = TRUE, epsilon = epsilon)
      ws.init <- ws <- stats::weights(cali2)
    }, error = function(e, printFlag = printFlag) {
      if(printFlag){message("ERROR :", conditionMessage(e), "\n", appendLF = FALSE)}
    })
  }
  
  if (length(newdat1) > 0) {
    tmp1 <- proc.time() - tmp
    tmp2 <- proc.time()
    
    condat2 <- condat1
    targets2 <- targets1
    
    ws <- my.qp(b.init = ws.init, X = t(condat2), a = targets2, Y = t(condat[, keep, drop = FALSE]), c = targets[keep], qpmeth = qpmeth, printFlag = printFlag)
    tmp2 <- proc.time() - tmp2
  }
  
  if (length(trim) > 0) {
    if (length(trim) == 1) {
      trim = c(trim, 1 - trim)
    }
    cali2 <- survey::trimWeights(cali2, upper = stats::quantile(ws, max(trim)), lower = stats::quantile(ws, min(trim)))
    ws <- stats::weights(cali2)
  }
  
  mse1 <- mean((colSums(ws * condat[, , drop = FALSE]) - targets)^2)
  wghts1 <- rep(NA, n)
  wghts1[!samp & use] <- ws
  wghts1[samp & use] <- mult
  wghts.init1 <- rep(NA, n)
  wghts.init1[!samp & use] <- ws.init
  wghts.init1[samp & use] <- mult
  
  out <- list(wghts = wghts1, wghts.init = wghts.init1, mse = mse1, scale.by = scale.by)
  
  return(out)
}



get.w.sub.par <- function (X, use.model, Int, int.val, colnam, rep.G, n, fin.boots, 
                           boots, newdat, newdat1, newdata, newdat1a, newdatb, newdat1b, 
                           end.pre, maxit, calfun, bounds, cal.epsilon, trim, 
                           qpmeth, scale.var, printFlag, cut.mse, use.backup, jack, 
                           back.state, back.state1, back.state2, back.state3, tmp.jack, 
                           boot, boot.lower, boot.upper, jack.lower, jack.upper, rownams, rownams1) {
  
  tmp.boot <- tmp.jack
  i <- X
  wghts.out <- Inter.out <- rep(NA, n)
  names(wghts.out) <- names(Inter.out) <- rownams
  mse.out <- rep(NA, 6)
  names(mse.out) <- rownams1
  
  use.model.i <- use.model
  tmp <- proc.time()
  if (i == 1) {
    samp <- Int == int.val
    use <- rep(TRUE, n)
  }
  else if (grepl("Jack", colnam[i])) {
    samp <- Int == int.val
    g <- as.numeric(gsub("Jack", "", colnam[i]))
    use <- rep.G != g
    rep.meth <- "jackknife"
  }
  else if (grepl("Perm", colnam[i])) {
    g <- as.numeric(gsub("Perm", "", colnam[i]))
    if (!fin.boots) {
      samp <- base::sample(1:n, sum(Int == int.val), replace = FALSE,
                           prob = NULL)
      samp <- is.element(1:n, samp)
    }
    else {
      samp <- is.element(1:n, boots[, g])
    }
    use <- rep(TRUE, n)
    rep.meth <- "permutation"
  }
  Inter.out[samp & use] <- TRUE
  Inter.out[!samp & use] <- FALSE
  if (use.model.i == 1) {
    mod.out <- "First"
    ws <- microsynth:::get.w.sub(newdat = newdat, newdat1 = newdat1,
                                 end.pre = end.pre, samp = samp, use = use,
                                 n = NROW(newdat), maxit = maxit, calfun = calfun,
                                 bounds = bounds, epsilon = cal.epsilon, trim = trim,
                                 qpmeth = qpmeth, scale.var = scale.var, printFlag = printFlag)
    if (ws$mse > cut.mse | is.na(ws$mse)) {
      if (use.backup) {
        cat.back <- ".  Using second model."
        cat.back1 <- ".  Used second model."
      }
      else {
        cat.back <- cat.back1 <- ".  Consider setting use.backup = TRUE."
      }
      if (i > 1) {
        back.state1 <- paste("First model was infeasible for ",
                             rep.meth, " group ", g, cat.back1, "\n",
                             sep = "")
      }
      else {
        back.state1 <- paste("First model is infeasible for primary weights",
                             cat.back, "\n", sep = "")
      }
      if (use.backup) {
        use.model.i <- 2
      }
    }
  }
  if (use.model.i == 2 & use.backup) {
    mod.out <- "Second"
    ws <- microsynth:::get.w.sub(newdat = newdata, newdat1 = newdat1a,
                                 end.pre = end.pre, samp = samp, use = use,
                                 n = NROW(newdat), maxit = maxit, calfun = calfun,
                                 bounds = bounds, epsilon = cal.epsilon, trim = trim,
                                 qpmeth = qpmeth, scale.var = scale.var, printFlag = printFlag)
    if (ws$mse > cut.mse | is.na(ws$mse)) {
      if (i > 1) {
        back.state2 <- paste("Second model was infeasible for ",
                             rep.meth, " group ", g, ".  Used third model.",
                             "\n", sep = "")
        back.state3 <- paste("First and second model were infeasible for ",
                             rep.meth, " group ", g, ".  Used third model.",
                             "\n", sep = "")
      }
      else {
        back.state2 <- paste("Second model is infeasible for primary weights.  Will use third model.",
                             "\n", sep = "")
        back.state3 <- paste("First and second model are infeasible for primary weights.  Will use third model.",
                             "\n", sep = "")
      }
      use.model.i <- 3
    }
  }
  if (use.model.i == 3 & use.backup) {
    mod.out <- "Third"
    ws <- microsynth:::get.w.sub(newdat = newdatb, newdat1 = newdat1b,
                                 end.pre = end.pre, samp = samp, use = use,
                                 n = NROW(newdat), maxit = maxit, calfun = calfun,
                                 bounds = bounds, epsilon = cal.epsilon, trim = trim,
                                 qpmeth = qpmeth, scale.var = scale.var, printFlag = printFlag)
  }
  if (back.state1 != "" & back.state2 != "") {
    back.state.final <- back.state3
  }
  else if (back.state1 != "" & back.state2 == "") {
    back.state.final <- back.state1
  }
  else if (back.state1 == "" & back.state2 != "") {
    back.state.final <- back.state2
  }
  else {
    back.state.final <- ""
  }
  back.state <- paste(back.state, back.state.final, sep = "")
  back.state1 <- back.state2 <- back.state3 <- ""
  mses <- microsynth:::get.mse(newdat, newdat1, samp, use, ws$wghts,
                               ws$wghts.init, ws$scale.by)
  msesa <- microsynth:::get.mse(newdata, newdat1a, samp, use, ws$wghts,
                                ws$wghts.init, ws$scale.by)
  msesb <- microsynth:::get.mse(newdatb, newdat1b, samp, use, ws$wghts,
                                ws$wghts.init, ws$scale.by)
  mse.out[1] <- mses$mse
  mse.out[2] <- mses$mse1
  mse.out[3] <- msesa$mse
  mse.out[4] <- msesa$mse1
  mse.out[5] <- msesb$mse
  mse.out[6] <- msesb$mse1
  wghts.out[] <- ws$wghts
  tmp <- proc.time() - tmp
  if (i == 1) {
    if (back.state != "")
      if(printFlag){message(back.state, appendLF = FALSE)}
    back.state <- ""
    use.model <- use.model.i
    if(printFlag){message("Created main weights for synthetic control: Time = ",
                          round(tmp[3], 2), "\n\n", sep = "", appendLF = FALSE)}
    if(printFlag){message("Matching summary for main weights:\n", appendLF = FALSE)}
    if (use.model.i == 1) {
      num.exact <- NCOL(newdat)
      if(length(newdat1) == 0) {
        num.prox <- 0
      } else {
        num.prox <- NCOL(newdat1)
      }
      printstuff <- mses$printstuff
      if(printFlag){message(paste0(utils::capture.output(round(printstuff,
                                                               4)), collapse = "\n"), appendLF = FALSE)}
      if(printFlag){message("\n", appendLF = FALSE)}
    }
    else if (use.model.i == 2) {
      num.exact <- NCOL(newdata)
      if(length(newdat1a) == 0) {
        num.prox <- 0
      } else {
        num.prox <- NCOL(newdat1a)
      }
      printstuff <- msesa$printstuff
      if(printFlag){message(paste0(utils::capture.output(round(printstuff,
                                                               4)), collapse = "\n"), appendLF = FALSE)}
      if(printFlag){message("\n", appendLF = FALSE)}
    }
    else if (use.model.i == 3) {
      num.exact <- NCOL(newdatb)
      if(length(newdat1b) == 0) {
        num.prox <- 0
      } else {
        num.prox <- NCOL(newdat1b)
      }
      printstuff <- msesb$printstuff
      if(printFlag){message(paste0(utils::capture.output(round(printstuff,
                                                               4)), collapse = "\n"), appendLF = FALSE)}
      if(printFlag){message("\n", appendLF = FALSE)}
    }
    if (jack > 0) {
      if(printFlag){message("Calculating weights for jackknife replication groups...\n",
                            appendLF = FALSE)}
      tmp.jack <- proc.time()
    }
    else if (boot > 0) {
      if(printFlag){message("Calculating weights for permutation groups...\n",
                            appendLF = FALSE)}
      tmp.boot <- proc.time()
    }
  }
  else if (i >= jack.lower & i <= jack.upper) {
    if (i == jack.lower) {
      if(printFlag){message("Completed weights for jackknife group:\n",
                            i - 1, sep = "", appendLF = FALSE)}
    }
    else if ((i - 1)%%20 != 1 & i != jack.upper) {
      if(printFlag){message(", ", i - 1, sep = "", appendLF = FALSE)}
    }
    else if ((i - 1)%%20 == 1 & i != jack.upper) {
      if(printFlag){message(", \n", i - 1, sep = "", appendLF = FALSE)}
    }
    else if ((i - 1)%%20 != 1 & i == jack.upper) {
      if(printFlag){message(", ", i - 1, "\n", sep = "", appendLF = FALSE)}
    }
    else {
      if(printFlag){message(", \n", i - 1, "\n", sep = "", appendLF = FALSE)}
    }
    if (i == jack.upper) {
      if (i == jack.lower) {
      }
      if (back.state != "")
        if(printFlag){message(back.state, appendLF = FALSE)}
      back.state <- ""
      tmp.jack <- proc.time() - tmp.jack
      if(printFlag){message("Completed weights for all jackknife replication groups: Time = ",
                            round(tmp.jack[3], 2), "\n\n", sep = "", appendLF = FALSE)}
      if (boot > 0) {
        if(printFlag){message("Calculating weights for permutation groups...\n",
                              appendLF = FALSE)}
        tmp.boot <- proc.time()
      }
    }
  }
  else if (i >= boot.lower & i <= boot.upper) {
    if (i == boot.lower) {
      if(printFlag){message("Completed weights for permutation group:\n",
                            i - jack - 1, sep = "", appendLF = FALSE)}
    }
    else if ((i - jack - 1)%%20 != 1 & i != boot.upper) {
      if(printFlag){message(", ", i - jack - 1, sep = "", appendLF = FALSE)}
    }
    else if ((i - jack - 1)%%20 == 1 & i != boot.upper) {
      if(printFlag){message(", \n", i - jack - 1, sep = "", appendLF = FALSE)}
    }
    else if ((i - jack - 1)%%20 != 1 & i == boot.upper) {
      if(printFlag){message(", ", i - jack - 1, "\n", sep = "", appendLF = FALSE)}
    }
    else {
      if(printFlag){message(", \n", i - jack - 1, "\n", sep = "",
                            appendLF = FALSE)}
    }
    if (i == boot.upper) {
      if (i == boot.lower) {
        if(printFlag){message("\n", appendLF = FALSE)}
      }
      tmp.boot <- proc.time() - tmp.boot
      if(printFlag){message("Completed weights for all permutation groups: Time = ",
                            round(tmp.boot[3], 2), "\n\n", sep = "", appendLF = FALSE)}
    }
  }
  return(list(wghts = wghts.out, Inter = Inter.out, mse = mse.out, mod = mod.out))
}



# Determine pre-intervention MSEs and create balance table
get.mse <- function(newdat, newdat1 = NULL, samp, use, ws, ws.init, scale.by) {
  
  mult <- sum(samp)/sum(use & samp)
  
  intdat <- newdat[samp & use, , drop = FALSE]
  condat <- newdat[!samp & use, , drop = FALSE]
  alldat <- newdat[, , drop = FALSE]
  targets <- colSums(intdat)
  targets <- mult * targets
  ws <- ws[!samp & use]
  mse <- mean((colSums(ws * condat) - targets)^2)
  
  mse1 <- NA
  if (length(newdat1) != 0) {
    intdat1 <- newdat1[samp & use, , drop = FALSE]
    condat1 <- newdat1[!samp & use, , drop = FALSE]
    alldat1 <- newdat1[, , drop = FALSE]
    targets1 <- colSums(intdat1)
    targets1 <- mult * targets1
    ws.init <- ws.init[!samp & use]
    mse1 <- mean((colSums(ws * condat1) - targets1)^2)
  }
  
  if (length(newdat1) == 0) {
    printstuff <- cbind(Targets = targets, Weighted.Control = colSums((ws) * (condat)), All.scaled = (scale.by) * colSums(alldat))
  } else {
    # printstuff <- cbind(Targets = c(targets, targets1), Initial.Weighted.Control = colSums(ws.init * cbind(condat, condat1)),
    # Final.Weighted.Control = colSums(ws * cbind(condat, condat1)), All.scaled = scale.by * colSums(cbind(alldat, alldat1)))
    printstuff <- cbind(Targets = c(targets, targets1), Final.Weighted.Control = colSums(ws * cbind(condat, condat1)), All.scaled = scale.by *
                          colSums(cbind(alldat, alldat1)))
  }
  
  return(list(mse = mse, mse1 = mse1, printstuff = printstuff))
}


# Function used to solve quadratic optimization when proximate constraints are invoked
my.qp <- function(b.init, X, Y, a, c, M = 10000, qpmeth = "LowRankQP", maxit = 1000, printFlag = TRUE) {
  q <- NCOL(X)
  n <- NROW(Y)
  
  a1 <- iso.nam(a, sep = ".")
  c1 <- iso.nam(c, sep = ".")
  
  a1[a1 == 0] <- 1
  c1[c1 == 0] <- 1
  
  a1 <- sqrt(a1)
  c1 <- sqrt(c1)
  
  X <- X/a1
  Y <- Y/c1
  a <- a/a1
  c <- c/c1
  
  if (qpmeth == "nleqslv") {
    requireNamespace("nleqslv", quietly = TRUE)
    XtX <- crossprod(X)
    
    M <- 1/M
    jacobian <- function(all) {
      b <- all[1:q]
      lambda <- all[(q + 1):(q + n)]
      lowerright <- matrix(0, n, n)
      upperleft <- diag(1/b) + XtX
      upperright <- t(Y)
      lowerleft <- Y
      out <- rbind(cbind(upperleft, upperright), cbind(lowerleft, lowerright))
      return((out))
    }
    
    f1 <- function(all) {
      b <- all[1:q]
      lambda <- all[(q + 1):(q + n)]
      out1 <- (M * log(b)) + crossprod(X, (X %*% b - a)) - crossprod(Y, lambda)
      out2 <- Y %*% b - c
      return(c(out1, out2))
    }
    b.init[b.init < 1e-06] <- 1e-06
    lambda.init <- solve(tcrossprod(Y), Y) %*% ((M * log(b.init)) + crossprod(X, (X %*% b.init - a)))
    all.init <- c(b.init, lambda.init)
    
    all.root <- nleqslv::nleqslv(all.init, fn = f1, jac = NULL, method = c("Broyden"), global = c("gline"), xscalm = c("fixed"),
                                 control = list(maxit = maxit, xtol = 1e-11, ftol = 1e-11, btol = 1e-06, cndtol = 1e-13))
    
    if(printFlag){message("Number of Jacobian evaluations = ", all.root$njcnt, ". \n", sep = "", appendLF = FALSE)}
    if(printFlag){message("Number of function evaluations = ", all.root$nfcnt, ". \n", sep = "", appendLF = FALSE)}
    if(printFlag){message("Number of iterations = ", all.root$iter, ". \n", sep = "", appendLF = FALSE)}
    all.root <- all.root$x
    
    b <- all.root[1:q]
    b[b < 0] <- 0
  } else if (qpmeth == "LowRankQP") {
    requireNamespace("LowRankQP", quietly = TRUE)
    
    if (NROW(Y) <= NCOL(Y)) {
      rem <- find.sing(tcrossprod(Y))
    } else {
      rem <- NULL
    }
    leave <- setdiff(1:NROW(Y), rem)
    sup.out <- utils::capture.output(all.root <- LowRankQP::LowRankQP(Vmat = t(X), dvec = -crossprod(X, a), Amat = Y[leave, , drop = FALSE],
                                                                      bvec = c[leave], uvec = rep(M, q), method = "SMW", verbose = FALSE, niter = maxit))
    b <- all.root$alpha
  } else if (qpmeth == "ipop") {
    requireNamespace("kernlab", quietly = TRUE)
    
    all.root <- kernlab::ipop(c = -crossprod(X, a), H = crossprod(X), b = c, A = Y, l = rep(0, q), r = rep(0, n), u = rep(M, q))
    b <- all.root$alpha
  } else {
    stop("qpmeth not recognized.")
  }
  
  X <- X * a1
  Y <- Y * c1
  a <- a * a1
  c <- c * c1
  
  return(b)
}

# Sub-function of my.qp()
iso.nam <- function(a, sep = "-") {
  nams <- names(a)
  locs <- regexpr(sep, nams, fixed = TRUE)
  outs <- nams
  outs[locs > 1] <- substr(outs[locs > 1], 1, locs[locs > 1] - 1)
  
  out <- rep(NA, length(outs))
  for (i in names(table(outs))) {
    out[outs == i] <- mean(a[outs == i])
  }
  return(out)
}
