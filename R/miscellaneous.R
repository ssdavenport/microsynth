
# Assign object to class "microsynth"
makemicrosynth <- function(x) {
  class(x) <- c("microsynth")
  x
}


# Clear global variable
maxit <- NULL


# Sub-function of get.w(), get.stats1(), get.stats1.sub(); Assign replication groups for jackknife procedures
assign.groups <- function(strata = NULL, n = length(strata), G = min(table(strata))) {

  states <- names(table(strata))

  Gs <- base::sample(1:G, G)
  Gs1 <- Gs

  rep.G <- rep(NA, n)

  for (i in 1:length(states)) {
    here <- which(strata == states[i])

    n <- length(here)
    samp <- base::sample(1:n, n)

    J <- floor(n/G)

    rep.G1 <- rep(NA, n)

    for (g in 1:G) {
      here1 <- samp[1:J]
      samp <- samp[-(1:J)]
      rep.G1[here1] <- g
    }

    if (length(is.na(rep.G1)) > 0) {
      here2 <- which(is.na(rep.G1))
      if (length(here2) > length(Gs1)) {
        Gs1 <- c(Gs1, Gs)
      }
      rep.G1[here2] <- Gs1[1:length(here2)]
      Gs1 <- Gs1[-(1:length(here2))]
    }

    rep.G[here] <- rep.G1
  }

  return(rep.G)
}


# Common sub-function; Find columns of a matrix that create singularity
find.sing <- function(X) {
  # X must be a square matrix
  X <- as.matrix(X)
  if (NROW(X) == 1) {
    Z <- X/X[1,1]
  } else {
    #  Scale the matrix to avoid computational issues when applying rref
    D <- diag((diag(X)^-.5))
    X <- D %*% X %*% D
    Z <- pracma::rref(X)
  }
  dimnames(Z) <- dimnames(X)

  rem <- NULL
  cont <- TRUE
  while (cont) {
    Z.tmp <- diag(Z)
    if (sum(Z.tmp == 0) == 0) {
      cont <- FALSE
    }
    if (cont) {
      Z.tmp <- min(which(Z.tmp == 0))
      Z <- Z[-NROW(Z), -Z.tmp]
      rem <- c(rem, Z.tmp + length(rem))
    }
  }
  return(rem)
}
