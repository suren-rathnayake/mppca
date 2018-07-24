mppca <- function(Y, g, q, clust, sigma_type = "common",
                              D_type = "common", ...) {
  
  if (class(clust) == "factor") {
    clust <- as.numeric(clust)
  }
  if (class(clust) != "numeric") {
    stop("clust must be factor or integer vector")
  }
  if (any(is.na(clust))) {
    stop("clust must not contain NA's")
  }
  if (any(is.na(Y))) {
    ERR <- "`Y' has missing value."
    class(ERR) <- "error"
    return(ERR)
  }
  
  if (!any(is.numeric(Y))) {
    ERR <- "`Y' has a non-numeric element."
    class(ERR) <- "error"
    return(ERR)
  }
  
  if (ncol(Y) <= q) {
    ERR <- "The number of factors must be less than the number of variables."
    class(ERR) <- "error"
    return(ERR)
  }
  if (q < 0) {
    ERR <- "q must be a positive integer."
    class(ERR) <- "error"
    return(ERR)
  }
  
  if (q != round(q)) {
    ERR <- "q must be a positive integer."
    class(ERR) <- "error"
    return(ERR)
  }
  if (g < 0) {
    ERR <- "g must be a positive integer."
    class(ERR) <- "error"
    return(ERR)
  }
  if (g != round(g)) {
    ERR <- "g must be a positive integer."
    class(ERR) <- "error"
    return(ERR)
  }
  if ((sigma_type == 'common') && (D_type == 'unique')) {
    ERR <- "D_type = 'unique' not available with sigma_type = 'common'."
    class(ERR) <- "error"
    return(ERR)
  }
  if ((sigma_type != "common") && (sigma_type != "unique")) {
    ERR <- "sigma_type needs to be either 'unique' or 'common'."
    class(ERR) <- "error"
    return(ERR)
  }
  if ((D_type != "common") && (D_type != "unique")) {
    ERR <- "D_type needs to be either 'unique' or 'common'."
    class(ERR) <- "error"
    return(ERR)
  }
  
  p <- ncol(Y)
  n <- nrow(Y)
  B <- array(NA, c(p, q, g))
  pivec <- array(NA, c(g, 1))
  mu <- array(NA, c(p, g))
  
  if (D_type == 'common') {
    
    D <- diag(diag(cov(Y)))
    Di.sqrt <- diag(sqrt(diag(D)))
    inv.Di.sqrt <- diag(1 / diag(Di.sqrt))
    
  } else if (D_type == 'unique') {
    
    D <- array(NA, c(p, p, g))
  }
  
  for(i in 1 : g) {
    
    indices <- which(clust == i)
    pivec[i] <- length(indices) / n
    mu[, i] <- apply(Y[indices, ], 2, mean)
    Si <- cov(Y[indices, ])
    if (D_type == 'unique'){
      D[,, i] <- diag(diag(Si))
      Di.sqrt <- diag(sqrt(diag(D[,, i])))
      inv.Di.sqrt <- diag(1 / diag(Di.sqrt))
    }
    
    eig.list <- try(eigen(inv.Di.sqrt %*% Si %*% inv.Di.sqrt), TRUE)
    if (class(eig.list) == "try-error")
      break
    
    H <- eig.list$vectors
    sort.lambda <- sort(eig.list$values, decreasing = TRUE,
                        index.return = TRUE)
    lambda <- sort.lambda$x
    ix.lambda   <- sort.lambda$ix
    if (q == p) {
      sigma2 <- 0
    } else {
      sigma2 <- mean(lambda[(q + 1) : p])
    }
    if (q == 1) {
      B[,, i] <- Di.sqrt %*% H[, ix.lambda[1 : q]] %*%
        diag((lambda[1 : q] - sigma2), q)
    } else {
      B[,, i] <- Di.sqrt %*% H[, ix.lambda[1 : q]] %*%
        diag((lambda[1 : q] - sigma2))
    }
  }
  if (sigma_type == 'common')
    B <- apply(B, c(1, 2), sum) / g
  
  model <- list(g = g, q = q, pivec = pivec, B = B, mu = mu, D = D,
                sigma_type = sigma_type, D_type = D_type)
  
  class(model) <- "mppca"
  return(model)
}
