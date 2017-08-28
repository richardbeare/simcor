## Original versions of functions

noisecor.orig <- function(cormat,
                          epsilon = .01,
                          eidim = 2) 
{
  ndim = dim(cormat)[1]
  diag(cormat) <- 1 - epsilon
  ### adding noise to the correlation matrix
  
  eivect <- c()
  for (i in 1:ndim) {
    ei <- runif(eidim,-1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei ^ 2)))
  }
  bigE <- t(eivect) %*% eivect
  cor.nz <- cormat + bigE
  return(cor.nz)
  
}

simcor.orig <- function (k = 6,
                         size = c(10, 5, 8, 2, 15, 50),
                         rho = c(0.7,
                                 0.7, 0.5, 0.9, 0.85, 0.4),
                         delta = 0.39,
                         epsilon = 0.99 -
                           max(rho),
                         eidim = 2)
{
  ndim <- sum(size)
  bigcor <- matrix(rep(delta, ndim * ndim), ncol = ndim)
  for (i in 1:k) {
    cor <- matrix(rep(rho[i], size[i] * size[i]), ncol = size[i])
    if (i == 1) {
      bigcor[1:size[1], 1:size[1]] <- cor
    }
    if (i != 1) {
      bigcor[(sum(size[1:(i - 1)]) + 1):sum(size[1:i]),
             (sum(size[1:(i - 1)]) + 1):sum(size[1:i])] <- cor
    }
  }
  diag(bigcor) <- 1 - epsilon
  
  eivect <- c()
  for (i in 1:ndim) {
    ei <- runif(eidim,-1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei ^ 2)))
  }
  bigE <- t(eivect) %*% eivect
  cor.nz <- bigcor + bigE
  return(cor.nz)
}

simcorTop.orig <-
  function(k = 6,
           size = c(10, 5, 8, 2, 15, 50),
           rho = c(.7, .7, .5, .9, .85, .4),
           epsilon = .01,
           eidim = 2) 
  {
    ndim <- sum(size)# dim of correlation matrix
    bigcor <- matrix(rep(0, ndim * ndim), ncol = ndim)
    
    ### generating the basic correlation matrix
    for (i in 1:k) {
      top <- c(1, rho[i] ^ (seq(1:(size[i] - 1))))
      cor <- toeplitz(top)
      
      if (i == 1) {
        bigcor[1:size[1], 1:size[1]] <- cor
      }
      if (i != 1) {
        bigcor[(sum(size[1:(i - 1)]) + 1):sum(size[1:i]),
               (sum(size[1:(i - 1)]) + 1):sum(size[1:i])] <- cor
      }
    }
    
    diag(bigcor) <- 1 - epsilon
    ### adding noise to the correlation matrix
    eivect <- c()
    for (i in 1:ndim) {
      ei <- runif(eidim,-1, 1)
      eivect <- cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei ^ 2)))
    }
    
    
    bigE <- t(eivect) %*% eivect
    cor.nz <- bigcor + bigE
    return(cor.nz)
    
  }
simcor.H.orig <- function(k = 6,
                          size = c(10, 5, 8, 7, 15, 50),
                          rho = rbind(c(.9, .7), c(.7, .7), c(.7, .2), c(.5, .3), c(.9, .85), c(.3, .2)),
                          power = 1,
                          epsilon = .08,
                          eidim = 2) 
{
  ndim <- sum(size)# dim of correlation matrix
  bigcor <- matrix(rep(0, ndim * ndim), ncol = ndim)
  
  ### generating the basic correlation matrix
  
  
  for (i in 1:(k)) {
    cor <- toeplitz(rho.func.orig(rho[i, 1], rho[i, 2], power, size[i]))
    
    if (i == 1) {
      bigcor[1:size[1], 1:size[1]] <- cor
    }
    if (i != 1) {
      bigcor[(sum(size[1:(i - 1)]) + 1):sum(size[1:i]),
             (sum(size[1:(i - 1)]) + 1):sum(size[1:i])] <- cor
    }
  }
  diag(bigcor) <- 1 - epsilon
  
  
  ### adding noise to the correlation matrix
  
  eivect <- c()
  for (i in 1:ndim) {
    ei <- runif(eidim,-1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei ^ 2)))
  }
  
  
  bigE <- t(eivect) %*% eivect
  cor.nz <- bigcor + bigE
  return(cor.nz)
}


rho.func.orig <- function(r.max, r.min, power, p) 
{
  rhovec <- c()
  
  rhovec[1] <- 1
  for (i in 2:p) {
    rhovec[i] <- r.max - ((i - 2) / (p - 2)) ^ power * (r.max - r.min)
  }
  return(rhovec)
}
