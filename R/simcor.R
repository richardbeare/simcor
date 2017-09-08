

#' @import stats

##########################################################
#### R code for adding noise to correlation matrices  ####
#### J. Hardin, S.R. Garcia, D. Golan                 ####
#### last updated 1/18/2013                           ####
##########################################################

#################
##   EXAMPLE   ##
#################


# adding noise to a 15x15 identity matrix
function()
{
  ## this stuff will be examples
noise.iden <- noisecor(diag(15), epsilon = .5, eidim = 2)


# adding noise to a 2 block constant correlation matrix
rho.blck = c(.9, .2)
eps.blck = 0.99 - max(rho.blck)
samp.blck <-
  simcor(
    k = 2,
    size = c(10, 5),
    rho = c(0.9, 0.2),
    delta = 0.39,
    epsilon = eps.blck,
    eidim = 2
  )


# adding noise to a 2 block correlation matrix with an AR(1) strucutre
rho.Top = c(.7, .9)
eps.Top = (1 - max(rho.Top)) / (1 + max(rho.Top)) - .01
samp.Top <-
  simcorTop(
    k = 2,
    size = c(10, 5),
    rho = c(.7, .9),
    epsilon = eps.Top,
    eidim = 2
  )


# adding noise to a 2 block correlation matrix with a Hub structure
rho.Hub = c(.9, .7)
tau.Hub = c((.9 - .7) / (10 - 2), (.7 - .6) / (5 - 2))
eps.Hub = min(1 - (rho.Hub) - 0.75 * (tau.Hub)) - .01
samp.Hub <-
  simcor.H(
    k = 2,
    size = c(10, 5),
    rho = rbind(c(.9, .7), c(.7, .6)),
    power = 1,
    epsilon = eps.Hub,
    eidim = 2
  )

}


########################
## The Basic Function ##
########################
#' Add noise to any user specified correlation matrix
#' @param cormat the correlation matrix to which noise is to be added.
#' @param epsilon epsilon maximum entry-wise random noise.
#' @param eidim dimension of the noise.
#' @description From Appendix A2. This is the core algorithm for adding noise to an
#' arbitrary positive semi-definite matrix. The choice of eidim (M in the paper)
#' provides flexibilty. If ei (M) is small (between 2 and 5) then many of the dot
#' product terms will be large, yielding a very noisy coefficient matrix. If M is large
#' the computation is still cheap. Other choices are discussed in the paper. This version
#' has been modified to remove loops.
#' @return The noise corrupted correlation matrix.
#' @export
#' @examples
#'  noise.iden <- noisecor(diag(15), epsilon = .5, eidim = 2)
noisecor <- function(cormat,
                          epsilon = .01,
                          eidim = 2) 
{
  ndim = dim(cormat)[1]
  diag(cormat) <- 1 - epsilon
  
  ### adding noise to the correlation matrix
  eivect <- runif(ndim*eidim, -1, 1)
  dim(eivect) <- c(eidim, ndim)
  denom <- sqrt(colSums(eivect^2))
  ## transpose due to rowwise ops - could use sweep.
  eivect <- sqrt(epsilon) * t(t(eivect)/denom)
  bigE <- t(eivect) %*% eivect
  cor.nz <- cormat + bigE
  return(cor.nz)
}



#######################################
##  Simulating different correlation ##
##  structures as in Hardin et al.   ##
#######################################


################################
##  Simulating the Constant Correlation
################################

# this function simulates a block correlation matrix
# The size and base correlation for each block is user specified
# There is an additional delta parameter for the off diagonal correlations


# k is the number of groups
# size is a vector of length k specifying the size of each group
# rho is a vector of length k specifying base correlation values
# epsilon <- 0.99 - max(rho)
# eidim is the space from which the noise is generated, the smaller the more noise
# delta is the correlation of the off diagonal blocks

#' Block correlation noise matrix
#' @param k Number of groups.
#' @param size Vector length k specifying the size of each group.
#' @param rho vector of length k specifying base correlation values
#' @param delta correlation of the off diagonal blocks (between group noise)
#' @param epsilon maximum entry-wise random noise.
#' @param eidim passed to noise corr (or should be)
#' @export
simcor <- function (k = 6,
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
  ss <- c(0,cumsum(size))
  for (i in 1:k) {
    cor <- matrix(rep(rho[i], size[i] * size[i]), ncol = size[i])
    ti <- (ss[i]+1):ss[i+1]
    bigcor[ti, ti] <- cor
  }
  cor.nz <- noisecor(bigcor, epsilon, eidim)
  return(cor.nz)
}



################################
##  Simulating the Toeplitz Matrix
################################

#' Toeplitz correlation structure
#' @details This function calculates an AR(1) Toeplitz matrix
#' with a block diagonal structure.  The size and base correlation for each
#' block is user specified
#' @param k is the number of groups
#' @param size is a vector of length k specifying the size of each group
#' @param rho is a vector of length k specifying base correlation values
#' @param maximum entry-wise random noise.
#' @param eidim is the space from which the noise is generated, the smaller the more noise
#' @export
simcorTop <-
  function(k = 6,
           size = c(10, 5, 8, 2, 15, 50),
           rho = c(.7, .7, .5, .9, .85, .4),
           epsilon = .01,
           eidim = 2) 
  {
    ndim <- sum(size)# dim of correlation matrix
    bigcor <- matrix(0, nrow = ndim, ncol = ndim)
    
    ### generating the basic correlation matrix
    ss <- c(0,cumsum(size))
    for (i in 1:k) {
      top <- c(1, rho[i] ^ (seq(1:(size[i] - 1))))
      cor <- toeplitz(top)
      ti <- (ss[i]+1):ss[i+1]
      bigcor[ti,ti] <- cor
    }
 
    cor.nz <- noisecor(bigcor, epsilon, eidim)
    return(cor.nz)
    
  }


################################
##  Simulating the Hub Matrix (entries filled in using Toeplitz structure)
################################

# this function calculates a Toeplitz matrix with values descending
# from a user specified maximum to minimum.  The matrix has a
# block diagonal structure.  The size and base correlation for each
# block is user specified.


#' Simulating the Hub Matrix (entries filled in using Toeplitz structure)
#' @param k is the number of groups
#' @param  size is a vector of length k specifying the size of each group
#' @param  rho is a vector of length k specifying base correlation values
#' @param  epsilon <- (1-min(rho) - 0.75*min(tau) ) - .01 - maximum entry-wise random noise.
# Not used
# @param  tau_k = (max(rho_k) - min(rho_k) )/ (size_k -2) - step size 
#' @param  eidim is the space from which the noise is generated, the smaller the more noise
#' @param power = 2 makes the correlations stay high, = 0.5 makes the correlations descent rapidly
#' @export
#' @examples 
#' # Figure 2 in the paper
#' vw <- function(Im)
#' {
#' image(Im[, ncol(Im):1], zlim=c(-0.4, 1), col=rev(heat.colors(64)))
#' }
#' k<-3
#' sz <- c(100, 50, 80)
#' rho <- matrix(c(0.7, 0.7, 0.4, 0, 0, 0), nrow=3)
#' epsilon <- 0.23
#' hTC1 <- simcor.H(k=k, size=sz, rho=rho, power=1, epsilon=0.23, eidim=2)
#' vw(hTC1)
#' rho <- matrix(c(0.7, 0.7, 0.4, 0.5, 0.6, 0.2), nrow=3)
#' 
#' hTC2 <- simcor.H(k=k, size=sz, rho=rho, power=1, epsilon=0.29, eidim=2)
#' vw(hTC2)
#' 
#' hTC3 <- simcor.H(k=k, size=sz, rho=rho, power=1, epsilon=0.29, eidim=25)
#' vw(hTC3)
#' 
#' hTC4 <- simcor.H(k=k, size=sz, rho=rho, power=1, epsilon=0.1, eidim=2)
#' vw(hTC4)
#' hTC5 <- simcor.H(k=k, size=sz, rho=rho, power=1, epsilon=0.25, eidim=2)
#' vw(hTC5)
#' rho <- matrix(c(0.8, 0.75, 0.7, 0, 0, 0), nrow=3)
#' hTC6 <- simcor.H(k=k, size=sz, rho=rho, power=1, epsilon=0.19, eidim=2)
#' vw(hTC6)
simcor.H <- function(k = 6,
                          size = c(10, 5, 8, 7, 15, 50),
                          rho = rbind(c(.9, .7), c(.7, .7), c(.7, .2), c(.5, .3), c(.9, .85), c(.3, .2)),
                          power = 1,
                          epsilon = .08,
                          eidim = 2) 
{
  ndim <- sum(size)# dim of correlation matrix
  bigcor <- matrix(0, nrow=ndim, ncol = ndim)
  
  ### generating the basic correlation matrix
  
  ss <- c(0,cumsum(size))
  for (i in 1:(k)) {
    cor <- toeplitz(rho.func(rho[i, 1], rho[i, 2], power, size[i]))
    ti <- (ss[i]+1):ss[i+1]
    bigcor[ti,ti] <- cor
  }
  cor.nz <- noisecor(bigcor, epsilon, eidim)
  return(cor.nz)
}


# rho.func is needed for filling in the rest of the structure of the Hub
# correlation matrix
# r.max is the maximum user specified correlation
# r.min is the minimum user specified correlation
# power is the power at which the correlations descend
# p is the size of the correlation block


rho.func <- function(r.max, r.min, power, p) 
{
  pp <- 2:p
  rhovec <- r.max - ((pp - 2) / (p - 2)) ^ power * (r.max - r.min)
  rhovec <- c(1, rhovec)
  return(rhovec)
}


