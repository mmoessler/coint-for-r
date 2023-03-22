
# see: app_coint_2_1/src/cregrs.src

source("./R/gauss_r.R")
source("./R/mhh_r.R")

source("./R/app_coint_2_1_base.R")
source("./R/app_coint_2_1_kernels.R")

FM_OLS <- function(y, x, d, l, ker_fun = c("parzen"), aband = 1, filter = 0, NoDet = 0, check = FALSE){
  
  # checks
  if (check == TRUE) {
    y <- y.mat[,c(1,2),drop=F]
    x <- y.mat[,c(3),drop=F]
    d <- NA
    l <- NULL
    ker_fun <- "parzen"
    aband <- 1
    filter <- 0
    NoDet <- 1
  }
  
  
  
  d <- as.matrix(d)
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  t <- rows(y) - 1
  m <- cols(x)
  n <- cols(y)
  
  if (NoDet != 1) {
    xd <- x - (d %*% inv(t(d) %*% d) %*% (t(d) %*% x))
    z  <- cbind(x,d)
  } else {
    xd <- x
    z <- x
  }
  
  ixx <- inv(t(z) %*% z)
  xy <- t(z) %*% y
  beta <- ixx %*% xy
  
  u <- y - (z %*% beta)
  e <- cbind(trimr(u, 1, 0), diff(xd, 1))
  
  del <- delta(e, l, aband = aband, filter = filter) # one-sided long-run variance matrix
  del <- trimr(t(trimr(t(del), n, 0)), 0, m) # lower left block
  sig <- lrvar(e, l, aband = aband, filter = filter) # long-run variance matrix
  sigxx <- trimr(t(trimr(sig, n, 0)), n, 0) # lower right block
  
  delxx <- delta(diff(xd, 1), l, aband = aband, filter = filter)
  
  true <- del %*% inv(sigxx)
  
  ys <- trimr(y,1,0) - diff(xd,1) %*% t(true)
  dels <- t(del) - (delxx %*% t(true))
  
  if ( NoDet != 1 ) {
    dels <- rbind(dels, zeros(cols(d), n))
  }
  
  xk <- trimr(z, 1, 0)
  ixx <- inv(t(xk) %*% xk)
  beta <- ixx %*% ((t(xk) %*% ys) - (rows(xk) * dels)) # /* Mutivariate FM */
  
  # /* Okay, compute the co-variance matrix.... */
  # tmp <- trimr(u[,1,drop=F], 1, 0) * xk
  tmp <- e_by_e_mul(trimr(u[,1,drop=F], 1, 0), xk)
  j <- 1
  while ( j < cols(u) ) {
    j <- j + 1
    # tmp <- cbind(tmp, (trimr(u[,j,drop=F], 1, 0) * xk))
    tmp <- cbind(tmp,(e_by_e_mul(trimr(u[,j,drop=F],1,0),xk)))
  }
  uxk <- tmp # for checks
  
  bige <- lrvar(tmp, l, aband = aband, filter = filter, ker_fun = ker_fun) # meat of the sandwich (MK x MK)
  tmp <- eye(cols(u)) %x% (sqrt(rows(xk))*ixx) # bread of the sandwich (MK x MK)
  
  var <- tmp %*% bige %*% tmp
  
  ret.lis <- list(beta = beta, var = var, checks = list(uxk = uxk, bige = bige, l = l, u = u, xk = xk))
  
  return(ret.lis)
  
}
