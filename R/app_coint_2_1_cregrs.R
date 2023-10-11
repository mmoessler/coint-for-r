
# see: app_coint_2_1/src/cregrs.src

source("./R/gauss_r.R")
source("./R/mhh_r.R")

source("./R/app_coint_2_1_base.R")
source("./R/app_coint_2_1_kernels.R")

FM_OLS <- function(y, x, d, l, ker_fun = c("parzen"), aband = 1, filter = 0, NoDet = 0, check = FALSE){
  
  # checks
  if (check == TRUE) {
    
    # see: https://www.pfaffikus.de/books/spex2/files/Rcode-4-3.R
    set.seed(12345)
    e1 <- rnorm(250, 0, 0.5)
    e2 <- rnorm(250, 0, 0.5)
    e3 <- rnorm(250, 0, 0.5)
    u1.ar1 <- arima.sim(model = list(ar = 0.75),
                        innov = e1, n = 250)
    u2.ar1 <- arima.sim(model = list(ar = 0.3),
                        innov = e2, n = 250)
    y3 <- cumsum(e3)
    y1 <- 0.8 * y3 + u1.ar1
    y2 <- -0.3 * y3 + u2.ar1
    y.mat <- data.frame(y1, y2, y3)
    
    y <- y.mat[,c(1,2),drop=F]
    x <- y.mat[,c(3),drop=F]
    d <- NA
    l <- NULL
    ker_fun <- "qs"
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
  
  del <- delta(e = e, v = l, aband = aband, filter = filter, ker_fun = ker_fun) # one-sided long-run variance matrix
  del <- trimr(t(trimr(t(del), n, 0)), 0, m) # lower left block
  sig <- lrvar(e = e, v = l, aband = aband, filter = filter, ker_fun = ker_fun) # long-run variance matrix
  sigxx <- trimr(t(trimr(sig, n, 0)), n, 0) # lower right block
  
  delxx <- delta(diff(xd, 1), l, aband = aband, filter = filter, ker_fun = ker_fun)
  
  true <- del %*% inv(sigxx)
  
  ys <- trimr(y, 1, 0) - diff(xd, 1) %*% t(true)
  dels <- t(del) - (delxx %*% t(true))
  
  if ( NoDet != 1 ) {
    dels <- rbind(dels, zeros(cols(d), n))
  }
  
  xk <- trimr(z, 1, 0)
  ixx <- inv(t(xk) %*% xk)
  beta <- ixx %*% ((t(xk) %*% ys) - (rows(xk) * dels)) # /* Mutivariate FM */
  
  # /* Okay, compute the co-variance matrix.... */
  # tmp <- trimr(u[,1,drop=F], 1, 0) * xk # use e_by_e_mul() instead
  tmp <- e_by_e_mul(trimr(u[,1,drop=F], 1, 0), xk)
  j <- 1
  while ( j < cols(u) ) {
    j <- j + 1
    # tmp <- cbind(tmp, (trimr(u[,j,drop=F], 1, 0) * xk)) # use e_by_e_mul() instead
    tmp <- cbind(tmp, (e_by_e_mul(trimr(u[,j,drop=F], 1, 0), xk)))
  }
  uxk <- tmp # for checks
  
  bige <- lrvar(tmp, l, aband = aband, filter = filter, ker_fun = ker_fun) # meat of the sandwich (MK x MK)
  tmp <- eye(cols(u)) %x% (sqrt(rows(xk))*ixx) # bread of the sandwich (MK x MK)
  
  var <- tmp %*% bige %*% tmp
  
  ret.lis <- list(beta = beta, var = var, checks = list(uxk = uxk, bige = bige, l = l, u = u, xk = xk))
  
  return(ret.lis)
  
}

# combination with cointreg package

# packages we will use in the course
pac <- c("cointReg")
# install and/or load packages
checkpac <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
  }
  require(x, character.only = TRUE)
}
# check if packages are install yet
suppressWarnings(sapply(pac, checkpac))

FM_OLS_cointreg <- function(y, x, d, l, ker_fun = c("parzen"), aband = 1, filter = 0, NoDet = 0, check = FALSE,
                            cointreg = list(kernel = c("qs"), bandwidth = c("and"), demeaning = FALSE, check = FALSE)){
  
  # checks
  if (check == TRUE) {
    
    # see: https://www.pfaffikus.de/books/spex2/files/Rcode-4-3.R
    set.seed(12345)
    e1 <- rnorm(250, 0, 0.5)
    e2 <- rnorm(250, 0, 0.5)
    e3 <- rnorm(250, 0, 0.5)
    u1.ar1 <- arima.sim(model = list(ar = 0.75),
                        innov = e1, n = 250)
    u2.ar1 <- arima.sim(model = list(ar = 0.3),
                        innov = e2, n = 250)
    y3 <- cumsum(e3)
    y1 <- 0.8 * y3 + u1.ar1
    y2 <- -0.3 * y3 + u2.ar1
    y.mat <- data.frame(y1, y2, y3)
    
    y <- y.mat[,c(1,2),drop=F]
    x <- y.mat[,c(3),drop=F]
    d <- NA
    l <- NULL
    ker_fun <- "qs"
    aband <- 1
    filter <- 0
    NoDet <- 1
    
    cointreg <- list(kernel = c("qs"), bandwidth = c("and"), demeaning = FALSE, check = FALSE)
    
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
  
  if(!is.null(cointreg)) {
    
    source("https://raw.githubusercontent.com/cran/cointReg/master/R/longrun-var.R")
    source("https://raw.githubusercontent.com/cran/cointReg/master/R/bandwidth.R")
    source("https://raw.githubusercontent.com/cran/cointReg/master/R/zzz.R")
    
    y.k <- ncol(y)
    x.k <- ncol(x)
    
    if (!is.numeric(cointreg$bandwidth)) {
      bw <- getBandwidth(e, bandwidth = cointreg$bandwidth, kernel = cointreg$kernel, check = FALSE)
      band <- switch(cointreg$bandwidth, and = "Andrews", nw = "Newey-West")
    } else {
      bw <- bandwidth
      band <- "set by user"
    }
    
    lr.var <- getLongRunVar(e, kernel = cointreg$kernel, bandwidth = bw, demeaning = cointreg$demeaning, check = FALSE)
    
    tmp <- lapply(lr.var, function(x) {
      out <- list()
      out[["all"]] <- x
      out[["uu"]] <- x[1:y.k, 1:y.k, drop = FALSE]
      out[["uv"]] <- x[1:y.k, (y.k + 1):(y.k + x.k), drop = FALSE]
      out[["vu"]] <- x[(y.k + 1):(y.k + x.k), 1:y.k, drop = FALSE]
      out[["vv"]] <- x[(y.k + 1):(y.k + x.k), (y.k + 1):(y.k + x.k), drop = FALSE]
      return(out)
    })
    
    Omega <- tmp[[1]]
    Delta <- tmp[[2]]
    
    sigxx <- Omega$vv
    del <- Delta$uv
    delxx <- Delta$vv
    
  } else {
    
    del <- delta(e = e, v = l, aband = aband, filter = filter, ker_fun = ker_fun) # one-sided long-run variance matrix
    del <- trimr(t(trimr(t(del), n, 0)), 0, m) # upper right block (Lambda_12)
    sig <- lrvar(e = e, v = l, aband = aband, filter = filter, ker_fun = ker_fun) # long-run variance matrix
    sigxx <- trimr(t(trimr(sig, n, 0)), n, 0) # lower right block (Omega_22)
    
    delxx <- delta(diff(xd, 1), v = l, aband = aband, filter = filter, ker_fun = ker_fun) # one-sided long-run variance matrix of v
    
  }

  true <- del %*% inv(sigxx)
  
  ys <- trimr(y,1,0) - diff(xd, 1) %*% t(true)
  dels <- t(del) - (delxx %*% t(true))
  
  if ( NoDet != 1 ) {
    dels <- rbind(dels, zeros(cols(d), n))
  }
  
  xk <- trimr(z, 1, 0)
  ixx <- inv(t(xk) %*% xk)
  beta <- ixx %*% ((t(xk) %*% ys) - (rows(xk) * dels)) # /* Mutivariate FM */
  
  # for check
  cointRegFM(x = x, y = y, deter = NULL, kernel = "qs", bandwidth = "and", demeaning = FALSE, check = TRUE)$theta
  beta
  
  # /* Okay, compute the co-variance matrix.... */
  # tmp <- trimr(u[,1,drop=F], 1, 0) * xk # use e_by_e_mul() instead
  tmp <- e_by_e_mul(trimr(u[,1,drop=F], 1, 0), xk)
  j <- 1
  while ( j < cols(u) ) {
    j <- j + 1
    # tmp <- cbind(tmp, (trimr(u[,j,drop=F], 1, 0) * xk)) # use e_by_e_mul() instead
    tmp <- cbind(tmp, (e_by_e_mul(trimr(u[,j,drop=F], 1, 0), xk)))
  }
  uxk <- tmp # for checks
  
  bige <- lrvar(tmp, l, aband = aband, filter = filter, ker_fun = ker_fun) # meat of the sandwich (MK x MK)
  tmp <- eye(cols(u)) %x% (sqrt(rows(xk))*ixx) # bread of the sandwich (MK x MK)
  
  var <- tmp %*% bige %*% tmp
  
  ret.lis <- list(beta = beta, var = var, checks = list(uxk = uxk, bige = bige, l = l, u = u, xk = xk))
  
  return(ret.lis)
  
}
