
# see: app_coint_2_1/src/base.src

diff <- function(y, k) {
  
  y <- as.matrix(y)
  t <- nrow(y)
  n <- ncol(y)
  
  if (k == 0) {
    return(y)
  }
  
  y <- trimr(y,k,0) - trimr(lagn(y,k),k,0)
  
  return(y)
  
}

xahat <- function(e) {
  
  t1 <- trimr(e,1,0)
  t2 <- trimr(lagn(e,1),1,0)
  a  <- inv(t(t2) %*% t2) %*% (t(t2) %*% t1)
  r  <- t1 - (t2 %*% a)
  
  ret.lis <- list(a=a,r=r)
  
  return(ret.lis)
  
}

detrend <- function(data, p){
  
  if (p == -1) {
    return(data)
  }
  nobs <- nrow(data)
  u    <- ones(nobs,1)
  if (p > 0) {
    timep <- zeros(nobs,p)
    t     <- seqa(1,1,nobs)/nobs
    m     <- 1
    while (m <= p) {
      timep[,m] <- t^m
      m <- m + 1
    }
    xmat <- cbind(u,timep)
  } else {
    xmat <- u
  }
  invx  <- inv(t(xmat) %*% xmat)
  beta  <- invx %*% (t(xmat) %*% data)
  resid <- data - xmat %*% beta
  
  return(resid)
  
}

# automatic bandwidth selection
kacf <- function(e, v, aband = 1, ker_fun = c("parzen"), check = FALSE) {
  
  # checks
  if (check == TRUE) {
    e <- e
    v <- v
    aband <- aband
    ker_fun <- ker_fun
    
  }
  
  
  
  if (aband == 1) {
    eb <- trimr(lagn(e, 1), 1, 0)
    ef <- trimr(e, 1, 0)
    ae <- sumc(eb*ef)/sumc(eb^2)
    # ee <- ef - eb*(t(ae))
    ee <- ef - eb*matrix(t(ae),nrow=nrow(eb),ncol=ncol(t(ae)),byrow=TRUE) # use e.g., e_by_e_mul 
    se <- meanc(ee^2)
    ad <- sumc((se/((1-ae)^2))^2)
    a1 <- 4*sumc((ae*se/(((1-ae)^3)*(1+ae)))^2)/ad
    a2 <- 4*sumc((ae*se/((1-ae)^4))^2)/ad
    nobs <- nrow(e)
    if (ker_fun == "qs") {
      v <- 1.3221*((a2*nobs)^.2)-1
    } else if (ker_fun == "parzen") {
      v <- 2.6614*((a2*nobs)^.2)-1
    } else if (ker_fun == "fejer") {
      v <- 1.1447*((a1*nobs)^.333)-1
    } else if (ker_fun == "tukham") {
      v <- 1.7462*((a2*nobs)^.2)-1
    }
  } 
  
  ret <- kernel(e, v, ker_fun = ker_fun)
  
  return(ret)
  
}

# long-run variance matrix
lrvar <- function(e, v, aband = 1, filter = 0, ker_fun = c("parzen"), check = FALSE) {
  
  # checks
  if (check == TRUE) {
    e <- tmp
    v <- l
    aband <- aband
    filter <- filter
    ker_fun <- ker_fun
  }
  
  
  
  if (filter == 1) { # use AR(1) filter
    xahat <- xahat(e)
    a <- xahat$a
    new_e <- xahat$r
    tmp <- inv(eye(ncol(e)) - a)
    io <- kacf(new_e, v, aband = aband, ker_fun = ker_fun)
    # MM added (for potentially io=0)
    io <- matrix(io, nrow = ncol(e), ncol = ncol(e))
    s  <- (t(new_e) %*% new_e)/nrow(new_e)
    lr <- t(tmp) %*% (s + io + t(io)) %*% tmp
  } else {
    io <- kacf(e, v, aband = aband, ker_fun = ker_fun)
    # MM added (for potentially io=0)
    io <- matrix(io, nrow = ncol(e), ncol = ncol(e))
    s  <- (t(e) %*% e)/nrow(e)
    lr <- s + io + t(io)
  }
  
  ret <- lr
  
  return(ret)
  
}

# one-sided long-run variance matrix
delta <- function(e, v, aband = 1, filter = 0) {
  
  if (filter == 1) {
    xahat <- xahat(e)
    a <- xahat$a
    new_e <- xahat$r
    tmp <- inv(eye(cols(e)) - a)
    io <- kacf(new_e, v, aband = aband)
    s  <- (t(new_e) %*% new_e)/rows(new_e)
    su <- (t(e) %*% e)/(rows(e)) 
    # lr <- su + t(tmp) %*% io %*% tmp + (t(tmp) %*% (t(a)) %*% su)
    lr <- su + (t(tmp) * io) %*% tmp + (t(tmp) %*% (t(a)) %*% su)
  } else {
    io <- kacf(e,v, aband = aband)
    lr <- ((t(e) %*% e)/rows(e)) + (io)
  }
  
  ret <- t(lr)
  
  return(ret)
  
}
