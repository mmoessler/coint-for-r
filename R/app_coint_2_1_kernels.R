
# see: app_coint_2_1/src/kernels.src

kernel <- function(uv, l, ker_fun = c("parzen"), check = FALSE) {
  
  # checks
  if (check == TRUE) {
    uv <- e
    l <- v
    ker_fun <- ker_fun
  }
  
  
  
  if (nrow(uv) <= 1) {
    print("ERROR: data must have more than one observation to estimate autocorrelations")
  }
  if (ker_fun == "qs") {
    ret <- -9999
  } else if (ker_fun == "parzen") {
    ret <- parzen(uv, l)
  } else if (ker_fun == "fejer") {
    ret <- -9999
  } else if (ker_fun == "tukham") {
    ret <- -9999
  }
  
  return(ret)
  
}

# /*
# ** Computes the Parzen window
# ** Brillinger (1981). p55
# */
parzen <- function(uv, k, check = FALSE) {
  
  # checks
  if (check == TRUE) {
    uv <- uv
    k <- l
  }
  
  
  
  if (k > nrow(uv)) {
    k <- nrow(uv) - 1
  }
  if (k >= 1) {
    weights <- zeros(k,1)
  }
  i <- 1
  a <- 0
  while (i < k/2) {
    m  <- 1.0 - 6*((i/(k+1))^2) + 6*((abs(i)/(k+1))^3)
    t1 <- detrend(trimr(uv,i,0),0)
    t2 <- detrend(trimr(lagn(uv,i),i,0),0)
    a  <- a + m * (t(t1) %*% t2)
    weights[i] <- m
    i <- i + 1
  }
  while (i < k) {
    m <- 2*(1 - (abs(i)/(k+1)))^3
    t1 <- detrend(trimr(uv,i,0),0)
    t2 <- detrend(trimr(lagn(uv,i),i,0),0)
    a  <- a + m * (t(t1) %*% t2)
    weights[i] <- m
    i <- i + 1
  }
  ret <- a/nrow(uv)

  parzen_weights <<- weights # store in global environment for checks
  
  return(ret)
  
}
