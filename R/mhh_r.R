
# see: Martin, Hurn, and Harris, 2013, EMTSUtil.R

#--------------------------------------------------------------------
# Returns a matrix (or vector) stripped of the specified rows
#
#   Inputs: 
#             x  = input matrix (or vector) (n x k)
#             rb = first n1 rows to strip
#             re = last  n2 rows to strip
#
#--------------------------------------------------------------------
trimr <- function(x,rb,re) {
  x <- cbind(x)
  n <- nrow(x)
  if ((rb+re) >= n) {
    stop('Attempting to trim too much')
  }
  z <- x[(rb+1):(n-re),,drop=F]
  return(z)  
}

# see: Martin, Hurn, and Harris, 2013, coint_reg.R

#----------------------------------------------------------------------------
#  Long-run variance
#----------------------------------------------------------------------------
lrvars <- function(z) {
  
  t <- nrow(z)
  k <- ncol(z)
  # p <- floor(8*(t/100)^(2/9)) # preliminary P based on Bartlett (see: MHH, 2013, table 9.1, p. 326)
  p <- floor(4*(t/100)^(2/9)) # preliminary P based on Bartlett (see: MHH, 2013, table 9.1, p. 326)
  
  J0 <- t(z) %*% z 
  J1 <- 0
  
  for (j in seq(p)) {
    Gj <- t(trimr(z,j,0)) %*% trimr(z,0,j)
    J0 <- J0 + Gj + t(Gj)
    J1 <- J1 + 2*j*Gj
  }
  
  i <- ones(k,1)
  v0 <- t(i) %*% J0 %*% i 
  v1 <- t(i) %*% J1 %*% i
  p <- min( rbind( floor(1.1447*((v1/v0)^2*t)^(1/3)), p ) ) # updated P based on Bartlett (see: MHH, 2013, table 9.1, p. 326) 
  
  Omega <- t(z) %*% z/t 
  Delta <- t(z) %*% z/t
  
  weights <- matrix(0, nrow = p)
  for (j in seq(p)) {
    Gj    <- t(trimr(z,j,0)) %*% trimr(z,0,j)/t
    Omega <- Omega + (1-j/(p+1))*(Gj + t(Gj)) # weighting scheme wi = (1-j/(p+1)) based on Bartlett (see: MHH, 2013, table 9.1, p. 326) 
    Delta <- Delta + (1-j/(p+1))*Gj
    weights[j] <- (1-j/(p+1))
  }
  
  return(list(Delta=Delta, Omega=Omega, weights=weights))
  
}

lrvars_parzen <- function(z) {
  
  t <- nrow(z)
  k <- ncol(z)
  p <- floor(8*(t/100)^(4/25)) # preliminary P based on Parzen (see: MHH, 2013, table 9.1, p. 326)
  
  J0 <- t(z) %*% z 
  J1 <- 0
  J2 <- 0 # for Parzen
  
  for (j in seq(p)) {
    Gj <- t(trimr(z,j,0)) %*% trimr(z,0,j)
    J0 <- J0 + Gj + t(Gj)
    J1 <- J1 + 2*j*Gj
    J2 <- J2 + 2*j^2*Gj # for Parzen (see: MHH, 2013, equation 9.27, p. 325)
  }
  
  i <- ones(k,1)
  v0 <- t(i) %*% J0 %*% i 
  v1 <- t(i) %*% J1 %*% i
  v2 <- t(i) %*% J2 %*% i
  p <- min( rbind( floor(2.6614*((v2/v0)^2*t)^(1/5)), p ) ) # updated P based on Parzen (see: MHH, 2013, table 9.1, p. 326) 
  
  Omega <- t(z) %*% z/t 
  Delta <- t(z) %*% z/t
  
  weights <- matrix(0, nrow = p)
  for (j in seq(p)) {
    Gj    <- t(trimr(z,j,0)) %*% trimr(z,0,j)/t
    # Omega <- Omega + (1-j/(p+1))*(Gj + t(Gj)) # weighting scheme wi = (1-j/(p+1)) based on Bartlett (see: MHH, 2013, table 9.1, p. 326) 
    # Delta <- Delta + (1-j/(p+1))*Gj
    if (j <= (p+1)/2) {
      Omega <- Omega + (1-6*(j/(p+1))^2-6*(j/(p+1))^3)*(Gj + t(Gj)) # weighting scheme wi based on Parzen (see: MHH, 2013, table 9.1, p. 326) 
      Delta <- Delta + (1-6*(j/(p+1))^2-6*(j/(p+1))^3)*Gj
      weights[j] <- (1-6*(j/(p+1))^2-6*(j/(p+1))^3)
    } else if (j > (p+1)/2) {
      Omega <- Omega + (2*(1-j/(p+1))^3)*(Gj + t(Gj)) # weighting scheme wi based on Parzen (see: MHH, 2013, table 9.1, p. 326)
      Delta <- Delta + (2*(1-j/(p+1))^3)*Gj
      weights[j] <- (2*(1-j/(p+1))^3)
    }
  }
  
  return(list(Delta=Delta, Omega=Omega, weights=weights))
  
}