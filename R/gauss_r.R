
# see: gauss/src/...

rows <- function(x) {
  
  nrow(x)
  
}

cols <- function(x) {
  
  ncol(x)
  
}

zeros <- function(n, m){
  
  matrix(0, nrow = n, ncol = m)
  
}

ones <- function(n,m){
  
  matrix(1,nrow = n, ncol = m)
  
}

eye <- function(n) {
  
  diag(n)
  
}

inv <- function(x) {
  
  solve(x)
  
}

lagn <- function(y, k) {
  
  y <- as.matrix(y)
  t <- nrow(y)
  n <- ncol(y)
  
  if (n == 1) {
    y <- rbind(matrix(rep(NA,k),nrow=k,ncol=1),trimr(y,0,k))
  } else {
    y <- rbind(matrix(rep(NA,k*n),nrow=k,ncol=n),trimr(y,0,k))
  }
  
  return(y)
  
}

sumc <- function(x) {
  
  x <- cbind(x)
  colSums(x)
  
}

meanc <- function(x) {
  
  colMeans(x)
  
}



e_by_e_mul <- function(x,y){
  
  # emulate Gauss element by element multiplication
  if (ncol(x) > ncol(y)){
    z <- x*matrix(y,nrow=nrow(x),ncol=ncol(x),byrow=FALSE)
  } else {
    z <- y*matrix(x,nrow=nrow(y),ncol=ncol(y),byrow=FALSE)
  }
  
  # what about same column dimension and different row dimension?
  
  z
  
}
