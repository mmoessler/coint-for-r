
# clear workspace
rm(list=ls())

# application based on simulation

source("./R/gauss_r.R")
source("./R/mhh_r.R")

source("./R/app_coint_2_1_base.R")
source("./R/app_coint_2_1_kernels.R")
source("./R/app_coint_2_1_cregrs.R")

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

# estimation

FM_OLS.res.01 <- FM_OLS(y = y.mat[,c(1,2),drop=F], x = y.mat[,c(3),drop=F],
                        d = NA, l = NULL, ker_fun = "parzen", aband = 1, filter = 0, NoDet = 1)

FM_OLS.res.01$beta
FM_OLS.res.01$var
head(FM_OLS.res.01$checks$uxk)




# check estimation results

FM_OLS.res.01$checks$bige

lrvars.res.01 <- lrvars(z = FM_OLS.res.01$checks$uxk)
lrvars.res.02 <- lrvars_parzen(z = FM_OLS.res.01$checks$uxk)

FM_OLS.res.01$checks$bige
lrvars.res.01$Omega
lrvars.res.02$Omega
# -> different

# steps: lrvars ----

z <- FM_OLS.res.01$checks$uxk

t <- nrow(z)
k <- ncol(z)
p <- floor(8*(t/100)^(2/9))

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
p <- min( rbind( floor(1.1447*((v1/v0)^2*t)^(1/3)), p ) )

Omega <- t(z) %*% z/t 
Delta <- t(z) %*% z/t

for (j in seq(p)) {
  Gj    <- t(trimr(z,j,0)) %*% trimr(z,0,j)/t
  Omega <- Omega + (1-j/(p+1))*(Gj + t(Gj))
  Delta <- Delta + (1-j/(p+1))*Gj
}

Omega
Delta

# steps: lrvar -> kacf -> kernel -> parzen ----

e <- FM_OLS.res.01$checks$uxk
v <- FM_OLS.res.01$l

# 1) bige <- lrvar(tmp, l, aband = aband, filter = filter) ----

# 2) io <- kacf(e, v, aband = aband, ker_fun = ker_fun) ----

e <- e
v <- v
aband <- 1
filter <- 0
ker_fun <- c("parzen")

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

# 3) ret <- kernel(e, v, ker_fun = ker_fun) ----

e <- e
v <- v

# 4) ret <- parzen(uv, k) ----

uv <- e
k <- v

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

io <- ret

# MM added (for potentially io=0)
io <- matrix(io, nrow = ncol(e), ncol = ncol(e))
s  <- (t(e) %*% e)/nrow(e)
lr <- s + io + t(io)

lr
FM_OLS.res.01$checks$bige
# -> same!
