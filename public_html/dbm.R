# dbm.R
n <- 1000
t <- 1
dt <- t/n
bm <- c(0, cumsum(rnorm(n,0,sqrt(dt))))
steps <- seq(0,t,length=n)
W <- diff(bm)/dt
plot(steps,W,type="l")
