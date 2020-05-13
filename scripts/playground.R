library(fda)
library(tidyverse)
# Chapter 1

# fig 1.3

# fig 1.4 - proceso de Wiener
b <- 100
nbasis <- 1000
Exp <- 1 / seq.int(100)
Wiener=cumsum(rnorm(b))/b # random walk on [0,K], K=10ˆ4
spline.basis=create.bspline.basis(rangeval=c(b), nbasis)
fourier.basis=create.fourier.basis(rangeval=c(0,b), nbasis)
plot.ts(x, xlab="", ylab="")
x.fd=smooth.basis(y=x, fdParobj=spline.basis)
x.fdf=smooth.basis(y=x, fdParobj=fourier.basis)
lines(x.fd, lwd=2, lty=2)
lines(x.fdf, lwd=1)

# Problema 1.1
COLS <- 2
ROWS <- 1
NBASIS <- 15
NORDER <- 4
df <- fda::pinch
J <- dim(df)[1]
N <- dim(df)[2]
c(min(df), max(df))
apply(df, COLS, max)
base <- create.bspline.basis(rangeval = J, nbasis = NBASIS, norder = NORDER)
df.suave <- smooth.basis(y = df, fdParobj =  base)

dev.off()
plot(df.suave$argvals, df.suave$y[,1], 'l')
for (i in seq.int(2, N)) {
  lines(df.suave$argvals, df.suave$y[,i], 'l', lwd = 0.5)
}
lines(df.suave)
mu <- apply(df, ROWS, mean)
sigma <- apply(df, ROWS, sd)
lines(mu, lwd = 1)
lines(sigma, lwd = 1)
plot.ts(df[,1:10])
