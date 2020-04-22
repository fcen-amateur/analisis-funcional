library(fda)
library(ggplot2)
# Chapter 1
spline.basis=create.bspline.basis(rangeval=c(0,10000), nbasis=20)
fourier.basis=create.fourier.basis(rangeval=c(0,10000), nbasis=13)

# fig 1.3
plot(spline.basis)
plot(fourier.basis)

# fig 1.4 - proceso de Wiener

Wiener=cumsum(rnorm(10000))/100 # random walk on [0,K], K=10ˆ4
plot.ts(Wiener, xlab="", ylab="")
B25.basis=create.bspline.basis(rangeval=c(0,10000), nbasis=25)
Wiener.fd=smooth.basis(y=Wiener, fdParobj=B25.basis)
Wiener.fdf=smooth.basis(y=Wiener, fdParobj=fourier.basis)
lines(Wiener.fd, lwd=2, lty=2)
lines(Wiener.fdf, lwd=1)
