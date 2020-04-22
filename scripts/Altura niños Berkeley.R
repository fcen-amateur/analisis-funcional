rm(list=ls())
library(fda)
data(growth)
demo(growth)

par(mfrow=c(1,1))

matplot(growth$age, growth$hgtm,type="b", lty=1, pch=16, lwd=2,ylab="Altura",xlab="Años")

matplot(agefine, velmhat,type="l", lty=1, pch=16, lwd=2,ylab="Velocidad",xlab="Años")


matplot(agefine, accmhat,type="l", lty=1, pch=16, lwd=2,ylab="Aceleracion",xlab="Años")



matplot(growth$age, growth$hgtf,type="b", lty=1, pch=16, lwd=2,ylab="Altura",xlab="Años")


matplot(agefine, velfhat,type="l", lty=1, pch=16, lwd=2,ylab="Velocidad",xlab="Años")


matplot(agefine, accfhat,type="l", lty=1, pch=16, lwd=2,ylab="Aceleracion",xlab="Años")



algunos= 1:10


matplot(growth$age, growth$hgtm[,algunos],type="l", lty=1, pch=16, lwd=2,ylab="Altura",xlab="Años")
matplot(growth$age, growth$hgtm[,algunos],type="p", lty=1, pch=16, lwd=2,add=T, cex=1.2)



matplot(growth$age, growth$hgtf[,algunos],type="l", lty=1, pch=16, lwd=2,ylab="Altura",xlab="Años")
matplot(growth$age, growth$hgtf[,algunos],type="p", lty=1, pch=16, lwd=2,add=T, cex=1.2)


matplot(agefine, velmhat[,algunos],type="l", lty=1, pch=16, lwd=2,ylab="Velocidad",xlab="Años")
###########################
# BOXPLOT FUNCIONAL
#################################


aa=fbplot(growth$hgtm, x = growth$age, ylab="Altura",xlab="Años",xlim=c(0,18), ylim=c( min(growth$hgtm),  max(growth$hgtm)))
aa=fbplot(growth$hgtf, x = growth$age, ylab="Altura",xlab="Años",xlim=c(0,18), ylim=c( min(growth$hgtf),  max(growth$hgtf)))

aa=fbplot(growth$hgtf, x = growth$age, prob=c(0.7,0.5,0.2), color=c("bisque1", "chocolate1", "chocolate4"), ylab="Altura",xlab="Años",xlim=c(0,18), ylim=c( min(growth$hgtm),  max(growth$hgtm)))

nombre=paste("fbplot-altura-mujeres-probs.pdf",sep="")
	pdf(nombre, bg='transparent')

aa=fbplot(growth$hgtf, x = growth$age, prob=c(0.9,0.5,0.2), color=c("bisque2", "magenta", "chocolate4"), ylab="Altura",xlab="Años",xlim=c(0,18), ylim=c( min(growth$hgtf),  max(growth$hgtf)))
dev.off()


nombre=paste("fbplot-altura-hombres-probs.pdf",sep="")
	pdf(nombre, bg='transparent')

aa=fbplot(growth$hgtm, x = growth$age, prob=c(0.9,0.5,0.2), color=c("bisque2", "magenta", "chocolate4"), ylab="Altura",xlab="Años",xlim=c(0,18), ylim=c( min(growth$hgtm),  max(growth$hgtm)))
dev.off()




 

aa$out


##
## 1. smooth the growth data for the Berkeley boys
##
# Specify smoothing weight
lambda.gr2.3 <- .03
# Specify what to smooth, namely the rate of change of curvature
Lfdobj.growth <- 2
# Set up a B-spline basis for smoothing the discrete data
nage <- length(growth$age)
norder.growth <- 6
nbasis.growth <- nage + norder.growth - 2
rng.growth <- range(growth$age)
wbasis.growth <- create.bspline.basis(rangeval=rng.growth,
nbasis=nbasis.growth, norder=norder.growth,
breaks=growth$age)
# Smooth the data
# in afda-ch06.R, and register to individual smooths:
cvec0.growth <- matrix(0,nbasis.growth,1)
Wfd0.growth <- fd(cvec0.growth, wbasis.growth)
growfdPar2.3 <- fdPar(Wfd0.growth, Lfdobj.growth, lambda.gr2.3)
hgtmfd.all <- with(growth, smooth.basis(age, hgtm, growfdPar2.3)$fd)
# Register the growth velocity rather than the
# growth curves directly
smBv <- deriv(hgtmfd.all, 1)
##
## 2. Register the first 2 Berkeley boys using the default basis
## for the warping function
##
# register.fd takes time, so use only 2 curves as an illustration
# to minimize compute time in these examples
nBoys <- 2
# Define the
Define the target function as the mean of the first nBoys records
smBv0 = mean.fd(smBv[1:nBoys])
# Register these curves. The default choice for the functional
# parameter object WfdParObj is used.
smB.reg.0 <- register.fd(smBv0, smBv[1:nBoys])
# plot each curve. Click on the R Graphics window to show each plot.
# The left panel contains:
# -- the unregistered curve (dashed blue line)
# -- the target function (dashed red line)
# -- the registered curve (solid blue line)
# The right panel contains:
# -- the warping function h(t)
# -- the linear function corresponding to no warping
plotreg.fd(smB.reg.0)
# Notice that all the warping functions all have simple shapes
# due to the use of the simplest possible basis
##
## 3. Define a more flexible basis for the warping functions
##
if(!CRAN()){
Wnbasis <- 4
Wbasis <- create.bspline.basis(rng.growth, Wnbasis)
Wfd0 <- fd(matrix(0,Wnbasis,1),Wbasis)
# set up the functional parameter object using only
# a light amount smoothing
WfdParobj <- fdPar(Wfd0, Lfdobj=2, lambda=0.01)
# register the curves
smB.reg.1 <- register.fd(smBv0, smBv[1:nBoys], WfdParobj)
plotreg.fd(smB.reg.1)
## 4. Change the target to the mean of the registered functions ...
## this should provide a better target for registration
##
smBv1 <- mean.fd(smB.reg.1$regfd)
# plot the old and the new targets
par(mfrow=c(1,1),ask=FALSE)
plot(smBv1)
lines(smBv0, lty=2)
# Notice how the new target (solid line) has sharper features and
# a stronger pubertal growth spurt relative to the old target
# (dashed line). Now register to the new target
smB.reg.2 <- register.fd(smBv1, smBv[1:nBoys], WfdParobj)
plotreg.fd(smB.reg.2)
# Plot the mean of these curves as well as the first and second targets
par(mfrow=c(1,1),ask=FALSE)
plot(mean.fd(smB.reg.2$regfd))
lines(smBv0, lty=2)
lines(smBv1, lty=3)
# Notice that there is almost no improvement over the age of the
# pubertal growth spurt, but some further detail added in the
# pre-pubertal region. Now register the previously registered
# functions to the new target.
smB.reg.3 <- register.fd(smBv1, smB.reg.1$regfd, WfdParobj)
plotreg.fd(smB.reg.3)
# Notice that the warping functions only deviate from the straight line
# over the pre-pubertal region, and that there are some small adjustments
# to the registered curves as well over the pre-pubertal region.
##}
} 

par(mfrow=c(1,1))
plot.fd( smBv[1:10],lwd=2, xlab="Edad", ylab="Velocidad", ylim=c(0,17))
smB.reg.0 <- register.fd(smBv0, smBv[1:10])
plot(smB.reg.0$regfd, lwd=2, xlab="Edad", ylab="Vel-warp", ylim=c(0,17))

lines(agefine, velmhat[,1],type="l", lty=2,  lwd=2,col="gray")


smB.reg.0 <- register.fd(smBv0, smBv[2])
plot(smB.reg.0$regfd, add=T,col="green4")

lines(agefine, velmhat[,2],type="l", lty=1,  lwd=2,col="green")





