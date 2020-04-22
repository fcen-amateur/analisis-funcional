library(fda.usc)
nt=100

nn=300

te   <- seq(0, 1, length = nt)

set.seed(1234)

xx <-   rproc2fdata(nn,t=te,sigma="brownian")$data 

xx_PB= xx
for (j in 1:nt){
xx_PB[ ,j]= xx[,j]-te[j]*xx[,nt]
}

set.seed(1234)

xx_OU <-   rproc2fdata(nn,t=te,sigma="OrnsteinUhlenbeck")$data 


set.seed(1234)

xx_VE <-   rproc2fdata(nn,t=te,sigma="vexponential")$data 





matplot(te, t(xx) , type="l", lwd=2, lty=1, xlab="t", ylab="W(t)")


matplot(te, t(xx_PB) , type="l", lwd=2, lty=1, xlab="t", ylab="X(t)")


matplot(te, t(xx_OU) , type="l", lwd=2, lty=1, xlab="t", ylab="Z(t)")


matplot(te, t(xx_VE) , type="l", lwd=2, lty=1, xlab="t", ylab="Y(t)")



algunas=1:10



matplot(te, t(xx[algunas,]) , type="l", lwd=2, lty=1, xlab="t", ylab="W(t)")


matplot(te, t(xx_PB[algunas,]) , type="l", lwd=2, lty=1, xlab="t", ylab="X(t)")


matplot(te, t(xx_OU[algunas,]) , type="l", lwd=2, lty=1, xlab="t", ylab="Z(t)")


matplot(te, t(xx_VE[algunas,]) , type="l", lwd=2, lty=1, xlab="t", ylab="Y(t)")


