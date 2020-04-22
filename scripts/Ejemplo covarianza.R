library(fda.usc)
nt=100

nn=300

grid<-te   <- seq(0, 1, length = nt)

set.seed(1234)

xx <-   rproc2fdata(nn,t=te,sigma="brownian")$data 

xx_PB= xx
for (j in 1:nt){
xx_PB[ ,j]= xx[,j]-te[j]*xx[,nt]
}


set.seed(1234)

xx_VE <-   rproc2fdata(nn,t=te,sigma="vexponential")$data 


cc_xx=colMeans(xx)

xx_center=t(t(xx) - cc_xx)

cov_xx=crossprod(xx_center) / nrow(xx_center)


cov_real=cov_xx

for (i in 1:nt){
 for (j in 1:nt){
	 	cov_real[i,j]=min(te[i],te[j])
	}
}

phi_real=matrix(NA, 4, nt)

for (jota in 1:4){
	phi_real[jota,]= sqrt(2)* sin((jota-0.5)*pi*te)
}

 

dato=phi_real[1,]
sum(dato*dato)/100


phi_est=phi_real

ese=svd(cov_xx) 
cov_xx.svd <- ese$u

for (jota in 1:4){

	phi_est[jota,] <- cov_xx.svd[,jota] 

      signo <- as.numeric(sign(phi_est[jota,]%*%phi_real[jota,]))

	phi_est[jota,] <- signo*phi_est[jota,] *10

}

dato=phi_est[1,]
sum(dato*dato)/100

M= max(phi_est)
m= min(phi_est)


for (jota in 1:4){

	nombre=paste("phi-real-hat-wiener-",jota,".pdf",sep="")
	pdf(nombre, bg='transparent')

	plot(te, phi_real[jota,], type="l", lty=2,lwd=4, xlab="t",ylab=expression(phi[j]),cex.lab=1.2,ylim=c(m,M))
	lines(te,phi_est[jota,], lwd=3, col="maroon")

	dev.off()


}

pdf('media-wiener.pdf', bg='transparent')

matplot(te, t(xx) ,   type="l", lwd=2, lty=1, xlab="t", ylab="W(t)")

lines(te, cc_xx, lwd=4,  col="black")
dev.off()




pdf('Cov-real-wiener.pdf', bg='transparent')
persp(grid, grid, cov_real, xlab="s", ylab="t", 
zlab=" ", theta=30,col="blue", ticktype="detailed" )

dev.off()


pdf('Cov-real-wiener-contour.pdf', bg='transparent')

filled.contour(grid,grid,cov_real)

dev.off()


pdf('Cov-real-wiener-contour2.pdf', bg='transparent')


contour(grid, grid, cov_real, lwd=2)

dev.off()


pdf('Cov-est-wiener.pdf', bg='transparent')

persp(grid, grid, cov_xx, xlab="s", ylab="t", 
zlab=" ", theta=30,col="blue", ticktype="detailed" )
dev.off()


pdf('Cov-est-wiener-contour.pdf', bg='transparent')

filled.contour(grid,grid,cov_xx)
dev.off()

pdf('Cov-est-wiener-contour2.pdf', bg='transparent')

contour(grid, grid, cov_xx, lwd=2)
dev.off()
 

#####################################
#PUENTE
#####################################

cc_xx_PB=colMeans(xx_PB)

xx_PB_center=t(t(xx_PB) - cc_xx_PB)

cov_xx_PB=crossprod(xx_PB_center) / nrow(xx_PB_center)


cov_PB_real=cov_xx_PB

for (i in 1:nt){
 for (j in 1:nt){
	 	cov_PB_real[i,j]=min(te[i],te[j])-te[i]*te[j]
	}
}

phi_PB_real=matrix(NA, 4, nt)

for (jota in 1:4){
	phi_PB_real[jota,]= sqrt(2)* sin( jota*pi*te)
}

 

dato=phi_PB_real[1,]
sum(dato*dato)/100


phi_PB_est=phi_PB_real

ese_PB=svd(cov_xx_PB) 
cov_xx_PB.svd <- ese_PB$u

for (jota in 1:4){

	phi_PB_est[jota,] <- cov_xx_PB.svd[,jota] 

      signo <- as.numeric(sign(phi_PB_est[jota,]%*%phi_PB_real[jota,]))

	phi_PB_est[jota,] <- signo*phi_PB_est[jota,] *10

}

dato=phi_PB_est[1,]
sum(dato*dato)/100

M= max(phi_PB_est)
m= min(phi_PB_est)


for (jota in 1:4){

	nombre=paste("phi-real-hat-puente-",jota,".pdf",sep="")
	pdf(nombre, bg='transparent')

	plot(te, phi_PB_real[jota,], type="l", lty=2,lwd=4, xlab="t",ylab=expression(phi[j]),cex.lab=1.2,ylim=c(m,M))
	lines(te,phi_PB_est[jota,], lwd=3, col="maroon")

	dev.off()


}

pdf('media-puente.pdf', bg='transparent')

matplot(te, t(xx_PB) ,   type="l", lwd=2, lty=1, xlab="t", ylab="W(t)")

lines(te, cc_xx_PB, lwd=4,  col="black")
dev.off()




pdf('Cov-real-puente.pdf', bg='transparent')
persp(grid, grid, cov_PB_real, xlab="s", ylab="t", 
zlab=" ", theta=30,col="blue", ticktype="detailed" )

dev.off()


pdf('Cov-real-puente-contour.pdf', bg='transparent')

filled.contour(grid,grid,cov_PB_real)

dev.off()

 

pdf('Cov-est-puente.pdf', bg='transparent')

persp(grid, grid, cov_xx_PB, xlab="s", ylab="t", 
zlab=" ", theta=30,col="blue", ticktype="detailed" )
dev.off()


pdf('Cov-est-puente-contour.pdf', bg='transparent')

filled.contour(grid,grid,cov_xx_PB)
dev.off()
 
 



#############################
#EXPONENTIAL
####################################
cc_xx_VE=colMeans(xx_VE)

xx_VE_center=t(t(xx_VE) - cc_xx_VE)

cov_xx_VE=crossprod(xx_VE_center) / nrow(xx_VE_center)


cov_real_VE=cov_xx_VE

for (i in 1:nt){
 for (j in 1:nt){
	u=abs(te[i]-te[j])/0.2
	cov_real_VE[i,j]=exp(-u)
	}
}

#######################
# VE
################################
pdf('media-vexp.pdf', bg='transparent')

matplot(te, t(xx_VE) ,   type="l", lwd=2, lty=1, xlab="t", ylab="W(t)")

lines(te, cc_xx_VE, lwd=4,  col="black")
dev.off()



pdf('Cov-real-vexp.pdf', bg='transparent')
persp(grid, grid, cov_real_VE, xlab="s", ylab="t", 
zlab=" ", theta=30,col="blue", ticktype="detailed" )

dev.off()


pdf('Cov-real-vexp-contour.pdf', bg='transparent')

filled.contour(grid,grid,cov_real_VE)

dev.off()


pdf('Cov-real-vexp-contour2.pdf', bg='transparent')


contour(grid, grid, cov_real_VE, lwd=2)

dev.off()


pdf('Cov-est-vexp.pdf', bg='transparent')

persp(grid, grid, cov_xx_VE, xlab="s", ylab="t", 
zlab=" ", theta=30,col="blue", ticktype="detailed" )
dev.off()


pdf('Cov-est-vexp-contour.pdf', bg='transparent')

filled.contour(grid,grid,cov_xx_VE)
dev.off()

pdf('Cov-est-vexp-contour2.pdf', bg='transparent')

contour(grid, grid, cov_xx_VE, lwd=2)
dev.off()


#####################################
# MOVIMIENTO SUAVIZADO
##################################
 

xx_smooth=xx

for (j in 1:nn){

	y=xx[j,]
	xx_smooth[ j, ]= ksmooth(te, y, kernel = "normal", bandwidth = 0.1, x.points=te)$y

}



cc_xx=colMeans(xx_smooth)

xx_center=t(t(xx_smooth) - cc_xx)

cov_xx=crossprod(xx_center) / nrow(xx_center)

 

ese=svd(cov_xx) 
cov_xx.svd <- ese$u

for (jota in 1:4){

	phi_est[jota,] <- cov_xx.svd[,jota] 

      signo <- as.numeric(sign(phi_est[jota,]%*%phi_real[jota,]))

	phi_est[jota,] <- signo*phi_est[jota,] *10

}

dato=phi_est[1,]
sum(dato*dato)/100

M= max(phi_est)
m= min(phi_est)


for (jota in 1:4){

	nombre=paste("phi-real-hat-wiener-smooth-",jota,".pdf",sep="")
	pdf(nombre, bg='transparent')

	plot(te, phi_real[jota,], type="l", lty=2,lwd=4, xlab="t",ylab=expression(phi[j]),cex.lab=1.2,ylim=c(m,M))
	lines(te,phi_est[jota,], lwd=3, col="maroon")

	dev.off()


}

pdf('media-wiener-smooth.pdf', bg='transparent')

matplot(te, t(xx_smooth) ,   type="l", lwd=2, lty=1, xlab="t", ylab="W(t)")

lines(te, cc_xx, lwd=4,  col="black")
dev.off()


 

pdf('Cov-est-wiener-smooth.pdf', bg='transparent')

persp(grid, grid, cov_xx, xlab="s", ylab="t", 
zlab=" ", theta=30,col="blue", ticktype="detailed" )
dev.off()


pdf('Cov-est-wiener-smooth-contour.pdf', bg='transparent')

filled.contour(grid,grid,cov_xx)
dev.off()

pdf('Cov-est-wiener-smooth-contour2.pdf', bg='transparent')

contour(grid, grid, cov_xx, lwd=2)
dev.off()
 

#####################################
# MOVIMIENTO SUAVIZADO SPLINES
##################################
 library(fda)


B25.basis=create.bspline.basis(rangeval=c(0,1), nbasis=25) #nbasis = nbreaks + norder - 2 y norder=4

W.fd=smooth.basis(argvals=te,y=t(xx), fdParobj=B25.basis)


plot(W.fd, ylab="", xlab="",col="gray",lty=1,lwd=3)



cc_xx=mean(W.fd$fd)
 

W.sd=std.fd(W.fd$fd)

W.cov=var.fd(W.fd$fd) # $fd extracts function values

cov_xx=eval.bifd(te, te, W.cov) 

ese=svd(cov_xx) 
cov_xx.svd <- ese$u

for (jota in 1:4){

	phi_est[jota,] <- cov_xx.svd[,jota] 

      signo <- as.numeric(sign(phi_est[jota,]%*%phi_real[jota,]))

	phi_est[jota,] <- signo*phi_est[jota,] *10

}

dato=phi_est[1,]
sum(dato*dato)/100

M= max(phi_est)
m= min(phi_est)


for (jota in 1:4){

	nombre=paste("phi-real-hat-wiener-BS-",jota,".pdf",sep="")
	pdf(nombre, bg='transparent')

	plot(te, phi_real[jota,], type="l", lty=2,lwd=4, xlab="t",ylab=expression(phi[j]),cex.lab=1.2,ylim=c(m,M))
	lines(te,phi_est[jota,], lwd=3, col="maroon")

	dev.off()


}

pdf('media-wiener-BS.pdf', bg='transparent')


plot(W.fd, ylab="", xlab="",lty=1,lwd=2)


lines(cc_xx, lwd=4,  col="black")
dev.off()


 

pdf('Cov-est-wiener-BS.pdf', bg='transparent')

persp(grid, grid, cov_xx, xlab="s", ylab="t", 
zlab=" ", theta=30,col="blue", ticktype="detailed" )
dev.off()


pdf('Cov-est-wiener-BS-contour.pdf', bg='transparent')

filled.contour(grid,grid,cov_xx)
dev.off()

pdf('Cov-est-wiener-BS-contour2.pdf', bg='transparent')

contour(grid, grid, cov_xx, lwd=2)
dev.off()
 

