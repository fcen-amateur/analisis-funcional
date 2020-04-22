library(fda)
data(daily)
names(daily)

dim(daily$tempav)

matplot(daily$tempav, type="l", lty=1, col="gray", lwd=3, xlab="Dia", ylab="Temperatura")


data(CanadianWeather)
names(CanadianWeather)

dim(CanadianWeather$dailyAv)

matplot(day.5, CanadianWeather$dailyAv[,, "Temperature.C"], type="l", lty=1, col="gray", lwd=3, xlab="Dia", ylab="Temperatura")
 

matplot(1:12, CanadianWeather$monthlyTemp, type="l", lty=1, col="gray", lwd=3, xlab="Mes", ylab="Temperatura")

stations <- c("Pr. Rupert", "Montreal", "Edmonton", "Resolute")

matplot(day.5, CanadianWeather$dailyAv[, stations, "Temperature.C"], type="l", axes=FALSE, xlab="", ylab="Temperatura",lty=1, lwd=3)
axis(2, las=1)
# Label the horizontal axis with the month names
axis(1, monthBegin.5, labels=FALSE)
axis(1, monthEnd.5, labels=FALSE)

axis(1, monthMid, monthLetters, tick=FALSE)
# Add the monthly averages
colores=c("gray", "maroon", "darkolivegreen4", "deepskyblue4")
matplot(monthMid, CanadianWeather$monthlyTemp[, stations], col=colores, type="l", lty=2, lwd=2, pch=16,cex=1.2, add=T)
# Add the names of the weather stations
mtext(stations, side=4,at=CanadianWeather$dailyAv[365, stations, "Temperature.C"],
las=1)




matplot(1:12, CanadianWeather$monthlyTemp, type="l", lty=1, col="gray", lwd=3, xlab="Mes", ylab="Temperatura")

stations <- c("Pr. Rupert", "Montreal", "Edmonton", "Resolute")

colores=c("black", "red", "green4", "blue")

matplot(1:12, CanadianWeather$monthlyTemp[,stations], type="l", lty=1, col=colores, lwd=3, add=T)
mtext(stations, side=4,at=CanadianWeather$dailyAv[365, stations, "Temperature.C"],
las=1, col=colores)


aa=fbplot(CanadianWeather$monthlyTemp, x = 1:12, ylim=c(-34,21), xlab="Mes", ylab="Temperatura")


#matplot(1:12, CanadianWeather$monthlyTemp, type="l", lty=1, col="gray", lwd=1,add=T)


aa$out

nombres=CanadianWeather$place[aa$out]

mtext(nombres[1:3], side=4,at=CanadianWeather$monthlyTemp[12,nombres[1:3]],las=1)


mtext(nombres[4:6], side=2,at=CanadianWeather$monthlyTemp[12,nombres[4:6]],las=1)

