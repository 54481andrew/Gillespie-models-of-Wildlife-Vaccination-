datname = 'DownPhase'
dat = read.table(file = datname, header = T, sep = " ")
names(dat)

require(deSolve)
source('Tools/Functions.r', local = TRUE)

parms = data.frame(b0 = 400, tb = 90, tv = 180, d = 0.004, Tb = 365, Bp = 0.0000005, gamp = 0.005, gamv = 0.07, NVacc = 5000)
times = seq(0,10000)
y0 = c(1000, 0, 100, 0, 00)
names(y0) = c('S','Iv','Ip','V','P')

maxtimes <- 365*10
timeseq  <- seq(0,maxtimes, by = 0.01)

VaccPer        <- 365
NVacc          <- parms$NVacc

tv <- parms$tv
vacctimes <- seq(tv, maxtimes,by = VaccPer)

out <- data.frame(lsoda(y = y0, times = times, 
                   func = rhs.fun, parms = parms,
                   events=list(func = vaccinate, time=vacctimes),
                   maxsteps = 100000))



pdf(file = paste(datname, ".pdf",sep = ""))
plot(S~time, out, xlim = c(0,365*5), ylim = c(0,50000), type = 'l', lwd = 4)
lines(Iv~time, out, col = 'green', lty = 1, lwd = 4)
lines(Ip~time, out, col = 'red', lty = 1, lwd = 4)
lines(V~time, out, col = 'blue', lty = 1, lwd = 4)
lines(P~time, out, col = 'pink', lty = 1, lwd = 4)

these <- seq(1,nrow(dat), by = 25)
points(S~time, dat[these,], col = 'black', cex = 1, pch = 3)
points(Iv~time, dat[these,], col = 'green', cex = 1, pch = 3)
points(Ip~time, dat[these,], col = 'red', cex = 1, pch = 3)
points(V~time, dat[these,], col = 'blue', cex = 1, pch = 3)
points(P~time, dat[these,], col = 'pink', cex = 1, pch = 3)


abline(v = vacctimes, lty = 3)

dev.off()





