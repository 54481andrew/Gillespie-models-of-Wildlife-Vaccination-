datname = 'test'
parmat = read.table(file = "Data/ParMat", header = T)

require(deSolve)
source('Tools/Functions.r', local = TRUE)

times = seq(0,10000)
y0 = c(1000, 0, 100, 0, 00)
names(y0) = c('S','Iv','Ip','V','P')

maxtimes <- 365*10
timeseq  <- seq(0,maxtimes, by = 0.01)


for(i in 1:length(parmat$NPar)){

    VaccPer        <- parmat$T[i]
    NVacc          <- parmat$Nv[i]

    
    tv <- parmat$tv[i]
    vacctimes <- seq(tv, maxtimes,by = VaccPer)
    
    print(paste("NPar:", i-1))
    datname = paste("Data/NPar_", i-1, sep = '')
    dat = read.table(file = datname, header = T, sep = " ")
    
    out <- data.frame(lsoda(y = y0, times = times, 
                            func = rhs.fun, parms = parmat[i,],
                            events=list(func = vaccinate, time=vacctimes),
                            maxsteps = 100000))
    
    pdf(file = paste(datname, ".pdf",sep = ""))
    plot(S~time, out, xlim = c(0,365*5), ylim = c(0,30000), type = 'l', lwd = 4)
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

    legend = c('S','Iv','Ip','V','P')
    legend(x = "topright", legend = legend, col = c('black', 'green', 'red', 'blue', 'pink'), lwd = 2)
    
    abline(v = vacctimes, lty = 3)

    dev.off()
}
   






