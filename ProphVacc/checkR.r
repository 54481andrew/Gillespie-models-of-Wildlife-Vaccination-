SimName = "Test"
parmat = read.table(file = paste("Data/ParMat_", SimName, sep=''), header = F)
names(parmat) = c('NPar','b0','d','Bp','Nv','tv','gamv','gamp','tb','T','IpInit','TPathInv')

require(deSolve)
source('Tools/Functions.r', local = TRUE)

times = seq(0,10000)
y0 = c(10000, 0, 0, 0, 00)
names(y0) = c('S','Iv','Ip','V','P')

maxtimes <- 365*10
timeseq  <- seq(0,maxtimes, by = 0.01)

for(i in 1:length(parmat$NPar)){

    VaccPer        <- parmat$T[i]
    NVacc          <- parmat$Nv[i]

    tv <- parmat$tv[i]
    vacctimes <- seq(tv, maxtimes,by = VaccPer)
   
    print(paste("NPar:", i-1))
    filename = paste("Data/", SimName, "/Par_", i-1, sep = '')
    dat = read.table(file = filename, header = T, sep = " ")

    tpathinv <- parmat$TPathInv[i]
    vacctimespre <- vacctimes[vacctimes < tpathinv]
    out1 <- data.frame(lsoda(y = y0, times = seq(0,tpathinv,by=1), 
                            func = rhs.fun, parms = parmat[i,],
                            events=list(func = vaccinate, time=vacctimespre),
                            maxsteps = 100000))
    y0 = out1[nrow(out1), -1]
    y0[3] = parmat$IpInit[i]
    y0 = as.numeric(y0)
    names(y0) = c('S','Iv','Ip','V','P')
    vacctimespost <- vacctimes[vacctimes >= tpathinv]
    out2 <- data.frame(lsoda(y = y0, times = seq(tpathinv, 10*365, by = 1), 
                            func = rhs.fun, parms = parmat[i,],
                            events=list(func = vaccinate, time=vacctimespost),
                            maxsteps = 100000))
    out = rbind(out1, out2[-1,])			

    #Big Picture
    pdf(file = paste(filename, "_big.pdf",sep = ""), height = 4, width = 6)
    par(mai = c(1,1,0.1,0.1))
    plot(S~time, out, xlim = c(0,365*10), ylim = c(0,35000), type = 'l', lwd = 1)
    lines(Iv~time, out, col = 'green', lty = 1, lwd = 1)
    lines(Ip~time, out, col = 'red', lty = 1, lwd = 1)
    lines(V~time, out, col = 'blue', lty = 1, lwd = 1)
    lines(P~time, out, col = 'pink', lty = 1, lwd = 1)

    ce = 0.15;
    these <- seq(1,nrow(dat), by = 200)
    points(S~time, dat[these,], col = 'black', cex = ce, pch = 1)
    points(Iv~time, dat[these,], col = 'green', cex = ce, pch = 1)
    points(Ip~time, dat[these,], col = 'red', cex = ce, pch = 1)
    points(V~time, dat[these,], col = 'blue', cex = ce, pch = 1)
    points(P~time, dat[these,], col = 'pink', cex = ce, pch = 1)

    legend = c('S','Iv','Ip','V','P')
    legend(x = "topright", legend = legend, col = c('black', 'green', 'red', 'blue', 'pink'), lwd = 2)
    
    abline(v = vacctimes, lty = 3)

    dev.off()

#########################

    #Small Picture
    pdf(file = paste(filename, "_small.pdf",sep = ""), height = 4, width = 6)
    par(mai = c(1,1,0.1,0.1))
    plot(S~time, out, xlim = c(0,10*365), ylim = c(0,1000), type = 'l', lwd = 4)
    lines(Iv~time, out, col = 'green', lty = 1, lwd = 4)
    lines(Ip~time, out, col = 'red', lty = 1, lwd = 4)
    lines(V~time, out, col = 'blue', lty = 1, lwd = 4)
    lines(P~time, out, col = 'pink', lty = 1, lwd = 4)
    
    these <- seq(1,nrow(dat), by = 100)
    points(S~time, dat[these,], col = 'black', cex = 1, pch = 3)
    points(Iv~time, dat[these,], col = 'green', cex = 1, pch = 3)
    points(Ip~time, dat[these,], col = 'red', cex = 1, pch = 3)
    points(V~time, dat[these,], col = 'blue', cex = 1, pch = 3)
    points(P~time, dat[these,], col = 'pink', cex = 1, pch = 3)

    legend = c('S','Iv','Ip','V','P')
    legend(x = "topright", legend = legend, col = c('black', 'green', 'red', 'blue', 'pink'), lwd = 2)
    
    abline(v = vacctimes, lty = 3)

    dev.off()

#########################

    #Tiny Picture
    pdf(file = paste(filename, "_tiny.pdf",sep = ""), height = 4, width = 6)
    par(mai = c(1,1,0.1,0.1))
    plot(S~time, out, xlim = c(tpathinv - 30,3500), ylim = c(0,20), type = 'l', lwd = 4)
    lines(Iv~time, out, col = 'green', lty = 1, lwd = 4)

    lines(V~time, out, col = 'blue', lty = 1, lwd = 4)
    #lines(P~time, out, col = 'pink', lty = 1, lwd = 4)
    
    these <- seq(1,nrow(dat), by = 1)
    #points(S~time, dat[these,], col = 'black', cex = 1, pch = 1)
    #points(Iv~time, dat[these,], col = 'green', cex = 1, pch = 1)
    points(Ip~time, dat[these,], col = 'red', cex = 1, pch = 1)
    #points(V~time, dat[these,], col = 'blue', cex = 1, pch = 1)
    #points(P~time, dat[these,], col = 'pink', cex = 1, pch = 3)

    lines(Ip~time, out, col = 'red', lty = 1, lwd = 4)
    
    legend = c('S','Iv','Ip','V','P')
    legend(x = "topright", legend = legend, col = c('black', 'green', 'red', 'blue', 'pink'), lwd = 2)
    
    abline(v = vacctimes, lty = 3)

    dev.off()

}
   






