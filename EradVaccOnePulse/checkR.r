SimName = "DeerMice_Base_LP_Ex";
parmat = read.table(file = paste("Data/", SimName, "/ParMat", sep=''), header = F)
names(parmat) = c('Par','b0','d','Bp','Nv','tv','gamv','gamp','tb','T','IpInit', 'TVaccStart')

if(grepl('Freq',SimName)){
	parmat$R0p = with(parmat, round( Bp/(d+gamp) ,2))
	ode.fun = rhs.freq.fun
}else{
	parmat$R0p = with(parmat, round( Bp*(b0*tb)/(T*d*(d+gamp)) ,2))
	ode.fun = rhs.fun
}

parmat$R0approx = with(parmat, R0p*(1-Nv/(b0*tb + Nv*exp(-d*(T-tv)))))
parmat$rho = with(parmat, round(Nv/(b0*tb/(d*T)),2))


require(deSolve)
source('Tools/Functions.r', local = TRUE)

times = seq(0,10000)

maxtimes <- 365*11
timeseq  <- seq(0,maxtimes, by = 0.01)
NPars = 8#length(parmat$Par)


for(i in 1:NPars){

y0 = c(1000, 0, 100, 0, 0)
names(y0) = c('S','Iv','Ip','V','P')

    VaccPer        <- parmat$T[i]
    NVacc          <- parmat$Nv[i]
    parmat$Nv[i] <- 0 #temporarily set to 0

    tv <- parmat$tv[i]
    vacctimes <- seq(tv, maxtimes,by = VaccPer)
   
    print(paste("NPar:", i-1))
    filename = paste("Data/", SimName, "/Par_", i-1, sep = '')
    dat = read.table(file = filename, header = T, sep = " ")

    tvaccstart <- parmat$TVaccStart[i] + tv 
    #vacctimespre <- vacctimes[vacctimes < tvaccstart] #Vaccination starts 5 years in
    out1 <- data.frame(lsoda(y = y0, times = seq(0,tvaccstart,by=1), 
                            func = ode.fun, parms = parmat[i,],
                            maxsteps = 100000))
    y0 = out1[nrow(out1), -1]
    y0 = as.numeric(y0)
    parmat$Nv[i] <- NVacc
    y0 = vaccinate(0,y0,parmat[i,])
    names(y0) = c('S','Iv','Ip','V','P')
    vacctimespost <- vacctimes[vacctimes >= tvaccstart]
    out2 <- data.frame(lsoda(y = y0, times = seq(tvaccstart, maxtimes, by = 1), 
                            func = ode.fun, parms = parmat[i,],
                            maxsteps = 100000))
    out = rbind(out1, out2[-1,])			

    #Big Picture
    pdf(file = paste(filename, "_big.pdf",sep = ""), height = 4, width = 6)
    par(mai = c(1,1,0.1,0.1))
    plot(NA, xlim = c(0,maxtimes), ylim = c(0,1000), type = 'l', lwd = 1, col = 'gray',xlab = 'Time', ylab = 'Density')
    lines(Iv~time, out, col = 'green', lty = 1, lwd = 1)

    lines(V~time, out, col = 'blue', lty = 1, lwd = 1)
    lines(P~time, out, col = 'pink', lty = 1, lwd = 1)
    
    these <- seq(1,nrow(dat), by = 100)
    points(S~time, dat[these,], col = 'black', cex = 0.25, pch = 1)
    points(Iv~time, dat[these,], col = 'green', cex = 0.25, pch = 1)
    points(Ip~time, dat[these,], col = 'red', cex = 0.25, pch = 1)
#    points(V~time, dat[these,], col = 'blue', cex = 0.25, pch = 1)
#    points(P~time, dat[these,], col = 'pink', cex = 0.25, pch = 1)
    lines(Ip~time, out, col = 'black', lty = 1, lwd = 1)
    legend = c('S','Iv','Ip','V','P')
    legend(x = "topright", legend = legend, col = c('gray', 'green', 'red', 'blue', 'pink'), lwd = 2)
    
    lines(S~time, out, type = 'l', lwd = 1, col = 'gray')

    abline(v = vacctimes, lty = 3)

    dev.off()

#########################

    #Small Picture
    pdf(file = paste(filename, "_small.pdf",sep = ""), height = 4, width = 6)
    par(mai = c(1,1,0.1,0.1))
    plot(NA, xlim = c(3000-10,3000+365), ylim = c(0,350), type = 'l', lwd = 1, col = 'gray', xlab = 'Time', ylab = 'Density')
    lines(Iv~time, out, col = 'green', lty = 1, lwd = 1)

    lines(V~time, out, col = 'blue', lty = 1, lwd = 1)
    lines(P~time, out, col = 'pink', lty = 1, lwd = 1)
    
    these <- seq(1,nrow(dat), by = 100)
    points(S~time, dat[these,], col = 'black', cex = 0.25, pch = 1)
    points(Iv~time, dat[these,], col = 'green', cex = 0.25, pch = 1)
    points(Ip~time, dat[these,], col = 'red', cex = 0.25, pch = 1)

    #points(V~time, dat[these,], col = 'blue', cex = 0.25, pch = 1)
    #points(P~time, dat[these,], col = 'pink', cex = 0.25, pch = 1)

    lines(S~time, out, type = 'l', lwd = 1, col = 'gray')
    lines(Ip~time, out, col = 'black', lty = 1, lwd = 1)
    legend = c('S','Iv','Ip','V','P')
    legend(x = "topright", legend = legend, col = c('gray', 'green', 'red', 'blue', 'pink'), lwd = 2)
    
    abline(v = vacctimes, lty = 3)

    dev.off()

#########################

    #Tiny Picture
    pdf(file = paste(filename, "_tiny.pdf",sep = ""), height = 4, width = 6)
    par(mai = c(1,1,0.1,0.1))
    plot(NA, xlim = c(0,maxtimes), ylim = c(0,100), type = 'l', lwd = 1, col = 'gray', xlab = 'Time', ylab = 'Density')
    lines(Iv~time, out, col = 'green', lty = 1, lwd = 1)

#    lines(V~time, out, col = 'blue', lty = 1, lwd = 1)
#    lines(P~time, out, col = 'pink', lty = 1, lwd = 1)
    
    these <- seq(1,nrow(dat), by = 100)
    points(S~time, dat[these,], col = 'black', cex = 0.25, pch = 1)
    points(Iv~time, dat[these,], col = 'green', cex = 0.25, pch = 1)
    points(Ip~time, dat[these,], col = 'red', cex = 0.25, pch = 1)
#    points(V~time, dat[these,], col = 'blue', cex = 0.25, pch = 1)
#    points(P~time, dat[these,], col = 'pink', cex = 0.25, pch = 1)

    lines(S~time, out, type = 'l', lwd = 1, col = 'gray')

    lines(Ip~time, out, col = 'black', lty = 1, lwd = 1)
    legend = c('S','Iv','Ip','V','P')
    legend(x = "topright", legend = legend, col = c('gray', 'green', 'red', 'blue', 'pink'), lwd = 2)
    
    abline(v = vacctimes, lty = 3)

    dev.off()

}
   






