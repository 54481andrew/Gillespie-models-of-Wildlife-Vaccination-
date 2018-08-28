
require(deSolve)
source('Tools/Functions.r', local = TRUE)
tvvals = seq(0,365,length.out = 10)
parms = expand.grid(b0 = 400, tb = 90, tv = tvvals, d = 0.004, T = 365, Bp = 0.0000005, gamp = 0.005, gamv = 0.07, NVacc = 5000, IpAvg = 0, MaxTimeIp = 0)
times = seq(0,10000)
y0 = c(1000, 0, 100, 0, 00)
names(y0) = c('S','Iv','Ip','V','P')

maxtimes <- 365*10
timeseq  <- seq(0,maxtimes, by = 0.01)

###Get statistic to graph (average)
NPars <- nrow(parms)
ExtTimesList = list()
for(NPar in 1:NPars){
    print(paste("NPar:", NPar))
    datname = paste("Data/NPar_", NPar-1, sep = '')
    dat = read.table(file = datname, header = T, sep = " ")
    
    tvval = parms$tv[NPar]    
    VaccPer <- 365
    NVacc   <- parms$NVacc[NPar]
    vacctimes <- seq(tvval, maxtimes,by = VaccPer)
    
    out <- data.frame(lsoda(y = y0, times = times, 
                   func = rhs.fun, parms = parms[NPar,],
                   events=list(func = vaccinate, time=vacctimes),
                   maxsteps = 100000))

    ###Calculate statistic (average)
    fode.Ip = approxfun(out$time, out$Ip)
        
    ###Recast stochastic data    
    startind = c(which(dat$time==0),nrow(dat))
    datlist = list()
    ExtTimes = c()
    for(i in 1:(length(startind)-1)){
    	  dati = dat[startind[i]:(startind[i+1] - 1),]
    	  datlist[[i]] = dati
	  ExtTimei = dati$time[max(which(dati$Ip > 0))]
	  ExtTimes =  c(ExtTimes,ExtTimei)
#	  parms$MaxTimeIp[i] = mean(ExtTimes[i])
    }
 
ExtTimesList[[NPar]] = ExtTimes

#pdf(file = paste(datname,'Tv_',round(tvval), ".pdf",sep = ""))
plot(S~time, out, xlim = c(0,365*5), ylim = c(0,50000), type = 'l', lwd = 4)
lines(Iv~time, out, col = 'green', lty = 1, lwd = 4)
lines(Ip~time, out, col = 'red', lty = 1, lwd = 4)
lines(V~time, out, col = 'blue', lty = 1, lwd = 4)
lines(P~time, out, col = 'pink', lty = 1, lwd = 4)

these <- seq(1,nrow(dat), by = 5 )
points(S~time, dat[these,], col = 'black', cex = 1, pch = 3)
points(Iv~time, dat[these,], col = 'green', cex = 1, pch = 3)
points(Ip~time, dat[these,], col = 'red', cex = 1, pch = 3)
points(V~time, dat[these,], col = 'blue', cex = 1, pch = 3)
points(P~time, dat[these,], col = 'pink', cex = 1, pch = 3)


abline(v = vacctimes, lty = 3)

#dev.off()

}



