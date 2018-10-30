require(parallel)
SimName = "LagImport"

parmat = expand.grid(NPar = 0, b0 = NA, Rp = c(1.5,2,2.5),d = 0.00548, rho = c(0.25, 0.5),
                     tv = seq(0.1,365,length.out = 20), gamp = c(0.01,0.03, 0.05), gamv =0.07, tb = 90,
                     T = 365, NAvg = 4500, VAvg = 0, Ipmin = 0, Ipmax = 0, ipmin = 0, ipmax = 0)

NPars = nrow(parmat)
parmat$NPar = 1:NPars
parmat$b0 = with(parmat, (d*T)*NAvg/tb) ###NAvg = b0 tb /(d*T)
parmat$Nv = with(parmat,rho*NAvg)
parmat$Bp = with(parmat, Rp/NAvg*(d+gamp))

require(deSolve)
source('Tools/Functions.r', local = TRUE)

maxtimes <- 365*40
tBurn    <- 365*38
times    <- seq(0,maxtimes, by = 0.01)

filename = paste("ODEData/", SimName, sep = '')
write.table(parmat[1,][-1,], file = filename, append = F)

mcl.fun = function(i){



    y0 = c(1000, 0, 100, 0, 0)
    names(y0) = c('S','Iv','Ip','V','P')

    T  <- parmat$T[i]
    NVacc <- parmat$Nv[i]
    tv <- parmat$tv[i]
    vacctimes <- seq(tv, maxtimes,by = T)
   
    print(paste("NPar:", i, ' / ', NPars))

    out  <- data.frame(lsoda(y = y0, times = seq(0, maxtimes, by=1), 
                            func = rhs.fun, parms = parmat[i,],
                            events=list(func = vaccinate, time=vacctimes),
                            maxsteps = 100000))
    out$NTot = out$S + out$Iv + out$Ip + out$V + out$P			 


    ###Store average serorpevalence
    these <- out$time > tBurn
    parmat$VAvg[i] <- mean(with(out, V/(S + Iv + V)))
    parmat$Ipmin[i] <- min(out$Ip[these])
    parmat$Ipmax[i] <- max(out$Ip[these]) 
    parmat$ipmin[i] <- min(out$Ip[these] / out$NTot)
    parmat$ipmax[i] <- max(out$Ip[these] / out$NTot)    
    


    write.table(parmat[i,], file = filename, append = T, col.names = F, row.names = F)
    
    return(i)
}

metaout <- mclapply(1:nrow(parmat), mcl.fun, mc.cores = 14)


refparmat <- read.table(file = filename, header = T)
refparmat <- refparmat[order(refparmat$NPar),]
write.table(refparmat, file = filename, append = F, row.names = F)

parmat <- read.table(file = filename, header = T)
parmat$tvstar <- NA


XName <- 'tv'
YName <- 'Rp'
XVals  <- unique(parmat[,XName])
YVals  <- unique(parmat[,YName])
nXVals <- length(XVals)
nYVals <- length(YVals)

OValName1 <- 'rho'
OVals1    <- c(0.25,0.5) 
OVal1 <- 0.5
OValName2 <- 'gamp'
OVals2    <- c(0.03) 
OVal2 <- 0.03


MinMat = matrix(NA,nrow = nXVals, ncol = nYVals)
MaxMat = matrix(NA,nrow = nXVals, ncol = nYVals)

for(xi in 1:nXVals){
    for(yi in 1:nYVals){
        wi <- XVals[xi]==parmat[,XName] & YVals[yi]==parmat[,YName] & 
	      OVal1==parmat[,OValName1] & OVal2==parmat[,OValName2]
        MinMat[xi,yi] <- parmat$Ipmin[wi]
        MaxMat[xi,yi] <- parmat$Ipmax[wi]
    }
}

matplot(MinMat, ylim = c(0, max(MinMat)))
matpoints(MaxMat)

#write.table(refparmat, file = filename, append = F, row.names = F)