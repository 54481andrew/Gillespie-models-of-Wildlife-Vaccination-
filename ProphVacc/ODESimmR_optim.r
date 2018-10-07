SimName = "LagImport"

parmat = expand.grid(NPar = 0, b0 = NA, d = 0.002, rho = seq(0.01,1,length.out = 5),
                     gamv =0.07, tb = seq(15,365,length.out = 5),
                     T = 365, NAvg = 10000, tvstar = 0, vavgstar = 0, vnull = 0)

NPars = nrow(parmat)
parmat$NPar = 1:NPars
parmat$b0 = with(parmat, (d*T)*NAvg/tb) ###NAvg = b0 tb /(d*T)
parmat$Nv = with(parmat,rho*NAvg)

require(deSolve)
source('Tools/Functions.r', local = TRUE)

maxtimes <- 365*20
times    <- seq(0,maxtimes, by = 0.01)

filename = paste("ODEData/", SimName, sep = '')
write.table(parmat[1,][-1,], file = filename, append = F)

tvvals <- matrix(seq(0.1,364.9,length.out = 20), ncol = 1)

mcl.fun = function(i){

    y0 = c(1000, 0, 0)
    names(y0) = c('S','Iv','V')
    TVal  <- parmat$T[i]
    
    cost.fxn <- function(tv,parms){
        NVacc <- parms$Nv
        vacctimes <- seq(tv, maxtimes,by = TVal)        

        out  <- data.frame(lsoda(y = y0, times = seq(0, maxtimes, by=1), 
                                 func = rhs.freq.jv.fun, parms = parms,
                                 events=list(func = vaccinate.jv, time=vacctimes),
                                 maxsteps = 100000))
        
###Return average serorpevalence
        VAvg <- mean(with(out, V/(S + Iv + V)))
        return(VAvg)
    }###End cost.fxn

    vnullvals <- apply(X = tvvals, MARGIN = 1, FUN = cost.fxn, parms = parmat[i,])
    tvstar <- tvvals[which.max(vnullvals)]
    uroot <- optimize(f=cost.fxn, interval = c(0.67*tvstar, min(1.33*tvstar,365)), parmat[i,], maximum=TVal)
    tvstar <- uroot$maximum ###Location
    vavgstar <- uroot$objective ###Maximum value

    parmat$tvstar[i] <- tvstar
    parmat$vavgstar[i] <- vavgstar
    parmat$vnull[i]  <- mean(vnullvals)
    write.table(parmat[i,], file = filename, append = T, col.names = F, row.names = F)
    print(paste("NPar:", i, ' / ', NPars))
    return(i)
}

metaout <- mclapply(X = 1:NPars, FUN = mcl.fun)

refparmat <- read.table(file = filename, header = T)
refparmat <- refparmat[order(refparmat$NPar),]
write.table(refparmat, file = filename, append = F, row.names = F)
parmat <- read.table(file = filename, header = T)

XName <- 'tb'
YName <- 'rho'
XVals  <- unique(parmat[,XName])
YVals  <- unique(parmat[,YName])
nXVals <- length(XVals)
nYVals <- length(YVals)

Import.Mat = matrix(NA,nrow = nXVals, ncol = nYVals)
Lag.Mat = matrix(NA,nrow = nXVals, ncol = nYVals)
for(xi in 1:nXVals){
    for(yi in 1:nYVals){
        wi <- XVals[xi]==parmat[,XName] & YVals[yi]==parmat[,YName]
        TvStar <- parmat$tvstar[wi]
        VAvgStar <- parmat$vavgstar[wi]
        VNull <-  parmat$vnull[wi]
        Lag.Mat[xi,yi] <- TvStar - unique(parmat$tb[wi])
        Import.Mat[xi,yi] <- VAvgStar/VNull
    }
}
