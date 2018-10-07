SimName = "LagImport"

parmat = expand.grid(NPar = 0, b0 = NA, d = 0.002, rho = seq(0,1,length.out = 10),
                     tv = seq(0.1,365,length.out = 20), gamv =0.07, tb = seq(15,365,length.out = 10),
                     T = 365, NAvg = 10000, VAvg = 0)

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

mcl.fun = function(i){

    y0 = c(1000, 0, 0)
    names(y0) = c('S','Iv','V')

    T  <- parmat$T[i]
    NVacc <- parmat$Nv[i]
    tv <- parmat$tv[i]
    vacctimes <- seq(tv, maxtimes,by = T)
   
    print(paste("NPar:", i, ' / ', NPars))

    out  <- data.frame(lsoda(y = y0, times = seq(0, maxtimes, by=1), 
                            func = rhs.freq.jv.fun, parms = parmat[i,],
                            events=list(func = vaccinate.jv, time=vacctimes),
                            maxsteps = 100000))

    ###Store average serorpevalence
    parmat$VAvg[i] <- mean(with(out, V/(S + Iv + V)))
    write.table(parmat[i,], file = filename, append = T, col.names = F, row.names = F)
    
    return(i)
}

metaout <- mclapply(1:NPars, mcl.fun)


refparmat <- read.table(file = filename, header = T)
refparmat <- refparmat[order(refparmat$NPar),]
write.table(refparmat, file = filename, append = F, row.names = F)

parmat <- read.table(file = filename, header = T)
parmat$tvstar <- NA

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
        wimax <- which.max(parmat$VAvg[wi])
        TvStar <- parmat$tv[wi][wimax]
        VAvgStar <- parmat$VAvg[wi][wimax]
        VNull <- mean(parmat$VAvg[wi])
        Lag.Mat[xi,yi] <- TvStar - unique(parmat$tb[wi])
        Import.Mat[xi,yi] <- VAvgStar/VNull
    }
}
