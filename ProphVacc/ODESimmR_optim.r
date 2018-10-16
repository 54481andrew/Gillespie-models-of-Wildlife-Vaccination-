require(parallel)

SimName = "LagImport.TMax40.1yr"

parmat = expand.grid(NPar = 0, b0 = NA, d = c(0.00274), 
       	             rho = seq(0.01,1,length.out = 50),
                     gamv =0.07, tb = seq(15,365,length.out = 50),
                     T = 365, NAvg = 10000, tvstar = 0, 
		     vavgstar = 0, vnull = 0)

NPars = nrow(parmat)
parmat$NPar = 1:NPars
parmat$b0 = with(parmat, (d*T)*NAvg/tb) ###NAvg = b0 tb /(d*T)
parmat$Nv = with(parmat,rho*NAvg)

require(deSolve)
source('Tools/Functions.r', local = TRUE)

nyears   <- 40
maxtimes <- 365*nyears
times    <- seq(0,maxtimes, by = 1)
witimes <- times >= (nyears-1)*365
yeartimes <- c(times[witimes]%%365)
yeartimes[length(yeartimes)] <- 365
tvvals <- matrix(seq(0.1,364.9,length.out = 25), ncol = 1)

filename = paste("ODEData/", SimName, sep = '')
write.table(parmat[1,][-1,], file = filename, append = F)



mcl.fun = function(i){

    y0 = c(1000, 0, 0)
    names(y0) = c('S','Iv','V')
    TVal  <- parmat$T[i]
    
    cost.fxn <- function(tv,parms){
        NVacc <- parms$Nv
        vacctimes <- seq(tv, maxtimes,by = TVal)        

        out  <- data.frame(lsoda(y = y0, times = times, 
                                 func = rhs.freq.jv.fun, parms = parms,
                                 events=list(func = vaccinate.jv, time=vacctimes),
                                 maxsteps = 100000))
        
###Return average serorpevalence
	vals    <- with(out[witimes,], V/(S + Iv + V))
	vals[1] <- vals[length(vals)]
	vapprox = splinefun(x = yeartimes, vals, method = "periodic")
        VAvg <- 1/365*(integrate(vapprox, lower = 0, upper = 365, stop.on.error = FALSE)$value)
        return(VAvg)
    }###End cost.fxn

    vnullvals <- apply(X = tvvals, MARGIN = 1, FUN = cost.fxn, parms = parmat[i,])
    tvstar <- tvvals[which.max(vnullvals)]

    wi.interval <- c(max(1, which.max(vnullvals) - 1), min(which.max(vnullvals) + 1, length(tvvals)))
    interval <- tvvals[wi.interval]

    uroot <- optimize(f=cost.fxn, interval = interval, parmat[i,], maximum=TRUE)
    tvstar <- uroot$maximum ###Location
    vavgstar <- uroot$objective ###Maximum value

    parmat$tvstar[i] <- tvstar
    parmat$vavgstar[i] <- vavgstar

    vnullfun = splinefun(x = c(tvvals,365.1), c(vnullvals,vnullvals[1]))
    parmat$vnull[i]  <- 1/365*integrate(vnullfun, lower = 0, upper = 365)$value

    write.table(parmat[i,], file = filename, append = T, col.names = F, row.names = F)
    print(paste("NPar:", i, ' / ', NPars))
    return(i)
}

metaout <- mclapply(X = 1:NPars, FUN = mcl.fun, mc.cores = detectCores()-2)

refparmat <- read.table(file = filename, header = T)
refparmat <- refparmat[order(refparmat$NPar),]
write.table(refparmat, file = filename, append = F, row.names = F)
parmat <- read.table(file = filename, header = T)

