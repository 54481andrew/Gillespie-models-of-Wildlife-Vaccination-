dat = read.table(file = '1', header = T, sep = " ")
names(dat)


rhs.fun = function(t,y,parms){
    S <- y[1]
    Iv <- y[2]
    Ip <- y[3]
    V <- y[4]
    P <- y[5]
    with(parms,{
        dS = b0*((t%%Tb) < tb) - Bp*S*Ip - d*S
        dIv = -(gamv + d)*Iv
        dIp = Bp*S*Ip - (gamp + d)*Ip
        dV = gamv*Iv - d*V
        dP = gamp*Ip - d*P
        
        return(list(c(dS,dIv,dIp,dV,dP)))
    })
}

require(deSolve)
parms = data.frame(b0 = 500, tb = 30, d = 0.001, Tb = 365, Bp = 0.00001, gamp = 0.001, gamv = 0.07)
times = seq(0,10000)
y0 = c(1000, 0, 100, 0, 0)
names(y0) = c('S','Iv','Ip','V','P')
out = data.frame(lsoda(y = y0, times = times, func = rhs.fun, parms = parms))



par(mfrow = c(1,2))
plot(S~time, out, xlim = c(0,365), ylim = c(0,20000), type = 'l', lwd = 2)
lines(Iv~time, out, col = 'green')
lines(Ip~time, out, col = 'red')
lines(V~time, out, col = 'blue')
lines(P~time, out, col = 'pink')

plot(S~time, dat, col = 'black', cex = 0.5, xlim = c(0, 365), ylim = c(0,25000))
points(Iv~time, dat, col = 'green', cex = 0.5)
points(Ip~time, dat, col = 'red', cex = 0.5)
points(V~time, dat, col = 'blue', cex = 0.5)
points(P~time, dat, col = 'pink', cex = 0.5)

par(mfrow = c(1,1))
pdf(file = "Comp.pdf")
plot(S~time, out, xlim = c(0,365*5), ylim = c(0,30000), type = 'l', lwd = 4)
lines(Iv~time, out, col = 'green', lty = 2, lwd = 4)
lines(Ip~time, out, col = 'red', lty = 2, lwd = 4)
lines(V~time, out, col = 'blue', lty = 2, lwd = 4)
lines(P~time, out, col = 'pink', lty = 2, lwd = 4)

these <- seq(1,nrow(dat), by = 10)
points(S~time, dat[these,], col = 'black', cex = 1)
points(Iv~time, dat[these,], col = 'green', cex = 1)
points(Ip~time, dat[these,], col = 'red', cex = 1)
points(V~time, dat[these,], col = 'blue', cex = 1)
points(P~time, dat[these,], col = 'pink', cex = 1)

dev.off()
