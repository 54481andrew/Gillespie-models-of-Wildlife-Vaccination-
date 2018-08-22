dat <- read.table(file = 'test', header = TRUE)


tiff(file = 'SIR.tif', width = 4, height = 4, units = 'in', res = 400)
plot(NA, xlim = c(0,100), ylim = c(0,max(dat[,c('S','Ip','Rp')])), xlab = 'Time', ylab = 'Number', main = "Endemic SIR model")
points(S~time, dat, col = 'blue')
points(Ip~time, dat, col = 'red')
points(Rp~time, dat, col = 'pink')
points(N~time, dat, col = 'black')

###Plot abundances from SIR model 
b <- 500
R0p <- 10
delp <- 0.1
d <- 0.5
Bp <- R0p*d*(d + delp) / b;

###Equilibrium values predicted by continuous-time SIR model
Sstar  <- b/d*1/R0p
Ipstar <-b/(d+delp)*(1-1/R0p)
Rpstar <-delp/d*Ipstar

abline(h = Sstar, lwd = 3, col = 'blue')
abline(h = Ipstar, lwd = 3, col = 'red')
abline(h = Rpstar, lwd = 3, col = 'pink')
legend(x = 'right', legend = c('S', 'Ip', 'Rp'), col = c("blue", "red", "pink"), lwd = 2, cex = 1)
dev.off()
