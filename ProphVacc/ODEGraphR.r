SimName = "LagImport.TMax40.1yr"

filename = paste("ODEData/", SimName, sep = '')
parmat <- read.table(file = filename, header = T)

XName <- 'tb'
YName <- 'rho'
XVals  <- unique(parmat[,XName])
YVals  <- unique(parmat[,YName])
nXVals <- length(XVals)
nYVals <- length(YVals)
FixName1 <- 'd'
FixVals1 <- unique(parmat[,FixName1])
nFixVals <- length(FixVals1)

FigFold = paste(SimName,'_Fig',sep='')
if(!dir.exists(FigFold)){dir.create(FigFold)}

for(F1 in 1:nFixVals){
FixVal <- FixVals1[F1]
Import.Mat = matrix(NA,nrow = nXVals, ncol = nYVals)
Lag.Mat = matrix(NA,nrow = nXVals, ncol = nYVals)
for(xi in 1:nXVals){
    for(yi in 1:nYVals){
        wi <- XVals[xi]==parmat[,XName] & YVals[yi]==parmat[,YName] & 
	      parmat[,FixName1]==FixVal
        TvStar <- parmat$tvstar[wi]
        VAvgStar <- parmat$vavgstar[wi]
        VNull <-  parmat$vnull[wi]
	Lag <- (TvStar - unique(parmat$tb[wi]))
	
	#if(Lag > 182.5){Lag <- Lag - 365}
	#if(Lag < -182.5){Lag <- Lag + 365}
        Lag.Mat[xi,yi] <- Lag
        Import.Mat[xi,yi] <- VAvgStar/VNull
    }
}

Import.Mat <- Import.Mat - 1

zmin.lag = floor(10*min(Lag.Mat))/10
zmax.lag = ceiling(10*max(Lag.Mat))/10
zmin.imp = 0
zmax.imp = ceiling(10*max(Import.Mat))/10
breaks.lag = seq(-185,400,by = 25)
breaks.imp = seq(zmin.imp, zmax.imp, by = 0.05)	
nbreaks.lag = length(breaks.lag)
nbreaks.imp = length(breaks.imp)
initcols = c('purple', 'orange')
cols.lag = colorRampPalette(initcols)(nbreaks.lag-1)
cols.imp = colorRampPalette(initcols)(nbreaks.imp-1)

            FileName = paste(FixName1, FixVals1[F1], '.png',sep='')      
            png(file = paste( FigFold, '/', FileName, sep=''), height = 3.25, width = 6, units = 'in', res = 400)
	    layout(mat = matrix(c(1,2),ncol = 2), widths = c(1,1))
            par(mai = c(0.75,0.15,0.15,0.15), omi = c(0,0.75,0,0))
	    image(x = XVals, y = YVals, z = Lag.Mat, col = cols.lag, breaks = breaks.lag, 
	    	    xlab = '', ylab = '',
		    xaxt = 'n', yaxt = 'n', xlim = c(0,365), ylim = c(0,1))
	    contour(x = XVals, y = YVals, z = Lag.Mat, levels = breaks.lag, add = T, 
	    	      labcex = 1)
	    axis(side = 1, labels = T, at = seq(0,365, by = 120))
	    axis(side = 2, labels = T, at = seq(0,1, by = 0.1))

	    mtext(side = 1, outer = F, text = 'Breeding Duration (days)', line = 2.25)
	    mtext(side = 2, outer = F, text = 'Vaccines per Host', line = 2.5)
###############
            par(mai = c(0.75,0.15,0.15,0.15))	    
	    image(x = XVals, y = YVals, z = Import.Mat, col = cols.imp, 
	    breaks = breaks.imp, xlab = '', ylab = '',
		    xaxt = 'n', yaxt = 'n', xlim = c(0,365), ylim = c(0,1))
	    contour(x = XVals, y = YVals, z = Import.Mat, levels = breaks.imp, 
	    	      add = T, labcex = 1)
	    axis(side = 1, labels = T, at = seq(0,365, by = 120))
	    axis(side = 2, labels = F, at = seq(0,1, by = 0.1))

	    mtext(side = 1, outer = F, text = 'Breeding Duration (days)', line = 2.25)

###############

dev.off()



###################
#### FIGURE 2 #####
###################

	    these <- c(20,30,40,50)#rho values to plot
            FileName = paste(FixName1, FixVals1[F1], 'Line.png',sep='')      
            png(file = paste( FigFold, '/', FileName, sep=''), height = 4, width = 5, units = 'in', res = 400)
            par(mai = c(0.75,1,0.15,0.15))	    
	    matplot(x = XVals, y = Lag.Mat[,these],
	    	    xlab = '', ylab = '',
		    xaxt = 'n', yaxt = 'n', type = 'l', ylim = c(-50, 50))
	    axis(side = 1, labels = T, at = seq(0,365, by = 120))
	    axis(side = 2, labels = T, at = seq(-182.5,182.5, by = 10))

	    mtext(side = 1, outer = F, text = 'Breeding Duration (days)', line = 2.25)

	    dev.off()











}