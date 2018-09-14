SimName = "Test_DeerMice_Base"
parmat = read.table(file = paste("ParMat_", SimName, sep=''), header = F)
names(parmat) = c('Par','b0','d','Bp','Nv','tv','gamv','gamp','tb','T','IpInit', 'TPathInv')
parmat$R0approx = with(parmat, Bp*(b0*tb)/(T*d*(d+gamp)))
parmat$R0star = with(parmat, R0approx*(1-Nv/(b0*tb + Nv*exp(-d*(T-tv)))))
NPars = nrow(parmat)
NTrials = 1000
VaccStartTime = 8*365 #Time at which vaccination is started

TExtMatFile = paste('TExtMat_',SimName, sep = '')
TExtMat = read.table(TExtMatFile, header = FALSE)


#####################################
######## EXTINCTION PLOT C ##########
#####################################
##Contour plot. Y-axis: TPathInv
##X-axis: tv;
XValName = 'tv'
XVals = unique(parmat[,XValName])
YValName = 'Nv'
YVals = unique(parmat[,YValName])
FixValName1 = 'Bp'
FixVals1 = unique(parmat[,FixValName1])#c(0.00001,0.00005,0.0001)
FixValName2 = 'b0'
FixVals2 = unique(parmat[,FixValName2])
FixValName3 = 'd'
FixVals3 = 0.004

FigFold = paste(SimName,'_Fig',sep='')
if(!dir.exists(FigFold)){dir.create(FigFold)}

nXVals <- length(unique(parmat[,XValName]))
nYVals <- length(unique(parmat[,YValName]))
nFixVals1 <- length(FixVals1)
nFixVals2 <- length(FixVals2)
nFixVals3 <- length(FixVals3)

TCrit <- 365*1 #This script finds how many sim's made it time TCrit past the 1st pulse vaccination

require(RColorBrewer)

ii = 1
vaccset = c()
i = 1
for(i1 in 1:nFixVals1){
    for(i2 in 1:nFixVals2){
    for(i3 in 1:nFixVals3){
        PExtMat = matrix(ncol = nYVals, nrow = nXVals)
        for(Yi in 1:nYVals){
            for(Xi in 1:nXVals){
                XVal = XVals[Xi]
                YVal = YVals[Yi]
                F1 = FixVals1[i1]
                F2 = FixVals2[i2]
                F3 = FixVals3[i3]
                wifix = parmat[,XValName]==XVal & parmat[,YValName]==YVal & parmat[,FixValName1]==F1 & 
		      parmat[,FixValName2]==F2 & parmat[,FixValName3]==F3 

		#Which trials had pathogen until time VaccStartTime
		wiTrialsToVacc = which(TExtMat[wifix,] > (VaccStartTime + XVal))
		NTrialsToVacc = length(wiTrialsToVacc)
                PExtMat[Xi,Yi] = sum(TExtMat[wifix,wiTrialsToVacc] > 
			       (VaccStartTime + XVal + TCrit))/NTrialsToVacc
	       vaccset[ii] = NTrialsToVacc 
			       ii = ii + 1
	    }}#End loops through XVals and YVals

zmin = floor(10*min(PExtMat))/10
zmax = 1
breaks = seq(zmin,zmax,by = 0.05)
nbreaks = length(breaks)
initcols = c('purple', 'orange')#brewer.pal(, "Blues")
cols = colorRampPalette(initcols)(nbreaks-1)

axistext <- function(GraphMe){
	switch(GraphMe, 
	'Nv' = 'Number Vaccines', 
	'tv' = 'Time of Vaccination')   
}

if(zmin < zmax){
            FileName = paste('PExt_',TCrit,FixValName1, F1, FixValName2, F2,FixValName3,F3,'.png',sep='') 
            png(file = paste( FigFold, '/', FileName, sep=''), height = 4, width = 5, units = 'in', 
	    	     res = 400)
            par(mai = c(1,1,0.25,0.25))
	    layout(mat = matrix(c(1,2),ncol = 2), widths = c(1,0.2))
	    
	    image(x = XVals, y = YVals, z = PExtMat, col = cols, breaks = breaks, 
		    xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
            axislabs = seq(0,365, by = 60)
            axislabs1 = YVals
	    axis(side = 1, labels = T, at = axislabs)
	    axis(side = 2, labels = T, at = seq(0,5000,length.out = 21))
	    mtext(text = axistext(XValName),side = 1, line = 3)
	    mtext(text = axistext(YValName),side = 2, line = 3)

	    #Build legend bar
	    Zmat = matrix(breaks, ncol = nbreaks)
	    ZmatCent = matrix(0.5*(Zmat[-1] + Zmat[-length(Zmat)]), ncol = nbreaks-1)

            par(mai = c(1,0.05,0.25,0.6))
	    image(x = 1, y = Zmat, z = ZmatCent, col = cols, breaks = breaks, 
	    	    ylim = range(breaks), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
	    axis(side = 4, las = 1)
            dev.off()

print(paste("Figure", i, "complete",sep = " " ))
}else{
print(paste("Figure", i, "fail",sep = " " ))
}
i=i+1
    }}}#End Loop through fixed vals




