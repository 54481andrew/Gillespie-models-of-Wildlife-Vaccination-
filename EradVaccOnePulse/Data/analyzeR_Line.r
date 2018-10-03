#Create a lineplot of extinction data

SimName = "DeerMice_Base_LP"
parmat = read.table(file = paste(SimName,"/ParMat", sep=''), header = F)
names(parmat) = c('Par','b0','d','Bp','Nv','tv','gamv','gamp','tb','T','IpInit', 'TPathInv')

if(grepl('Freq',SimName)){
	parmat$R0p = with(parmat, round( Bp/(d+gamp) ,2))
}else{
	parmat$R0p = with(parmat, round( Bp*(b0*tb)/(T*d*(d+gamp)) ,2))
}

parmat$R0approx = with(parmat, R0p*(1-Nv/(b0*tb + Nv*exp(-d*(T-tv)))))
parmat$rho = with(parmat, round(Nv/(b0*tb/(d*T)),2))

NPars = nrow(parmat)
NTrials = 1000
VaccStartTime = 8*365 #Time at which vaccination is started

TExtMatFile = paste(SimName,'/TExtMat', sep = '')
TExtMat = read.table(TExtMatFile, header = FALSE)


#####################################
######## EXTINCTION PLOT C ##########
#####################################
##Contour plot. Y-axis: TPathInv
##X-axis: tv;
XValName = 'tv'
XVals = unique(parmat[,XValName])
YValName = 'R0p'
YVals = c(1.1, 1.53, 1.97, 2.4)
#YVals = c(1.24, 1.48, 2.03, 2.5)#unique(parmat[,YValName])
ZValName = 'ProbExt'

FixValName1 = 'rho'
FixVals1 = unique(parmat[,FixValName1])#c(0.00001,0.00005,0.0001)
FixValName2 = 'tb'
FixVals2 = unique(parmat[,FixValName2])
FixValName3 = 'gamp'
FixVals3 = unique(parmat[,FixValName3])

FigFold = paste(SimName,'_Fig',sep='')
if(!dir.exists(FigFold)){dir.create(FigFold)}

nXVals <- length(unique(parmat[,XValName]))
nYVals <- length(YVals)#length(unique(parmat[,YValName]))
nFixVals1 <- length(FixVals1)
nFixVals2 <- length(FixVals2)
nFixVals3 <- length(FixVals3)

round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}
TCrit <- 365*1.5 #This script finds how many sim's made it time TCrit past the 1st pulse vaccination
TExtMat <- round2(TExtMat,2)

require(RColorBrewer)

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
		wiTrialsToVacc = which(TExtMat[wifix,] > round2(VaccStartTime + XVal,2))
		NTrialsToVacc = length(wiTrialsToVacc)
                PExtMat[Xi,Yi] = sum(TExtMat[wifix,wiTrialsToVacc] <  
			       (VaccStartTime + XVal + TCrit))/max(NTrialsToVacc,1)
	    }}#End loops through XVals and YVals


zmin = floor(10*min(PExtMat[!is.na(PExtMat)]))/10
zmax = 1
nbreaks = length(YVals)
initcols = c('purple', 'orange')#brewer.pal(, "Blues")
cols = colorRampPalette(initcols)(nbreaks)

axistext <- function(GraphMe){
	switch(GraphMe, 
	'Nv' = 'Number Vaccines', 
	'tv' = 'Time of Vaccination',
	'rho' = 'Scaled Number Vaccines',
	'R0p' = 'R0p')
}

if(zmin < zmax){
            FileName = paste('PExtLine_',TCrit,FixValName1, F1, FixValName2, F2,FixValName3,F3,'.png',sep='') 
            png(file = paste( FigFold, '/', FileName, sep=''), height = 4, width = 5, units = 'in', 
	    	     res = 400)
            par(mai = c(1,1,0.25,0.25))
	    PExtMat[is.na(PExtMat)] <- -1
	    matplot(x = XVals, y = PExtMat, col = cols, 
		    xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', ylim = c(0,1.5), 
		    type = 'l', lwd = 2)
            axislabs = seq(0,365, by = 60)
            axislabs1 = YVals
	    axis(side = 1, labels = T, at = axislabs)
	    axis(side = 2, labels = T)
	    mtext(text = axistext(XValName),side = 1, line = 3)
	    mtext(text = "Probability of Extinction",side = 2, line = 3)
	    abline(v = parmat$tb[wifix], lwd = 3, lty = 3)
	    legend(x = 'topright', legend = paste(axistext(YValName),'=',round(YVals,1)), col = cols, lwd = 2) 
           dev.off()

print(paste("Figure", i, "complete",sep = " " ))
}else{
print(paste("Figure", i, "fail",sep = " " ))
}
i=i+1
    }}}#End Loop through fixed vals




