SimName = "Sim_DeerMice"
parmat = read.table(file = paste("ParMat_", SimName, sep=''), header = F)
parmat$R0approx = with(parmat, Bp*(b0*tb)/(T*d*(d+gamp)))
parmat$R0star = with(parmat, R0approx*(1-Nv/(b0*tb + Nv*exp(-d*(T-tv)))))
names(parmat) = c('Par','b0','d','Bp','Nv','tv','gamv','gamp','tb','T','IpInit', 'TPathInv')
NPars = nrow(parmat)
NTrials = 50

TExtMatFile = paste('TExtMat_',SimName, sep = '')
TExtMat = read.table(TExtMatFile, header = FALSE)


#####################################
######## EXTINCTION PLOT C ##########
#####################################
##Contour plot. Y-axis: TPathInv
##X-axis: tv; Z-axis: F_Ext?
XValName = 'tv'
YValName = 'TPathInv'
YVals = unique(parmat[,YValName])
FixValName1 = 'Bp'
FixVals1 = c(0.00005)
FixValName2 = 'IpInit'
FixVals2 = c(1)

FigFold = paste(SimName,'_Fig',sep='')
if(!dir.exists(FigFold)){dir.create(FigFold)}

nXVals <- length(unique(parmat[,XValName]))
nYVals <- length(unique(parmat[,YValName]))
nFixVals1 <- length(FixVals1)
nFixVals2 <- length(FixVals2)

TCrit <- 3*365 

for(i1 in 1:nFixVals1){
    for(i2 in 1:nFixVals2){
        PExtMat = matrix(ncol = nYVals, nrow = nXVals)
        for(Yi in 1:nYVals){
            for(Xi in 1:nXVals){
                XVal = XVals[Xi]
                YVal = YVals[Yi]
                F1 = FixVals1[i1]
                F2 = FixVals2[i2]
                wifix = parmat[,XValName]==XVal & parmat[,YValName]==YVal & parmat[,FixValName1]==F1 & parmat[,FixValName2]==F2 
                PExtMat[Xi,Yi] = sum(TExtMat[wifix,] > TCrit+YVal)/NTrials	                       
            }
	}	
            FileName = paste('PExt_',TCrit,FixValName1, F1, FixValName2, F2,'.png',sep='')      
            #png(file = paste( FigFold, '/', FileName, sep=''), height = 5, width = 5, units = 'in', res = 400)
            #par(mai = c(1,1,0.1,0.1))
#            contour(PExtMat, xlab = '', ylab = '', ylim = range(YVals), xlim = range(XVals))

#	    contour(PExtMat)
filled.contour(PExtMat)
            #dev.off()
    }
}#End Loop through fixed vals




