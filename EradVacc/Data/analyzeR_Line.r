SimName = "DeerMice_Ha_1"
parmat = read.table(file = paste("ParMat_", SimName, sep=''), header = F)
names(parmat) = c('Par','b0','d','Bp','Nv','tv','gamv','gamp','tb','T','IpInit', 'TPathInv')
NPars = nrow(parmat)
NTrials = 200

TExtMatFile = paste('TExtMat_',SimName, sep = '')
TExtMat = read.table(TExtMatFile, header = FALSE)


#####################################
######## EXTINCTION PLOT L ##########
#####################################
##Line plot. Y-axis: proportion of NTrials that are extinct 
##after TFoc years. X-axis, XValName.  
XValName = 'tv'
XVals = unique(parmat[,XValName])
ZValName = '' #Parameter varied across line plots
ZVals = unique(parmat[,ZValName])
FixValName1 = 'Bp'
FixVals1 = c(0.005)
FixValName2 = 'Nv'
FixVals2 = unique(parmat[,FixValName2])

FigFold = paste(SimName,'_Fig',sep='')
if(!dir.exists(FigFold)){dir.create(FigFold)}

nXVals <- length(unique(parmat[,XValName]))
nZVals <- length(unique(parmat[,ZValName]))
nFixVals1 <- length(FixVals1)
nFixVals2 <- length(FixVals2)

TCrit <- 1*365

for(i1 in 1:nFixVals1){
    for(i2 in 1:nFixVals2){
        PExtMat = matrix(ncol = nZVals, nrow = nXVals)
        for(Yi in 1:nZVals){
            for(Xi in 1:nXVals){
                XVal = XVals[Xi]
                ZVal = ZVals[Yi]
                F1 = FixVals1[i1]
                F2 = FixVals2[i2]
                wifix = parmat[,XValName]==XVal & parmat[,ZValName]==ZVal & parmat[,FixValName1]==F1 & 
		      parmat[,FixValName2]==F2 
                PExtMat[Xi,Yi] = sum(TExtMat[wifix,] > TCrit)/NTrials	                       
            }
            FileName = paste('PExt_',TCrit,FixValName1, F1, FixValName2, F2,'.png',sep='')      
            png(file = paste( FigFold, '/', FileName, sep=''), height = 5, width = 5, units = 'in', res = 400)
            par(mai = c(1,1,0.1,0.1))
            matplot(PExtMat, xlab = '', ylab = '')
            dev.off()
        }
    }
}#End Loop through fixed vals

