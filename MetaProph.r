BuildDataFlag = FALSE
FoldName = "MetaProphFreq_Fig"
SimType = "ProphVacc"
#SimNameList = paste(SimType,"/Data/DeerMice_Base_varNv_redo",sep='')
SimNameList = paste(SimType,"/Data/DeerMice_Base_Freq_varNv",sep='')

#Collect all data
#for(i in 1:length(SimNameList)){

ii = 1

SimName = SimNameList[ii]
parmat = read.table(file = paste(SimName,"/ParMat", sep=''), header = F)
names(parmat) = c('Par','b0','d','Bp','Nv','tv','gamv','gamp','tb','T','IpInit', 'TPathInv')
parmat$R0approx = with(parmat, Bp*(b0*tb)/(T*d*(d+gamp)))
parmat$R0star = with(parmat, R0approx*(1-Nv/(b0*tb + Nv*exp(-d*(T-tv)))))
NPars = nrow(parmat)
NTrials = 1000

parmat$NAvg = with(parmat, b0*tb/(d*T))
parmat$rho = with(parmat, Nv/(b0*tb/(d*T)))

#####################################
######## EXTINCTION PLOT C ##########
#####################################
##Contour plot. Y-axis: TPathInv
##X-axis: tv; Z-axis: F_Ext?
XValName = 'tv'
XVals = unique(parmat[,XValName])
YValName = 'TPathInv'
YVals = unique(parmat[,YValName])
FixValName1 = 'Bp'
FixVals1 = unique(parmat[,FixValName1])[1:2]  #c(0.00001,0.00005,0.0001)
FixValName2 = 'IpInit'
FixVals2 = 10#unique(parmat[,FixValName2])
FixValName3 = 'Nv'
FixVals3 = c(100,250,500)#unique(parmat[,FixValName3])

FigFold = paste(FoldName,sep='')
if(!dir.exists(FigFold)){dir.create(FigFold)}

nXVals <- length(unique(parmat[,XValName]))
nYVals <- length(unique(parmat[,YValName]))
nFixVals1 <- length(FixVals1)
nFixVals2 <- length(FixVals2)
nFixVals3 <- length(FixVals3)

TCrit <- 365*1 

require(RColorBrewer)
require(MASS)


if(BuildDataFlag){

TExtMatFile = paste(SimName,'/TExtMat', sep = '')
TExtMat = read.table(TExtMatFile, header = FALSE)


MetaDat = matrix(nrow = nFixVals1*nFixVals2*nFixVals3 + 1, ncol = length(XVals) + 3)
MetaDat[1,] = c(1,2,3,XVals)

i = 2
for(i1 in 1:nFixVals1){
    for(i2 in 1:nFixVals2){
    	   
    	   F1 = FixVals1[i1]
	   F2 = FixVals2[i2]

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
                PExtMat[Xi,Yi] = sum(TExtMat[wifix,] >= TCrit+YVal)/NTrials	                       
            }
	}

TPathAvg = rowMeans(PExtMat)
MetaDat[i,] = c(F1,F2,F3,TPathAvg)	    
i = i + 1

print(paste("Finished i3: ", i3, sep = " " ))

}     #End loop through FixedVal3

}}#End loop through FixedVals

write.matrix(MetaDat, file = paste(FigFold, '/', SimType,'_metadat.txt',sep=''))
}else{
	MetaDat = as.matrix(read.table(file = paste(FigFold,'/', SimType,'_metadat.txt',sep=''), header = FALSE))
}


###Graph MetaDat
nbreaks = length(FixVals3)
initcols = c('purple', 'orange')#brewer.pal(, "Blues")
cols = c('black','green','red')#colorRampPalette(initcols)(nbreaks-1)

for(i1 in 1:nFixVals1){
    for(i2 in 1:nFixVals2){
	  	   
    	   F1 = FixVals1[i1]
	   F2 = FixVals2[i2]
    	   wipar = MetaDat[,1]== F1 & MetaDat[,2]==F2

            FileName = paste('TPathAvg_',SimType, '_',FixValName1, F1, FixValName2, F2,'.png',sep='')      
            png(file = paste( FigFold, '/', FileName, sep=''), height = 4, width = 5, units = 'in', res = 400)
            par(mai = c(1,1,0.25,0.25))
	    
	    if(SimType=="ProphVacc"){xrange <- c(0,365)}else{xrange <- c(8*365, 9*365)}

	    matplot(matrix(MetaDat[1,-(1:3)], ncol=1),t(MetaDat[wipar,-(1:3)]), xlim = xrange, ylim = c(0,1), xaxt = 'n', yaxt = 'n', 
	    	     xlab = '', ylab = '', pch = 1, lwd = 2, cex = 1, col = cols, main = SimType)

#	    lines(x = XVals%%365, y = TPathAvg, col = cols[i3], lwd = 2)

            axislabs = seq(0,365, by = 60)
            axislabs1 = YVals
	    if(SimType=="ProphVacc"){axlabs = seq(0,365, by = 60)}else{axlabs = seq(8*365, 9*365, by = 60)}
	    axis(side = 1, labels = axislabs, at = axlabs)

	    axis(side = 2, labels = T, at = seq(0,1,by =0.1))

	    mtext(side = 1, text = 'Time of Vaccination', line = 3)
	    mtext(side = 2, text = 'Probability of Establishment', line = 3)
	    abline(v = parmat$tb[1] + 8*365, lwd = 3, lty = 3)
	    #legendtext = paste('Nv = ', FixVals3, sep='')
	    legendtext = sapply(round(FixVals3/(parmat$NAvg[1]),1), function(x) as.expression(substitute(rho == B,
                    list(B = as.name(x)))))
	    #legendtext = expression(bquote(rho ~ '= ' ~ .(FixVals3[1]/(parmat$NAvg[1]))))
	    legend(x = 'bottomright', legend = legendtext, pch = 1,pt.lwd = 2, col = cols)
            dev.off()

}}