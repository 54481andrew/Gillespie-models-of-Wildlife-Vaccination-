SimName = "Sim_DeerMice"
parmat = read.table(file = paste("ParMat_", SimName, sep=''), header = F)
names(parmat) = c('Par','b0','d','Bp','Nv','tv','gamv','gamp','tb','T','IpInit', 'TPathInv')
NPars = nrow(parmat)
NTrials = 50

#Loop through parameters, read data into matrix. 
#For now, just find the extinction time for each
#trial, of each parameter set. 
#ExtTimeMat = matrix(nrow = NTrials, ncol = NPars)
#for(j in 1:NPars){
#    print(paste("NPar:", j-1))
#    filename = paste(SimName, "/NPar_", j-1, sep = '')
#    dat = read.table(file = filename, header = T, sep = " ")
#    wiStarts <- c(which(dat$time==0), nrow(dat))
#    wiEnds   <- wiStarts[-1] - 1
#    ExtTimeMat[,j] = dat$time[wiEnds] - parmat$TPathInv[j]
#}
TExtMatFile = paste('TExtMat_',SimName, sep = '')
TExtMat = read.table(TExtMatFile, header = FALSE)

#####################################
######## SIMPLE HISTOGRAM ###########
#####################################
XValName = 'tv'
FixValName1 = 'IpInit'
FixVals1 = c(1)
FixValName2 = 'TPathInv'
parmat[,FixValName2] = round(parmat[,FixValName2]%%365,2)
FixVals2 = 29.2
FixValName3 = 'Bp'
FixVals3 = c(0.00005)
FixValName4 = 'tv'
parmat[,FixValName4] = round(parmat[,FixValName4],2)
FixVals4 = c(1,59.24, 117.48, 175.72, 233.96, 306.76)

FigFold = paste(SimName,'_Fig',sep='')
if(!dir.exists(FigFold)){dir.create(FigFold)}


nFixVals1 <- length(FixVals1)
nFixVals2 <- length(FixVals2)
nFixVals3 <- length(FixVals3)
nFixVals4 <- length(FixVals4)

breaks = seq(0,12*365, 12)

FileName = paste(FixValName1, F1, FixValName2, F2, FixValName3, F3,FixValName4,F4,'.png',sep='')      
png(file = paste( FigFold, '/', FileName, sep=''), height = 5, width = 5, units = 'in', res = 400)
par(mai = c(1,1,0.1,0.1))
plot(NA, xlim = c(0*365, 11*365), ylim = c(0,0.15))

cols = c('black', 'green', 'red', 'blue', 'pink')

for(i1 in 1:nFixVals1){
for(i2 in 1:nFixVals2){
for(i3 in 1:nFixVals3){
for(i4 in 1:nFixVals4){

#if(i1*i2*i3*i4==1){
#    hist(unlist(HistVals), freq = F, xlab = '', ylab = '', add = F, breaks = breaks, ylim = c(0,1))
#}

       F1 = FixVals1[i1]
       F2 = FixVals2[i2]
       F3 = FixVals3[i3]
       F4 = FixVals4[i4]
       wifix = parmat[,FixValName1]==F1 & parmat[,FixValName2]==F2 &
       	       parmat[,FixValName3]==F3 & parmat[,FixValName4]==F4
       HistVals = TExtMat[wifix, ]	       
 
    hist(unlist(HistVals), freq = F, xlab = '', ylab = '', add = T, breaks = breaks, col = cols[i4])
    #legend = c('S','Iv','Ip','V','P')
    #legend(x = "topright", legend = legend, col = c('black', 'green', 'red', 'blue', 'pink'), lwd = 2)
}
}}}#End Loop through fixed vals
dev.off()


#####################################
######## EXTINCTION PLOT L ##########
#####################################
##Line plot. Y-axis: proportion of NTrials that are extinct 
##after TFoc years. X-axis, XValName.  
XValName = 'tv'
XVals = unique(parmat[,XValName])
YValName = 'IpInit'
YVals = unique(parmat[,YValName])
FixValName1 = 'Bp'
FixVals1 = c(0.00005)
FixValName2 = 'TPathInv'
FixVals2 = unique(parmat[,FixValName2])
FixVals2 = c(29.2)

FigFold = paste(SimName,'_Fig',sep='')
if(!dir.exists(FigFold)){dir.create(FigFold)}

nXVals <- length(unique(parmat[,XValName]))
nYVals <- length(unique(parmat[,YValName]))
nFixVals1 <- length(FixVals1)
nFixVals2 <- length(FixVals2)

TCrit <- 2*365

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




