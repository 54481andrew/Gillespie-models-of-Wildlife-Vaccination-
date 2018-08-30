SimName = "Sim_1"
parmat = read.table(file = paste("ParMat_", SimName, sep=''), header = T)
NPars = nrow(parmat)
NTrials = 50

#Loop through parameters, read data into matrix. 
#For now, just find the extinction time for each
#trial, of each parameter set. 
ExtTimeMat = matrix(nrow = NTrials, ncol = NPars)
for(j in 1:NPars){
    print(paste("NPar:", j-1))
    filename = paste(SimName, "/NPar_", j-1, sep = '')
    dat = read.table(file = filename, header = T, sep = " ")
    wiStarts <- c(which(dat$time==0), nrow(dat))
    wiEnds   <- wiStarts[-1] - 1
    ExtTimeMat[,j] = dat$time[wiEnds] - parmat$TPathInv[j]
}


#####################################
######## SIMPLE HISTOGRAM ###########
#####################################
XValName = 'tv'
FixValName1 = 'PInit'
FixVals1 = c(1)
FixValName2 = 'TPathInv'
FixVals2 = c(30)
FixValName3 = 'Bp'
FixVals3 = c(0.00005)

FigFold = paste(SimName,'_Fig',sep='')
if(!dir.exists(FigFold)){dir.create(FigFold)}


nFixVals1 <- length(FixVals1)
nFixVals2 <- length(FixVals2)
nFixVals3 <- length(FixVals3)

for(i1 in 1:nFixVals1){
for(i2 in 1:nFixVals2){
for(i3 in 1:nFixVals3){

       F1 = FixVals1[i1]
       F2 = FixVals2[i2]
       F3 = FixVals3[i3]
       wifix = parmat[,FixValName1]==F1 && parmat[,FixValName2]==F2 && 
       	       parmat[,FixValName3]==F3 
       HistVals = ExtTimeMat[,wifix]	       
 
    FileName = paste(FixValName1, F1, FixValName2, F2, FixValName3,'_', F3,'.png',sep='')      
    png(file = paste( FigFold, '/', FileName, sep=''), height = 5, width = 5, units = 'in', res = 400)
    par(mai = c(1,1,0.1,0.1))
    hist(HistVals, freq = F, xlab = '', ylab = '')
    #legend = c('S','Iv','Ip','V','P')
    #legend(x = "topright", legend = legend, col = c('black', 'green', 'red', 'blue', 'pink'), lwd = 2)
    dev.off()

}}}#End Loop through fixed vals

#####################################
######## EXTINCTION PLOT L ##########
#####################################
##Line plot. Y-axis: proportion of NTrials that are extinct 
##after TFoc years. X-axis, XValName.  
XValName = 'tv'
YValName = 'PInit'
YVals = unique(parmat[,YValName])
FixValName1 = 'Bp'
FixVals1 = c(0.00005)
FixValName2 = 'TPathInv'
FixVals2 = c(30)

FigFold = paste(SimName,'_Fig',sep='')
if(!dir.exists(FigFold)){dir.create(FigFold)}

nXVals <- length(unique(parmat[,XValName]))
nYVals <- length(unique(parmat[,YValName]))
nFixVals1 <- length(FixVals1)
nFixVals2 <- length(FixVals2)

TCrit <- 365

for(i1 in 1:nFixVals1){
    for(i2 in 1:nFixVals2){
        PExtMat = matrix(ncol = nYVals, nrow = nXVals)
        for(Yi in 1:nYVals){
            for(Xi in 1:nXVals){
                
                YVal = YVals[Yi]
                F1 = FixVals1[i1]
                F2 = FixVals2[i2]
                wifix = parmat[,XValName]==XVal && parmat[,YValName]==YVal && parmat[,FixValName1]==F1 && parmat[,FixValName2]==F2 
                PExtMat[Xi,Yi] = sum(ExtTimeMat[,wifix] > TCrit)/NTrials	                       
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
######## EXTINCTION PLOT C ##########NOT DONE
#####################################
##Contour plot. Y-axis: TPathInv
##X-axis: tv; Z-axis: F_Ext?
XValName = 'tv'
YValName = 'TPathInv'
YVals = unique(parmat[,YValName])
FixValName1 = 'Bp'
FixVals1 = c(0.00005)
FixValName2 = 'PInit'
FixVals2 = c(10)

FigFold = paste(SimName,'_Fig',sep='')
if(!dir.exists(FigFold)){dir.create(FigFold)}

nXVals <- length(unique(parmat[,XValName]))
nYVals <- length(unique(parmat[,YValName]))
nFixVals1 <- length(FixVals1)
nFixVals2 <- length(FixVals2)

TCrit <- 365

for(i1 in 1:nFixVals1){
    for(i2 in 1:nFixVals2){
        PExtMat = matrix(ncol = nYVals, nrow = nXVals)
        for(Yi in 1:nYVals){
            for(Xi in 1:nXVals){
                
                YVal = YVals[Yi]
                F1 = FixVals1[i1]
                F2 = FixVals2[i2]
                wifix = parmat[,XValName]==XVal && parmat[,YValName]==YVal && parmat[,FixValName1]==F1 && parmat[,FixValName2]==F2 
                PExtMat[Xi,Yi] = sum(ExtTimeMat[,wifix] > TCrit)/NTrials	                       
            }
            FileName = paste('PExt_',TCrit,FixValName1, F1, FixValName2, F2,'.png',sep='')      
            png(file = paste( FigFold, '/', FileName, sep=''), height = 5, width = 5, units = 'in', res = 400)
            par(mai = c(1,1,0.1,0.1))
            contour(PExtMat, xlab = '', ylab = '')
            dev.off()
        }
    }
}#End Loop through fixed vals




