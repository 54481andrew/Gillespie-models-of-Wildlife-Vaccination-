SimName = "Sim_DeerMice"
parmat = read.table(file = paste("ParMat_", SimName, sep=''), header = F)
names(parmat) = c('Par','b0','d','Bp','Nv','tv','gamv','gamp','tb','T','IpInit', 'TPathInv')
NPars = nrow(parmat)
NTrials = 50

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
