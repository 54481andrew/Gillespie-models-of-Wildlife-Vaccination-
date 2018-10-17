round2 = function(x, n) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5
    z = trunc(z)
    z = z/10^n
    z*posneg
}

SimNameVals = c("B_1")
FoldNames = c("Data/")

for(FoldName in FoldNames){
for(SimName in SimNameVals){
    print('*********************')

    FileName <- paste(FoldName, SimName,sep='')
    print(paste("Working on ", FileName))
    
    parmat = read.table(file = paste(FileName,'/ParMat',sep=''), header = F)
    names(parmat) = c('Par','b0','d','Bp','Nv','tv','gamv','gamp','tb','T','IpInit', 'TPathInv')

    if(grepl('Freq',SimName)){
	parmat$R0p = with(parmat, round( Bp/(d+gamp) ,2))
    }else{
	parmat$R0p = with(parmat, round( Bp*(b0*tb)/(T*d*(d+gamp)) ,2))
    }

    parmat$R0approx = with(parmat, R0p*(1-Nv/(b0*tb + Nv*exp(-d*(T-tv)))))
    parmat$rho = with(parmat, round(Nv/(b0*tb/(d*T)),2))
    
    NPars = nrow(parmat)
    
    print('Reading TExtMat...')
    TExtMatFile = paste(FileName,'/TExtMat', sep = '')
    TExtMat = read.table(TExtMatFile, header = FALSE)
    
    print('Adding columns to parmat')
    refparmat = parmat
    pcritnamepre = 'pc'
    pcritvals <- c(0.99)
    ii <- 16 #New rows PExt
    for(pcrit in pcritvals){
        parmat[,ii] <- NA #Create column
        names(parmat)[ii] <- paste(pcritnamepre,pcrit,sep='')#name it
        ii <- ii + 1
    }

    TExtMat <- round2(TExtMat,2)
    nrowmat <- nrow(TExtMat)
    for(i in 1:nrowmat){
        tp <- parmat$TPathInv[i]
###Which trials had pathogen until time VaccStartTime       
        wiTrialsToVacc = 1:ncol(TExtMat)

###Loop through pcritvals
        for(pcrit in pcritvals){
	    TCrit <- with(parmat[i,], -1/(d+gamp)*log(1-pcrit))
            pcritname <- paste(pcritnamepre,pcrit,sep='')
            if(length(wiTrialsToVacc)==0)
            {
                refparmat[i,pcritname] <- NA
            }else{
                NTrialsToVacc = length(wiTrialsToVacc)
                refparmat[i,pcritname] <- sum(TExtMat[i,wiTrialsToVacc] <  
                                              round2(TCrit + tp, 2))/max(NTrialsToVacc,1)
            }
        }###Loop through 
        if(i%%1000==0){print(paste('i = ', i, '/', nrowmat))}
    
    }#Loop through parmat
    print('Writing refparmat')
    write.table(refparmat,file = paste(FileName,'/Ref_ParMat',sep=''), row.names = FALSE, col.names = TRUE)
}}
