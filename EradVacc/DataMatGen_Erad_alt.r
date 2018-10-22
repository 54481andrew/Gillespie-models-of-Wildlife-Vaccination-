round2 = function(x, n) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5
    z = trunc(z)
    z = z/10^n
    z*posneg
}

SimNameVals = c("D_Freq")
FoldNames = c("Data/")

for(FoldName in FoldNames){
for(SimName in SimNameVals){
    print('*********************')

    FileName <- paste(FoldName, SimName,sep='')
    print(paste("Working on ", FileName))
    
    parmat = read.table(file = paste(FileName,'/ParMat',sep=''), header = F)
    names(parmat) = c('Par','b0','d','Bp','Nv','tv','gamv','gamp','tb','T','IpInit', 
    'TVacc')

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
    tcritnamepre = 'tc'
    pcritvals <- c(0.99)
    tcritvals  <- ceiling(c(0.5,1,2,3,4)*365)
    ii <- 16 #New rows PExt

    for(tcrit in tcritvals){
        parmat[,ii] <- NA #Create column
        names(parmat)[ii] <- paste(tcritnamepre,tcrit,sep='')#name it
        ii <- ii + 1
    }

    TExtMat <- round2(TExtMat,2)
    nrowmat <- nrow(TExtMat)
    for(i in 1:nrowmat){
        tvacc <- parmat$TVacc[i] + parmat$tv[i]
###Which trials had pathogen until time VaccStartTime       
        wiTrialsToVacc = TExtMat[i,] > round2(tvacc, 2)

###Loop through pcritvals
        for(tcrit in tcritvals){
	    TCrit <- tcrit  #with(parmat[i,], -1/(d+gamp)*log(1-pcrit))
            tcritname <- paste(tcritnamepre,tcrit,sep='')
            if(length(wiTrialsToVacc)==0)
            {
                refparmat[i,tcritname] <- NA
            }else{
                NTrialsToVacc = length(wiTrialsToVacc)
                refparmat[i,tcritname] <- sum(TExtMat[i,wiTrialsToVacc] <  
                                              round2(TCrit + tvacc, 2))/max(NTrialsToVacc,1)
            }
        }###Loop through 
        if(i%%1000==0){print(paste('i = ', i, '/', nrowmat))}
    
    }#Loop through parmat
    print('Writing refparmat')
    write.table(refparmat,file = paste(FileName,'/Ref_ParMat',sep=''), row.names = FALSE, col.names = TRUE)
}}
