round2 = function(x, n) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5
    z = trunc(z)
    z = z/10^n
    z*posneg
}
VaccStartTime = 5*365 #Time at which vaccination is started
SimNameVals = "A_Freq"
FoldNames = c("Data/")

for(FoldName in FoldNames){
for(SimName in SimNameVals){
    print('*********************')

    FileName <- paste(FoldName, SimName,sep='')
    print(paste("Working on ", FileName))
    
    parmat = read.table(file = paste(FileName,'/ParMat',sep=''), header = F)
    names(parmat) = c('Par','b0','d','Bp','Nv','tv','gamv','gamp','tb','T','IpInit', 'TPathInv', 'NPeak')

    if(grepl('Freq',SimName)){
	parmat$R0p = with(parmat, round( Bp/(d+gamp) ,2))
    }else{
	parmat$R0p = with(parmat, round( Bp*(b0*tb)/(T*d*(d+gamp)) ,2))
    }

    parmat$rho = with(parmat, round(Nv/NPeak,2))
    
    NPars = nrow(parmat)
    
    print('Reading TExtMat...')
    TExtMatFile = paste(FileName,'/TExtMat', sep = '')
    TExtMat = read.table(TExtMatFile, header = FALSE)
    
    print('Adding columns to parmat')
    refparmat = parmat
    tcritnamepre = 'time'
    tcritvals <- c(365)
    ii <- ncol(parmat) + 1 #New rows PExt
    for(tcrit in tcritvals){
        parmat[,ii] <- NA ##Num trials that survived to vacc
	parmat[,ii+1] <- NA ##Num trials eradicated
	parmat[,ii+2] <- NA ##Proportion trials eradicated
	n1 <- paste(tcritnamepre,tcrit,'num',sep='')
	n2 <- paste(tcritnamepre,tcrit,'numerad',sep='')#name it
	n3 <- paste(tcritnamepre,tcrit,'perad',sep='')#name it
        names(parmat)[ii] <- n1
        names(parmat)[ii+1] <- n2
        names(parmat)[ii+2] <- n3
        ii <- ii + 3
    }
        
    TExtMat <- round2(TExtMat,2)
    nrowmat <- nrow(TExtMat)
    for(i in 1:nrowmat){
        tv <- parmat$tv[i]
###Which trials had pathogen until time VaccStartTime
        wiTrialsToVacc = which(TExtMat[i,] > round2(VaccStartTime + tv,2))

###Loop through tcritvals
        for(tcrit in tcritvals){
		  n1 <- paste(tcritnamepre,tcrit,'num',sep='')
		  n2 <- paste(tcritnamepre,tcrit,'numerad',sep='')#name it
		  n3 <- paste(tcritnamepre,tcrit,'perad',sep='')#name it

            if(length(wiTrialsToVacc)==0)
            {
                refparmat[i,n1] <- 0
		refparmat[i,n2] <- 0
		refparmat[i,n3] <- NA
            }else{
                NTrialsToVacc = length(wiTrialsToVacc)
		refparmat[i,n1] <- NTrialsToVacc					      
                refparmat[i,n2] <- sum(TExtMat[i,wiTrialsToVacc] <  
                                              round2(VaccStartTime + tv + tcrit,2))
		refparmat[i,n3] <- refparmat[i,n2]/refparmat[i,n1] 
            }
        }#End for loop
        if(i%%1000==0){print(paste('i = ', i, '/', nrowmat))}

    }###Loop through i   
    print('Writing refparmat')
    write.table(refparmat,file = paste(FileName,'/Ref_ParMat',sep=''), row.names = FALSE, col.names = TRUE)    
}}


wi <- with(parmat, R0p==5 & gamp==0.033)
plot(time365perad~tv,refparmat[wi,])
abline(lm(time365perad~tv,refparmat[wi,]))

Pr = with(refparmat,exp(-(d+gamp)*(365-tv)))
NPop = with(refparmat,NPeak/R0p*d/(d+gamp))

(1 - Pr)^(NPop)
