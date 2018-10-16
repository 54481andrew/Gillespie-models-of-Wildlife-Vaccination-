#Create a lineplot of extinction data reading modified parmat in 
SimName = "B_lam"
FoldName = c("Data/")
FileName <- paste(FoldName, SimName,sep='')

parmat = read.table(file = paste(FileName,'/ParMat',sep=''), header = F)
names(parmat) = c('Par','b0','d','Bp','Nv','tv','gamv','gamp','tb',
	      'T','IpInit', 'TPathInv', 'lam', 'ExtRate')

if(grepl('Freq',SimName)){
	parmat$R0p = with(parmat, round( Bp/(d+gamp) ,2))
}else{
	parmat$R0p = with(parmat, round( Bp*(b0*tb)/(T*d*(d+gamp)) ,2))
}

parmat$R0approx = with(parmat, R0p*(1-Nv/(b0*tb + Nv*exp(-d*(T-tv)))))
parmat$rho = with(parmat, round(Nv/(b0*tb/(d*T)),2))

NPars = nrow(parmat)

#####################################
######## EXTINCTION PLOT C ##########
#####################################
##Contour plot. Y-axis: TPathInv
##X-axis: tv;
XValName = 'tv'
XVals = unique(parmat[,XValName])
YValName = 'R0p'
YVals = unique(parmat[,YValName])
ZValName = 'ExtRate'

FixValName1 = 'rho'
FixVals1 = c(0.5,1,1.5) #unique(parmat[,FixValName1])#c(0.00001,0.00005,0.0001)
FixValName2 = 'gamp'
FixVals2 = c(0.01,0.03, 0.07) 
FixValName3 = 'lam'
FixVals3 = c(0.005, 0.01) 
FixValName4 = 'd'
FixVals4 = c(0.00274)
F4 = 0.00274

FigFold = paste('Data/',SimName,'_Fig',sep='')
if(!dir.exists(FigFold)){dir.create(FigFold)}

nXVals <- length(unique(parmat[,XValName]))
nYVals <- length(YVals)#length(unique(parmat[,YValName]))
nFixVals1 <- length(FixVals1)
nFixVals2 <- length(FixVals2)
nFixVals3 <- length(FixVals3)

require(RColorBrewer)

i = 1
for(i1 in 1:nFixVals1){
    for(i2 in 1:nFixVals2){
    for(i3 in 1:nFixVals3){

        ExtMat = matrix(ncol = nYVals, nrow = nXVals)
        for(Yi in 1:nYVals){
                YVal = YVals[Yi]
                F1 = FixVals1[i1]
                F2 = FixVals2[i2]
                F3 = FixVals3[i3]
                
                wivals = which(parmat[,YValName]==YVal & 
		       parmat[,FixValName1]==F1 & 
		       parmat[,FixValName2]==F2 & 
		       parmat[,FixValName3]==F3 & 
		       parmat[,FixValName4]==F4)
                
###sort
                wivals = wivals[order(parmat[wivals,XValName])]
                parvals = parmat[wivals,ZValName]
                ExtMat[,Yi] = parvals
	    }#End loops through XVals and YVals

zmin = 0 #floor(10*min(ExtMat[!is.na(ExtMat)]))/10
zmax = 0.01 #1.5*max(ExtMat[!is.na(ExtMat)])
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
            FileName = paste('ExtLine_', FixValName1, F1, 
	    	     FixValName2, F2, FixValName3, F3, 
		     FixValName4,F4,'.png',sep='') 
            png(file = paste( FigFold, '/', FileName, sep=''), height = 4, 
	    	     width = 5, units = 'in', res = 400)
            par(mai = c(1,1,0.25,0.25))
	    ExtMat[is.na(ExtMat)] <- NA
	    matplot(x = XVals, y = ExtMat, col = cols, 
		    xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', 
		    ylim = c(0,zmax), type = 'l', lwd = 2, lty = 1)
            axislabs = seq(0,365, by = 60)
            axislabs1 = YVals
	    axis(side = 1, labels = T, at = axislabs)
	    axis(side = 2, labels = T)
	    mtext(text = axistext(XValName),side = 1, line = 3)
	    mtext(text = "Rate of Extinction",side = 2, line = 3)
	    abline(v = parmat$tb[wivals][1], lwd = 3, lty = 3)
	    legend(x = 'topright', legend = paste(axistext(YValName),'=',round(YVals,1)), col = cols, lwd = 2) 

dev.off()

print(paste("Figure", i, "complete",sep = " " ))
}else{
print(paste("Figure", i, "fail",sep = " " ))
print(head(ExtMat))
#print(parmat[wivals,])
}
i=i+1
    }}}#End Loop through fixed vals
