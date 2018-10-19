#Create a 1x2 lineplot of extinction data reading modified parmat in 
SimName = "C_lam"
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

#

FixValName1 = 'd'
FixVals1 = c(0.00274,0.00548) #unique(parmat[,FixValName1])#c(0.00001,0.00005,0.0001)
FixValName2 = 'tb'
FixVals2 = c(60,90) 
FixValName3 = 'gamp'
FixVals3 = c(0.07, 0.01) 
FixValName4 = 'rho'
F4 = 1

FigFold = paste('Data/',SimName,'_Fig',sep='')
if(!dir.exists(FigFold)){dir.create(FigFold)}

nXVals <- length(unique(parmat[,XValName]))
nYVals <- length(YVals)#length(unique(parmat[,YValName]))
nFixVals1 <- length(FixVals1)
nFixVals2 <- length(FixVals2)
nFixVals3 <- length(FixVals3)

require(RColorBrewer)

colfun <- function(name, percent){
       rgb.val <- col2rgb("gray") 
       ## Make new color using input color as base and alpha set by transparency
       t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
	invisible(t.col)
}       
	       

i = 1
for(i1 in 1:nFixVals1){

                F1 = FixVals1[i1]

##################
###Begin figure###
##################

            FileName = paste('ExtLine_2x2_', FixValName1, F1, 
		     FixValName4,F4,'.png',sep='') 
            png(file = paste( FigFold, '/', FileName, sep=''), height = 5.5, 
	    	     width = 5.5, units = 'in', res = 400)
            par(mai = c(0.0,0.1,0.1,0.1), omi = c(0.65, 0.6, 0, 0), mfrow = c(2,2))


    for(i2 in 1:nFixVals2){
    for(i3 in 1:nFixVals3){
                
		F2 = FixVals2[i2]

        ExtMat = matrix(ncol = nYVals, nrow = nXVals)
        for(Yi in 1:nYVals){
                YVal = YVals[Yi]
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
	    ExtMat[is.na(ExtMat)] <- NA
	    matplot(x = XVals, y = ExtMat, col = cols, 
		    xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', 
		    ylim = c(0,zmax), type = 'l', lwd = 2, lty = 1)


tbval <- unique(parmat$tb[wivals])
polygon(x = c(0,tbval, tbval, 0), y = c(-1,-1,1,1), col = colfun('gray', 75), border =NA)

            axislabs = seq(0,365, by = 120)
            axislabs1 = YVals

	    axis(side = 1, labels = i2==2, at = axislabs, padj = -0.8)
	    axis(side = 2, labels = i3==1, padj = 0.6)

	    if(i2==2){mtext(text = axistext(XValName),side = 1, line = 2)}
	    if(i3==1){mtext(text = "Rate of Extinction", side = 2, line = 2.25)}

	    abline(v = parmat$tb[wivals][1], lwd = 3, lty = 3)
	    if(i2==1 & i3==1){legend(x = 'bottomright', legend = paste(axistext(YValName),'=',round(YVals,1)), col = cols, lwd = 2) }

#Variable name varied across panels
#legendtext = c(paste(round(1/FixVals3[i3]),'day infection'), paste(FixVals2[i2], 'vaccines / host'))
legendtext = c(paste(round(1/FixVals3[i3]),'day infection'), paste(FixVals2[i2], 'breeding season'))
legend(x = 120, y = 0.01075, legend = legendtext, bty = 'n')


print(paste("Pane", i, "complete",sep = " " ))
}else{
print(paste("Pane", i, "fail",sep = " " ))
print(head(ExtMat))
#print(parmat[wivals,])
}
i=i+1
    }}##End loop through figure panes

dev.off()

}#End Loop through fixed vals
