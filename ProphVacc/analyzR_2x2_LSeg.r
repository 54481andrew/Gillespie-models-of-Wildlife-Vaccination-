#Create a 1x2 lineplot of extinction data reading modified parmat in 
SimName = "B_1"
FoldName = c("Data/")
FileName <- paste(FoldName, SimName,sep='')

parmat = read.table(file = paste(FileName,'/Ref_ParMat',sep=''), header = T)
parmat$R0approx = with(parmat, R0p*(1-Nv/(b0*tb + Nv*exp(-d*(T-tv)))))
parmat$rho = with(parmat, round(Nv/(b0*tb/(d*T)),2))

NPars = nrow(parmat)

CloseToMax <- 0.05

#####################################
######## EXTINCTION PLOT C ##########
#####################################
##Contour plot. Y-axis: TPathInv
##X-axis: tv;
XValName = 'tv'
XVals = unique(parmat[,XValName])
YValName = 'R0p'
YVals = unique(parmat[,YValName])
ZValName = 'pc0.99'

#

FixValName1 = 'rho'
FixVals1 = c(0.5,1.5) #unique(parmat[,FixValName1])#c(0.00001,0.00005,0.0001)
FixValName2 = 'tb'
FixVals2 = c(30,90) 
FixValName3 = 'gamp'
FixVals3 = c(0.07, 0.01) 
FixValName4 = 'IpInit'
F4 = 5

TPVals <- unique(parmat$TPathInv)

FigFold = paste('Data/',SimName,'_Fig',sep='')
if(!dir.exists(FigFold)){dir.create(FigFold)}

nXVals <- length(unique(parmat[,XValName]))
nYVals <- length(YVals)#length(unique(parmat[,YValName]))
nFixVals1 <- length(FixVals1)
nFixVals2 <- length(FixVals2)
nFixVals3 <- length(FixVals3)

require(RColorBrewer)

colfun <- function(name, percent = 0.95){
       rgb.val <- col2rgb(name) 
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

            FileName = paste('LSeg_ExtLine_2x2_', FixValName1, F1, 
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
                
                wivalsind = parmat[,YValName]==YVal & 
		       parmat[,FixValName1]==F1 & 
		       parmat[,FixValName2]==F2 & 
		       parmat[,FixValName3]==F3 & 
		       parmat[,FixValName4]==F4
                wivals <- which(wivalsind)
###sort
		parvals <- 1:nXVals
		for(Xi in 1:nXVals){
		       wixvals <- wivalsind & parmat[,XValName]==XVals[Xi]
		       parvals[Xi] <- mean(parmat[wixvals,ZValName])
		}
                ExtMat[,Yi] = parvals
	    }#End loops through XVals and YVals

zmin = 0 #floor(10*min(ExtMat[!is.na(ExtMat)]))/10
zmax = 1 #1.5*max(ExtMat[!is.na(ExtMat)])
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
		    ylim = c(0,zmax), type = 'l', lwd = 0.5, lty = 1)

       xii = apply(ExtMat, MARGIN = 2, FUN = which.max)
for(ii in 1:ncol(ExtMat)){
       XMax = XVals[xii[ii]]
       MaxVal = ExtMat[xii[ii],ii]
       CompValsInd = which( abs(ExtMat[, ii] - MaxVal) < CloseToMax )
       CompValsX = XVals[CompValsInd]
       CompValsY = ExtMat[CompValsInd,ii]
##       segments(x0 = min(CompValsX), y0 = MaxVal, x1 = max(CompValsX), y1 = MaxVal, col = cols[ii])

	 polygon(x = c(min(CompValsX), min(CompValsX), max(CompValsX), max(CompValsX)), 
	 y = c(min(CompValsY),max(CompValsY), max(CompValsY), min(CompValsY)), col = 
	 paste(cols[ii],"75",sep=''))

##	 polygon(x = c(min(CompValsX), min(CompValsX), max(CompValsX), max(CompValsX)), 
##	 y = CloseToMax*c(-1, 0, 0, -1) + MaxVal, border = cols[ii], lwd = 3)

       matpoints(XMax, MaxVal, pch = 8, lwd = 1, col = "black")##, col = cols[ii])
}


tbval <- unique(parmat$tb[wivals])
polygon(x = c(0,tbval, tbval, 0), y = c(-1,-1,1,1), col = colfun('gray', 75), border =NA)

            axislabs = seq(0,365, by = 30)
            axislabs1 = YVals

	    axis(side = 1, labels = if(i2==2){c(0,NA,NA,NA,120,NA,NA,NA,240, NA,NA,NA,360)}else{F}, at = axislabs, padj = -0.8)
	    axis(side = 2, labels = i3==1, padj = 0.6)

	    if(i2==2){mtext(text = axistext(XValName),side = 1, line = 2)}
	    if(i3==1){mtext(text = "Probability of Extinction", side = 2, line = 2.25)}

	    abline(v = parmat$tb[wivals][1], lwd = 3, lty = 3)
	    if(i2==1 & i3==1){legend(x = 'topright', legend = 
	    paste(axistext(YValName),'=',round(YVals,1)), 
	    col = cols, lwd = 4, cex = 0.65) }

#Variable name varied across panels
##legendtext = c(paste(round(1/FixVals3[i3]),'day infection'), paste(FixVals2[i2], 'vaccines/host'))
legendtext = c(paste(round(1/FixVals3[i3]),'day infection'), paste(FixVals2[i2], 'breeding season'))
legend(x = -60, y = 1.075, legend = legendtext, fill = 'white', box.col = 'black', cex = 1)

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
