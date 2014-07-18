# TODO: Add comment
# 
# Author: yaping
###############################################################################
#parameter: work directory, file_name_prefix, input gch, hcg/wcg files, step pace to draw line, scale, smooth or not


for (e in commandArgs(TRUE)) {
	ta = strsplit(e,"=",fixed=TRUE)
	if(! is.na(ta[[1]][2])) {
		if(ta[[1]][1] == "wd"){
			wd<-ta[[1]][2]  ## directory to make plot file
			
		}
		if(ta[[1]][1] == "prefix"){
			prefix<-ta[[1]][2]
		}
		if(ta[[1]][1] == "gchfn"){
			gchfn<-ta[[1]][2]
		}
		if(ta[[1]][1] == "step"){
			step<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "scale"){
			scale<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "axistep"){
			axistep<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "smooth"){
			smooth<-as.numeric(ta[[1]][2])  ## 0: no smooth; 1: smooth step;
		}
		if(ta[[1]][1] == "y_scale_min"){
			y_scale_min<-as.numeric(ta[[1]][2])  
		}
		if(ta[[1]][1] == "y_scale_max"){
			y_scale_max<-as.numeric(ta[[1]][2])  
		}
		if(ta[[1]][1] == "y_step"){
			y_step<-as.numeric(ta[[1]][2])  
		}
		if(ta[[1]][1] == "color"){
			line_color<-ta[[1]][2]
		}
#		if(ta[[1]][1] == "sumFeatureFn"){
#			sumFeatureFn<-ta[[1]][2]
#		}
#		if(ta[[1]][1] == "numFeature"){
#			numFeature<-ta[[1]][2]
#		}
#		if(ta[[1]][1] == "legendName"){
#			legendName<-ta[[1]][2]
#		}
	}
}
setwd(wd)
gch<-read.table(gchfn,sep="\t",header=F)
#sumFeature<-read.table(sumFeatureFn,sep="\t",header=F)
#sumFeature<-sum(sumFeature,na.rm=T)

axisSeqForPlot<-seq(0-scale, scale, by=axistep)
if(smooth == 0){ ## no smooth
	axisSeq<-seq(0-scale, scale, by=step)
	dataSeq<-seq(5+((length(gch[1,])-5)/2)-scale, 5+((length(gch[1,])-5)/2)-scale+scale * 2, by=step)
}
if(smooth != 0){ ## smooth
	axisSeq<-seq(0-scale, scale, by=smooth)
	dataSeq<-seq(5+((length(gch[1,])-5)/2)-scale, 5+((length(gch[1,])-5)/2)-scale+scale * 2, by=smooth)
	
}
valueGch<-array()

for(i in dataSeq){
	valueGch<-cbind(valueGch,mean(colMeans(gch[,(i-step/2):(i+step/2-1)], na.rm=T), na.rm=T))
	
	
}
valueGch<-valueGch[,2:length(valueGch[1,])]

numElemGch<-length(gch[,scale+5])
#numHaveValueElemGch<-length(gch[!is.na(mean(gch[,((length(gch[1,])-4)/2):((length(gch[1,])-4)/2+step-1)], na.rm=T)),][,(length(gch[1,])-4)/2])

#numHaveValueElemHcg<-length(hcg[!is.na(mean(hcg[,dataSeq[length(dataSeq)/2]:dataSeq[length(dataSeq)/2]+step-1], na.rm=T)),][,dataSeq[length(dataSeq)/2]])
mainTitle=paste(prefix,"smooth:", smooth," \n-numElemCenterToAlign:",numElemGch)


#pdf(paste(gchfn,"pdf",sep="."),height=24, width=24, pointsize=12)
#postscript(paste(gchfn,"eps",sep="."), height=15, width=8)
#par(oma = c(4, 0, 0, 0))
pdf(paste(gchfn,"pdf",sep="."), paper="special", height=4, width=4)
par(oma=c(1, 1, 1, 1))
par(mar=c(1, 1, 3, 1))
#plot(axisSeq,valueGch*10,type="l",axes=FALSE,ylim=c(0.00,1.0),xlab="",ylab="",col=c("black"),lty=1,font=2,lwd=8)
plot(axisSeq,valueGch,type="l",axes=FALSE,xlab="",ylab="",ylim=c(y_scale_min,y_scale_max),col=line_color,lty=1,font=2,lwd=3)

axis(1,at=axisSeqForPlot,lty=1,font=2,cex.axis=1.2,cex.lab=1.2,font.lab=2,lwd=2)

axis(2,at=seq(y_scale_min,y_scale_max,by=y_step),lty=1,font=2,cex.axis=1.2,cex.lab=1.2,font.lab=2,lwd=2)
title(mainTitle, cex.main = 0.6, font.main= 4, col.main= "black",xlab="Distance to elements (bp)", ylab="Methylation level")
legend("topright","Methylation level", col=line_color,lty=1,cex=1.2,lwd=2)
abline(v=0)
dev.off()

