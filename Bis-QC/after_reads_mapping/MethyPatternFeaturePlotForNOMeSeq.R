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
		if(ta[[1]][1] == "hcgfn"){
			hcgfn<-ta[[1]][2]
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

	}
}
setwd(wd)
gch<-read.table(gchfn,sep="\t",header=F)
hcg<-read.table(hcgfn,sep="\t",header=F)

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
valueHcg<-array()
for(i in dataSeq){
	valueGch<-cbind(valueGch,mean(colMeans(gch[,(i-step/2):(i+step/2-1)], na.rm=T), na.rm=T))
	valueHcg<-cbind(valueHcg,mean(colMeans(hcg[,(i-step/2):(i+step/2-1)], na.rm=T), na.rm=T))
	
}
valueGch<-valueGch[,2:length(valueGch[1,])]
valueHcg<-valueHcg[,2:length(valueHcg[1,])]
numElemGch<-length(gch[,scale+5])
#numHaveValueElemGch<-length(gch[!is.na(mean(gch[,((length(gch[1,])-4)/2):((length(gch[1,])-4)/2+step-1)], na.rm=T)),][,(length(gch[1,])-4)/2])
numElemHcg<-length(hcg[,scale+5])
#numHaveValueElemHcg<-length(hcg[!is.na(mean(hcg[,dataSeq[length(dataSeq)/2]:dataSeq[length(dataSeq)/2]+step-1], na.rm=T)),][,dataSeq[length(dataSeq)/2]])
mainTitle=paste(prefix,"\nsmooth:", smooth," ---numElem:",numElemGch,",","numElemHcg:",numElemHcg)

#pdf(paste(gchfn,"pdf",sep="."),height=24, width=24, pointsize=12)
#postscript(paste(gchfn,"eps",sep="."), height=15, width=8)
pdf(paste(gchfn,"pdf",sep="."), paper="special", height=4, width=4)
par(oma=c(1, 1, 1, 1))
par(mar=c(1, 1, 3, 1))
#par(oma = c(4, 0, 0, 0))

plot(axisSeq,valueGch,type="l",axes=FALSE,ylim=c(y_scale_min,y_scale_max),xlab="",ylab="",col=c("#00CC99"),lty=1,font=2,lwd=3)
par(new=T)
#par(new=T)
plot(axisSeq,valueHcg,type="l",axes=FALSE,xlab="",ylab="",ylim=c(y_scale_min,y_scale_max),col=c("black"),lty=1,font=2,lwd=3)
#par(new=T)
par(new=T)
#par(mgp = c(0, 1, 0))

axis(1,at=axisSeqForPlot,lty=1,font=2,cex.axis=1.2,cex.lab=1.2,font.lab=2,lwd=2)
#text(axisSeqForPlot,rep(-0.06,9),labels=axisSeqForPlot,font=2,cex=4)
#par(mgp = c(0, 1, 0))
axis(2,at=seq(y_scale_min,y_scale_max,by=y_step),lty=1,font=2,cex.axis=1.2,cex.lab=1.2,font.lab=2,lwd=2)
title(mainTitle, cex.main = 0.6, font.main= 4, col.main= "black",xlab="Distance to elements (bp)", ylab="Methylation")
#sub = "subTitle",cex.sub = 0.75, font.sub = 3, col.sub = "red")
legend("topright",c("GCH","HCG"),col=c("#00CC99","black"),lty=1,cex=1.2,lwd=2)
abline(v=0)
dev.off()

