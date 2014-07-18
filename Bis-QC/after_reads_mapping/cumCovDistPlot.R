# cumCovDistPlot.R
# 
# Author: yaping
# Oct 30, 2013
# 10:51:04 PM
###############################################################################

Colors<-c("red","blue","orange","purple","darkred","black","yellow","cyan","brown","darkgreen","pink","gold","darkcyan","darkmagenta","navy","silver","aqua","grey")

for (e in commandArgs(TRUE)) {
	ta = strsplit(e,"=",fixed=TRUE)
	if(! is.na(ta[[1]][2])) {
		if(ta[[1]][1] == "wd"){
			wd<-ta[[1]][2]  ## directory to make plot file
			
		}
		if(ta[[1]][1] == "prefix"){
			prefix<-ta[[1]][2]
		}
		if(ta[[1]][1] == "input"){
			input<-ta[[1]][2]
		}
		if(ta[[1]][1] == "step"){
			step<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "xmin"){
			xmin<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "xmax"){
			xmax<-as.numeric(ta[[1]][2])
		}
	}
}
setwd(wd)
d <-read.table(input,sep="\t",header=T)
d<-d[order(d[,1]),]

axisSeqForPlot<-seq(xmin-1, xmax, by=step)
xmin<-xmin+1
xmax<-xmax+1
axisSeq<-seq(xmin, xmax, by=1)
mainTitle=prefix

pdf(paste(prefix,"pdf",sep="."), paper="special", height=4, width=4)
#par(oma=c(1, 1, 1, 1), mar=c(2, 2, 3, 0.1))
for(i in 1:dim(d)[1]){
	plot(axisSeq,1-d[i,axisSeq],type="l",xlim = c(xmin,xmax), ylim=c(0,1),xaxt="n",ann=F,yaxt="n",xlab="Percentile",ylab="Coverage",col=Colors[i],lty=1,font=2,lwd=2)
	par(new=T)
}
axis(1,at=seq(axisSeq[1],axisSeq[length(axisSeq)],length.out=length(axisSeqForPlot)),labels=axisSeqForPlot,lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.axis=2,font.lab=2,lwd=1)
axis(2,at=seq(0,1,by=0.2),labels=seq(0,100,by=20), lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.axis=2,font.lab=2,lwd=1, las=1)

title(xlab="Percentile", ylab="Coverage",cex.lab=1)
legend("bottomright",legend=d[,1], col=Colors[1:dim(d)[1]],lty=1,cex=0.7,lwd=2)
dev.off()
