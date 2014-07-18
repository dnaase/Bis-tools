# methyBiasDistPlot.R
# 
# Author: yaping
# Jul 15, 2014
# 2:05:44 PM
###############################################################################


for (e in commandArgs(TRUE)) {
	ta = strsplit(e,"=",fixed=TRUE)
	if(! is.na(ta[[1]][2])) {

		if(ta[[1]][1] == "input"){
			input<-ta[[1]][2]
		}
		if(ta[[1]][1] == "output"){
			output<-ta[[1]][2]
		}
	}
}
d <-read.table(input,sep="\t",header=F)
xmin<-d[1,1]
xmax<-d[length(d[,1]),1]

axisSeq<-seq(xmin, xmax, by=1)

pdf(output, paper="special", height=4, width=4)
#par(oma=c(0, 0, 0, 0), mar=c(3, 3, 3, 0.1))

plot(d[,1],d[,4],type="b",xlim = c(xmin,xmax), ylim=c(0,100),xaxt="n",ann=F,yaxt="n",xlab="Cycle(bp)",ylab="Methylation",lty=1,font=2,lwd=1,cex=0.8,pch=16)

axis(1,at=seq(xmin,xmax,length.out=3),labels=c(xmin,0,xmax),lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.axis=2,font.lab=2,lwd=1)
axis(2,at=seq(0,100,by=20),labels=seq(0,100,by=20), lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.axis=2,font.lab=2,lwd=1, las=1)

title(main="M-Bias",xlab="Cycle(bp)", ylab="Methylation",cex.lab=1)
dev.off()
