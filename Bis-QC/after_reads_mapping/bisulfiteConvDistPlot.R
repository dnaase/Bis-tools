# bisulfiteConvDistPlot.R
# 
# Author: yaping
# Jul 15, 2014
# 3:08:19 PM
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
data<-100*d[,1]/d[,2]

pdf(output, paper="special", height=4, width=4)
#par(oma=c(0, 0, 0, 0), mar=c(3, 3, 3, 0.1))

hist(data,breaks=seq(0,100,by=1),freq=F,xlim = c(0,100),xlab="Methylation",ylab="Frequency",main="Ditribution of bisulfite conversion\n rate within read")

dev.off()

