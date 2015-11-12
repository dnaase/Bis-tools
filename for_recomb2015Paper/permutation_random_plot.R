library(fields)
fudgeit <- function(){
	xm <- get('xm', envir = parent.frame(1))
	ym <- get('ym', envir = parent.frame(1))
	z  <- get('dens', envir = parent.frame(1))
	colramp <- get('colramp', parent.frame(1))
	image.plot(xm,ym,z, col = colramp(256), legend.only = T, add =F)
}
Lab.palette <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))


random_hics<-NULL
for (e in commandArgs(TRUE)) {
	ta = strsplit(e,"=",fixed=TRUE)
	if(! is.na(ta[[1]][2])) {
		if(ta[[1]][1] == "wd"){
			wd<-ta[[1]][2]  ## directory to make plot file
		}
		if(ta[[1]][1] == "prefix"){
			prefix<-ta[[1]][2]
		}
		if(ta[[1]][1] == "meqtl_hic"){
			meqtl_hic<-ta[[1]][2]
		}
		if(ta[[1]][1] == "random_hics"){
			random_hics<-c(random_hics, ta[[1]][2])
		}
		
	}
}

setwd(wd)
meqtl<-read.table(meqtl_hic,sep="\t",header=F)
rownames(meqtl)<-meqtl[,4]

random_array<-read.table(random_hics,sep="\t",header=F)
rownames(random_array)<-random_array[,4]
#for(random_hic in random_hics){
#	a<-read.table(random_hic,sep="\t",header=F)
#	rownames(a)<-a[,4]
#	if(is.null(random_array)){
#		random_array<-as.matrix(a[,length(a[1,])])
#		rownames(random_array)<-a[,4]
		
#	}else{
#		common<-intersect(rownames(random_array),a[,4])
#		b<-a[common,]
#		random_array<-cbind(random_array[common,],b[,length(b[1,])])
#	}

#}

common<-intersect(rownames(meqtl),rownames(random_array))
#x<-cbind(meqtl[common,c(1:3,length(meqtl[1,]))],rowMeans(random_array[common,],na.rm=T))
x<-cbind(meqtl[common,c(1:3,length(meqtl[1,]))],random_array[common,length(random_array[1,])])
rm(meqtl)
rm(random_array)
rm(common)

pdf(paste("recombination_rate",prefix,"pdf",sep="."), paper="special", height=4, width=4)
par(plt = c(0.1171429 ,0.8200000, 0.1457143, 0.8828571))
par(mar = c(5,4,4,5) + .1)
plot(0,0,xlim=c(0,2.0),ylim=c(0,2.0),xlab="",ylab="")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="#00007F")
par(new=T,mar = c(5,4,4,5) + .1)
print(dim(x))
print(summary(log10(x[,4:5]+1.1)))
smoothScatter(cbind(log10(x[,4]+1.1),log10(x[,5]+1.1)),xlim=c(0,2.0),ylim=c(0,2.0),ylab="log10(contact_freq+1),Random_pair",xlab="log10(contact_freq+1),pair",main=prefix,colramp = Lab.palette,postPlotHook = fudgeit)
abline(a=0,b=1)
#abline(lm(log10(x[,5]+1.1) ~ log10(x[,4]+1.1)),col="purple")
#x<-x[rowSums(is.na(x))==0,]
#lines(lowess(log10(x[,4]+1.1),log10(x[,5]+1.1)), col="purple")
dev.off()

####make boxplot and significant p value around different scale

smallRange<-x[x[,3]-x[,2]<=10000,]
midRange<-x[x[,3]-x[,2]>10000 & x[,3]-x[,2]<=100000,]
largeRange<-x[x[,3]-x[,2]>100000,]

main_tilte=""
if(dim(smallRange)[1]>1){
	p.smallRange<-t.test(smallRange[,4],smallRange[,5], paired=T, alternative="less")$p.value
	print(p.smallRange)
	main_tilte<-paste(main_tilte,format(p.smallRange,digits=3),sep="")
}
if(dim(midRange)[1]>1){
	p.midRange<-t.test(midRange[,4],midRange[,5], paired=T, alternative="less")$p.value
	print(p.midRange)
	main_tilte<-paste(main_tilte,format(p.midRange,digits=3),sep=",")
}

if(dim(largeRange)[1]>1){
	p.largeRange<-t.test(largeRange[,4],largeRange[,5], paired=T, alternative="less")$p.value
	print(p.largeRange)
	main_tilte<-paste(main_tilte,format(p.largeRange,digits=3),sep=",")
}

if(dim(smallRange)[1]>1){
	p.smallRange<-wilcox.test(smallRange[,4],smallRange[,5], paired=T, alternative="less")$p.value
	print(p.smallRange)
	main_tilte<-paste(main_tilte,format(p.smallRange,digits=3),sep="\n")
}
if(dim(midRange)[1]>1){
	p.midRange<-wilcox.test(midRange[,4],midRange[,5], paired=T, alternative="less")$p.value
	print(p.midRange)
	main_tilte<-paste(main_tilte,format(p.midRange,digits=3),sep=",")
}

if(dim(largeRange)[1]>1){
	p.largeRange<-wilcox.test(largeRange[,4],largeRange[,5], paired=T, alternative="less")$p.value
	print(p.largeRange)
	main_tilte<-paste(main_tilte,format(p.largeRange,digits=3),sep=",")
}




s<-log10(smallRange[,4]+1)
s_r<-log10(smallRange[,5]+1)
m<-log10(midRange[,4]+1)
m_r<-log10(midRange[,5]+1)
l<-log10(largeRange[,4]+1)
l_r<-log10(largeRange[,5]+1)

z <- c("s", "s_r", "m", "m_r","l", "l_r")
dataList <- lapply(z, get, envir=environment())
names(dataList) <- z

pdf(paste("recombRate_boxplot",prefix,"pdf",sep="."), paper="special", height=4, width=4)
boxplot(dataList,ylim=c(0,1.0), outline=F,col=c("red","grey"),main=main_tilte)
dev.off()


meqtl_col<-rgb(col2rgb("red")[1]/255,col2rgb("red")[2]/255,col2rgb("red")[3]/255,alpha=0.3)
random_col<-rgb(col2rgb("black")[1]/255,col2rgb("black")[2]/255,col2rgb("black")[3]/255,alpha=0.3)

##calculate average Hi-C contact frequency along distance (smoothed version)
distance<-seq(0,1000000,by=1000)
prev=0
meQTL_avg<-NULL
random_avg<-NULL
meQTL_sd<-NULL
random_sd<-NULL

for(i in distance){
	start<-i
	end<-i+10000
	meQTL_avg<-c(meQTL_avg,mean(x[(x[,3]-x[,2])<end & (x[,3]-x[,2])>=start,][,4],na.rm=T))
	random_avg<-c(random_avg,mean(x[(x[,3]-x[,2])<end & (x[,3]-x[,2])>=start,][,5],na.rm=T))
	meQTL_sd<-c(meQTL_sd,sd(x[(x[,3]-x[,2])<end & (x[,3]-x[,2])>=start,][,4],na.rm=T))
	random_sd<-c(random_sd,sd(x[(x[,3]-x[,2])<end & (x[,3]-x[,2])>=start,][,5],na.rm=T))
	prev=i
}
	
meqtl.low<-meQTL_avg-1.96*meQTL_sd/10
meqtl.high<-meQTL_avg+1.96*meQTL_sd/10
random.low<-random_avg-1.96*random_sd/10
random.high<-random_avg+1.96*random_sd/10

plot_range<-c(2:(length(distance)-2))
na_collection<-c(which(is.na(meqtl.high)),which(is.na(meqtl.low)),which(is.nan(meqtl.high)),which(is.nan(meqtl.low)))
if(length(na_collection)>0){
	na_start_index<-min(na_collection,na.rm=T)-1
	plot_range<-c(2:na_start_index)
}

plot_range_random<-c(2:(length(distance)-2))
na_collection_random<-c(which(is.na(random.high)),which(is.na(random.low)),which(is.nan(random.high)),which(is.nan(random.low)))
if(length(na_collection_random)>0){
	na_start_index<-min(na_collection_random,na.rm=T)-1
	plot_range_random<-c(2:na_start_index)
}


pdf(paste("recombination_rate",prefix,"freq_vs_Log10distance.ave_plot.pdf",sep="."), paper="special", height=4, width=4)
plot(log10(distance),log10(meQTL_avg+1.1),xlab="",ylab="",main="",col="red",xlim=c(3,6),ylim=c(0,0.6),axes=FALSE,type="l")
par(new=T)
plot(log10(distance),log10(random_avg+1.1),axes=FALSE,xlab="",ylab="",main=,col="black",xlim=c(3,6),ylim=c(0,0.6),cex.main=0.7,type="l")
#lines(log10(distance), log10(meqtl.low+1.1), col = meqtl_col)
#lines(log10(distance), log10(meqtl.high+1.1), col = meqtl_col)
polygon(c(log10(distance[plot_range]), rev(log10(distance[plot_range]))), c(log10(meqtl.high[plot_range]+1.1), rev(log10(meqtl.low[plot_range]+1.1))), col = meqtl_col, border = NA)
#lines(log10(distance), random.low, col = random_col)
#lines(log10(distance), random.high, col = random_col)
polygon(c(log10(distance[plot_range_random]), rev(log10(distance[plot_range_random]))), c(log10(random.high[plot_range_random]+1.1), rev(log10(random.low[plot_range_random]+1.1))), col = random_col, border = NA)

axis(2,at=seq(0,0.6,by=0.3),lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.lab=2,lwd=2)
axis(1,at=seq(3,6,by=1),labels=c("1kb","10kb","100kb","1Mb"),lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.lab=2,lwd=2)
legend("topright",c("pair","random_pair"), col=c("red","black"),lty=1,cex=0.8,lwd=2)
title(paste("recombination_rate.freq_vs_Log10distance",prefix, sep="\n"), cex.main = 0.6, font.main= 1, col.main= "black",xlab="log10(Distance)", ylab="log10(Recombination rate+1)")
dev.off()

pdf(paste("recombination_rate",prefix,"freq_vs_Log10distance.ave_plot_noLog10.pdf",sep="."), paper="special", height=4, width=4)
plot(log10(distance),meQTL_avg,xlab="",ylab="",main="",col="red",xlim=c(3,6),ylim=c(0,3),axes=FALSE,type="l")
par(new=T)
plot(log10(distance),random_avg,axes=FALSE,xlab="",ylab="",main=,col="black",xlim=c(3,6),ylim=c(0,3),cex.main=0.7,type="l")
#lines(log10(distance), meqtl.low, col = meqtl_col)
#lines(log10(distance), meqtl.high, col = meqtl_col)
polygon(c(log10(distance[plot_range]), rev(log10(distance[plot_range]))), c(meqtl.high[plot_range], rev(meqtl.low[plot_range])), col = meqtl_col, border = NA)
#lines(log10(distance), random.low, col = random_col)
#lines(log10(distance), random.high, col = random_col)
polygon(c(log10(distance[plot_range_random]), rev(log10(distance[plot_range_random]))), c(random.high[plot_range_random], rev(random.low[plot_range_random])), col = random_col, border = NA)

axis(2,at=seq(0,3,by=1),lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.lab=2,lwd=2)
axis(1,at=seq(3,6,by=1),labels=c("1kb","10kb","100kb","1Mb"),lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.lab=2,lwd=2)
legend("topright",c("pair","random_pair"), col=c("red","black"),lty=1,cex=0.8,lwd=2)
title(paste("recombination_rate.freq_vs_Log10distance",prefix, sep="\n"), cex.main = 0.6, font.main= 1, col.main= "black",xlab="log10(Distance)", ylab="Recombination rate")
dev.off()


##calculate average Hi-C contact frequency along distance (not smoothed version)
distance<-seq(0,1000000,by=1000)
prev=0
meQTL_avg<-NULL
random_avg<-NULL
meQTL_sd<-NULL
random_sd<-NULL

for(i in distance){
	start<-i
	end<-i+1000
	meQTL_avg<-c(meQTL_avg,mean(x[(x[,3]-x[,2])<end & (x[,3]-x[,2])>=start,][,4],na.rm=T))
	random_avg<-c(random_avg,mean(x[(x[,3]-x[,2])<end & (x[,3]-x[,2])>=start,][,5],na.rm=T))
	meQTL_sd<-c(meQTL_sd,sd(x[(x[,3]-x[,2])<end & (x[,3]-x[,2])>=start,][,4],na.rm=T))
	random_sd<-c(random_sd,sd(x[(x[,3]-x[,2])<end & (x[,3]-x[,2])>=start,][,5],na.rm=T))

	prev=i
}

meqtl.low<-meQTL_avg-1.96*meQTL_sd/10
meqtl.high<-meQTL_avg+1.96*meQTL_sd/10
random.low<-random_avg-1.96*random_sd/10
random.high<-random_avg+1.96*random_sd/10

plot_range<-c(2:(length(distance)-2))
na_collection<-c(which(is.na(meqtl.high)),which(is.na(meqtl.low)),which(is.nan(meqtl.high)),which(is.nan(meqtl.low)))
if(length(na_collection)>0){
	na_start_index<-min(na_collection,na.rm=T)-1
	plot_range<-c(2:na_start_index)
}

plot_range_random<-c(2:(length(distance)-2))
na_collection_random<-c(which(is.na(random.high)),which(is.na(random.low)),which(is.nan(random.high)),which(is.nan(random.low)))
if(length(na_collection_random)>0){
	na_start_index<-min(na_collection_random,na.rm=T)-1
	plot_range_random<-c(2:na_start_index)
}


pdf(paste("recombination_rate",prefix,"freq_vs_Log10distance.ave_plot.notSmoothed.pdf",sep="."), paper="special", height=4, width=4)
plot(log10(distance),log10(meQTL_avg+1.1),xlab="",ylab="",main="",col="red",xlim=c(3,6),ylim=c(0,0.6),axes=FALSE,type="l")
par(new=T)
plot(log10(distance),log10(random_avg+1.1),axes=FALSE,xlab="",ylab="",main=,col="black",xlim=c(3,6),ylim=c(0,0.6),cex.main=0.7,type="l")
#lines(log10(distance), log10(meqtl.low+1.1), col = meqtl_col)
#lines(log10(distance), log10(meqtl.high+1.1), col = meqtl_col)
polygon(c(log10(distance[plot_range]), rev(log10(distance[plot_range]))), c(log10(meqtl.high[plot_range]+1.1), rev(log10(meqtl.low[plot_range]+1.1))), col = meqtl_col, border = NA)
#lines(log10(distance), random.low, col = random_col)
#lines(log10(distance), random.high, col = random_col)
polygon(c(log10(distance[plot_range_random]), rev(log10(distance[plot_range_random]))), c(log10(random.high[plot_range_random]+1.1), rev(log10(random.low[plot_range_random]+1.1))), col = random_col, border = NA)

axis(2,at=seq(0,0.6,by=0.3),lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.lab=2,lwd=2)
axis(1,at=seq(3,6,by=1),labels=c("1kb","10kb","100kb","1Mb"),lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.lab=2,lwd=2)
legend("topright",c("pair","random_pair"), col=c("red","black"),lty=1,cex=0.8,lwd=2)
title(paste("recombination_rate.freq_vs_Log10distance",prefix, sep="\n"), cex.main = 0.6, font.main= 1, col.main= "black",xlab="log10(Distance)", ylab="log10(Recombination rate+1)")
dev.off()

pdf(paste("recombination_rate",prefix,"freq_vs_Log10distance.ave_plot_noLog10.notSmoothed.pdf",sep="."), paper="special", height=4, width=4)
plot(log10(distance),meQTL_avg,xlab="",ylab="",main="",col="red",xlim=c(3,6),ylim=c(0,3),axes=FALSE,type="l")
par(new=T)
plot(log10(distance),random_avg,axes=FALSE,xlab="",ylab="",main=,col="black",xlim=c(3,6),ylim=c(0,3),cex.main=0.7,type="l")
#lines(log10(distance), meqtl.low, col = meqtl_col)
#lines(log10(distance), meqtl.high, col = meqtl_col)
polygon(c(log10(distance[plot_range]), rev(log10(distance[plot_range]))), c(meqtl.high[plot_range], rev(meqtl.low[plot_range])), col = meqtl_col, border = NA)
#lines(log10(distance), random.low, col = random_col)
#lines(log10(distance), random.high, col = random_col)
polygon(c(log10(distance[plot_range_random]), rev(log10(distance[plot_range_random]))), c(random.high[plot_range_random], rev(random.low[plot_range_random])), col = random_col, border = NA)

axis(2,at=seq(0,3,by=1),lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.lab=2,lwd=2)
axis(1,at=seq(3,6,by=1),labels=c("1kb","10kb","100kb","1Mb"),lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.lab=2,lwd=2)
legend("topright",c("pair","random_pair"), col=c("red","black"),lty=1,cex=0.8,lwd=2)
title(paste("recombination_rate.freq_vs_Log10distance",prefix, sep="\n"), cex.main = 0.6, font.main= 1, col.main= "black",xlab="log10(Distance)", ylab="Recombination rate")
dev.off()



