# TODO: Add comment
# 
# Author: yaping
###############################################################################

fns<-NULL
colors<-NULL
legendNames<-NULL
line_types<-NULL
bin_size_align=1
y_scale_min<-NULL
y_scale_max<-NULL
lowCov=FALSE
for (e in commandArgs(TRUE)) {
	ta = strsplit(e,"=",fixed=TRUE)
	if(! is.na(ta[[1]][2])) {
		if(ta[[1]][1] == "wd"){
			wd<-ta[[1]][2]  ## directory to make plot file
			
		}
		if(ta[[1]][1] == "prefix"){
			prefix<-ta[[1]][2]
		}
		if(ta[[1]][1] == "fn"){
			fns<-c(fns, ta[[1]][2])
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
		if(ta[[1]][1] == "bin_size_align"){
			bin_size_align<-as.numeric(ta[[1]][2])
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
			colors<-c(colors,ta[[1]][2])
		}		
		if(ta[[1]][1] == "legendName"){
			legendNames<-c(legendNames,ta[[1]][2])
		}
		if(ta[[1]][1] == "autoScale"){
			autoScale<-ta[[1]][2]
		}
		if(ta[[1]][1] == "line_types"){
			line_types<-c(line_types,as.numeric(ta[[1]][2]))
		}
		if(ta[[1]][1] == "lowCov"){
			lowCov<-as.logical(ta[[1]][2])
		}
	}
}
if(is.null(line_types)){
	line_types<-rep(length(fns),1)
}
setwd(wd)

axisSeqForPlot<-seq(0-scale, scale, by=axistep)

pdf(paste(fns[1],"pdf",sep="."), paper="special", height=4, width=4)
par(oma=c(1, 1, 1, 1))
par(mar=c(1, 1, 3, 1))

data_max=NA
data_min=NA


if(autoScale){
	for(fn in fns){
		gch<-read.table(fn,sep="\t",header=F)	
		#gch[gch[,4]%in%"-",5:length(gch[1,])]=0-gch[gch[,4]%in%"-",5:length(gch[1,])]
		if(smooth == 0){ ## no smooth
			axisSeq<-seq(0-scale, scale, by=step)
			dataSeq<-seq(5+((length(gch[1,])-5)/2)-as.integer(scale/bin_size_align), 5+((length(gch[1,])-5)/2)+as.integer(scale/bin_size_align), by=as.integer(step/bin_size_align))
		}
		if(smooth != 0){ ## smooth
			axisSeq<-seq(0-scale, scale, by=smooth)
			dataSeq<-seq(5+((length(gch[1,])-5)/2)-as.integer(scale/bin_size_align), 5+((length(gch[1,])-5)/2)+as.integer(scale/bin_size_align), by=as.integer(smooth/bin_size_align))
			
		}
		
		valueGch<-NULL	
		
		
		for(i in dataSeq){
			if(floor(i+step/(2*bin_size_align)-1)-ceiling(i-step/(2*bin_size_align)) <= 0){
				valueGch<-cbind(valueGch,mean(gch[,ceiling(i-step/(2*bin_size_align))], na.rm=T))	
			}else{
				valueGch<-cbind(valueGch,mean(colMeans(gch[,ceiling(i-step/(2*bin_size_align)):floor(i+step/(2*bin_size_align)-1)], na.rm=T), na.rm=T))	
			}
			
		}
		if(is.null(y_scale_max)){
			data_max<-max(data_max,max(valueGch,na.rm =T),na.rm =T)
		}else{
			data_max<-y_scale_max
		}
		
		if(is.null(y_scale_min)){
			data_min<-min(data_min,min(valueGch,na.rm =T),na.rm =T)
		}else{
			data_min<-y_scale_min
		}
		
	}
}

if(lowCov){
	fileOrder=seq(1,length(fns),by=2)
}else{
	fileOrder=c(1:length(fns))
}

for(j in fileOrder){
	gch<-read.table(fns[j],sep="\t",header=F)	
	#gch[gch[,4]%in%"-",5:length(gch[1,])]=0-gch[gch[,4]%in%"-",5:length(gch[1,])]
	if(smooth == 0){ ## no smooth
		axisSeq<-seq(0-scale, scale, by=step)
		dataSeq<-seq(5+((length(gch[1,])-5)/2)-as.integer(scale/bin_size_align), 5+((length(gch[1,])-5)/2)+as.integer(scale/bin_size_align), by=as.integer(step/bin_size_align))
	}
	if(smooth != 0){ ## smooth
		axisSeq<-seq(0-scale, scale, by=smooth)
		dataSeq<-seq(5+((length(gch[1,])-5)/2)-as.integer(scale/bin_size_align), 5+((length(gch[1,])-5)/2)+as.integer(scale/bin_size_align), by=as.integer(smooth/bin_size_align))
		
	}
	valueGch<-array()
	if(lowCov){
		##the even number of file is to provide the number of total reads..	
		gch_cov<-read.table(fns[j+1],sep="\t",header=F)
		c_cov<-gch*gch_cov
		for(i in dataSeq){
			if(floor(i+step/(2*bin_size_align)-1)-ceiling(i-step/(2*bin_size_align)) <= 0){
				cov_sum<-sum(gch_cov[,ceiling(i-step/(2*bin_size_align))],na.rm=T)
				c_cov_sum<-sum(c_cov[,ceiling(i-step/(2*bin_size_align))],na.rm=T)
				valueGch<-cbind(valueGch,c_cov_sum/cov_sum)
				
			}else{
				cov_sum<-sum(gch_cov[,ceiling(i-step/(2*bin_size_align)):floor(i+step/(2*bin_size_align)-1)],na.rm=T)
				c_cov_sum<-sum(c_cov[,ceiling(i-step/(2*bin_size_align)):floor(i+step/(2*bin_size_align)-1)],na.rm=T)
				valueGch<-cbind(valueGch,c_cov_sum/cov_sum)
				
				#cov_sum<-sum(gch_cov[,ceiling(i-step/(2*bin_size_align)):floor(i+step/(2*bin_size_align)-1)],na.rm=T)
				#c_cov_sum<-sum(c_cov[,ceiling(i-step/(2*bin_size_align)):floor(i+step/(2*bin_size_align)-1)],na.rm=T)
				#valueGch<-cbind(valueGch,c_cov_sum/cov_sum)	
			}
			
		}
		valueGch<-valueGch[,2:length(valueGch[1,])]
		numElemGch<-length(gch[,as.integer(scale/bin_size_align)+5])
		

	}else{
		for(i in dataSeq){
			if(floor(i+step/(2*bin_size_align)-1)-ceiling(i-step/(2*bin_size_align)) <= 0){
				valueGch<-cbind(valueGch,mean(gch[,ceiling(i-step/(2*bin_size_align))], na.rm=T))	
			}else{
				valueGch<-cbind(valueGch,mean(colMeans(gch[,ceiling(i-step/(2*bin_size_align)):floor(i+step/(2*bin_size_align)-1)], na.rm=T), na.rm=T))	
			}
			
		}
		valueGch<-valueGch[,2:length(valueGch[1,])]
		numElemGch<-length(gch[,as.integer(scale/bin_size_align)+5])
	}	
		
	
	
	if(autoScale){
		plot(axisSeq,valueGch,type="l",axes=FALSE,xlab="",ylab="",ylim=c(data_min,data_max),col=colors[j],lty=line_types[j],font=2,lwd=3)
	}else{
		plot(axisSeq,valueGch,type="l",axes=FALSE,xlab="",ylab="",ylim=c(y_scale_min,y_scale_max),col=colors[j],lty=line_types[j],font=2,lwd=3)
	}
	par(new=T)
	
}


#numElemGch<-length(gch[,scale+5])
#numHaveValueElemGch<-length(gch[!is.na(mean(gch[,((length(gch[1,])-4)/2):((length(gch[1,])-4)/2+step-1)], na.rm=T)),][,(length(gch[1,])-4)/2])

#numHaveValueElemHcg<-length(hcg[!is.na(mean(hcg[,dataSeq[length(dataSeq)/2]:dataSeq[length(dataSeq)/2]+step-1], na.rm=T)),][,dataSeq[length(dataSeq)/2]])
mainTitle=paste(prefix,"\nsmooth:", smooth," -numElemCenterToAlign:",numElemGch)
if(autoScale){
	axis(2,at=seq(data_min,data_max,by=(data_max-data_min)/4),labels=format(seq(data_min,data_max,by=(data_max-data_min)/4),digits=3),lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.lab=2,lwd=2)
}else{
	axis(2,at=seq(y_scale_min,y_scale_max,by=y_step),labels=format(seq(y_scale_min,y_scale_max,by=y_step),digits=3),lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.lab=2,lwd=2)
}

axis(1,at=axisSeqForPlot,lty=1,font=2,cex.axis=1.2,cex.lab=1.2,font.lab=2,lwd=2)


title(mainTitle, cex.main = 0.6, font.main= 4, col.main= "black",xlab="Distance to elements (bp)")
if(!is.null(legendNames)){
	legend("topright",legendNames, col=colors,lty=line_types,cex=0.5,lwd=2)
}
abline(v=0)
dev.off()

