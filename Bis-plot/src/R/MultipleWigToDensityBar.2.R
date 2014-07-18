# TODO: Add comment
# 
#base on these two command in R:
#image(y, 1:4, cbind(test,-test,test,-test), col = jet.colors(100), axes = FALSE,zlim=c(-200,200),xlab="UMR",ylab="")
#axis(2, 1:4,c("1","2","3","4"), las = 2, tick = FALSE, font=2) add axis and label at left of y-axis
#axis(4, 1:4,c("1","2","3","4"), las = 2, tick = FALSE, font=2) add axis and label at right side of y-axis
#axis(1, 100,c("1"), las = 1, tick = FALSE, font=2) add triangle at the center.. tick would add tick there..
# Author: yaping
# 2013-7-27
# 
#TO DO: use the right order of experiment list name
#

###############################################################################
library(colorRamps)
library(RColorBrewer)
jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

file_name_lists<-NULL
category_names<-NULL
sample_names<-NULL
experiment_names<-NULL
rep_num_experiment<-NULL  ##should be the same order as wiggle order
addLegendEachExp=FALSE
enrichScoreMax<-NULL
enrichScoreMin=0
logscale=FALSE
capLimitPerc=-1
medianMode=TRUE
for (e in commandArgs(TRUE)) {
	ta = strsplit(e,"=",fixed=TRUE)
	if(! is.na(ta[[1]][2])) {
		if(ta[[1]][1] == "wd"){
			wd<-ta[[1]][2]  ## directory to make plot file
			
		}
		if(ta[[1]][1] == "prefix"){ ##prefix of the output file
			prefix<-ta[[1]][2]
		}
		if(ta[[1]][1] == "file_name_lists"){ ##each file name list should contain a list of experiment in this sample, aligned to one location; the second column should mentioned it is belong to "percentage" or "enrichment" system.
			file_name_lists<-c(file_name_lists,ta[[1]][2])
		}
		if(ta[[1]][1] == "category_names"){ ##number of different location bed files
			category_names<-c(category_names,ta[[1]][2])
		}
		if(ta[[1]][1] == "sample_names"){
			sample_names<-c(sample_names,ta[[1]][2])
		}
		if(ta[[1]][1] == "experiment_names"){
			experiment_names<-c(experiment_names,ta[[1]][2])
		}
		if(ta[[1]][1] == "rep_num_experiment"){
			rep_num_experiment<-c(rep_num_experiment,as.numeric(ta[[1]][2]))
		}
		if(ta[[1]][1] == "step"){
			step<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "scale"){
			scale<-as.numeric(ta[[1]][2])
		}
		
		if(ta[[1]][1] == "bin_size_align"){
			bin_size_align<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "logscale"){
			logscale<-as.logical(ta[[1]][2])
		}
		if(ta[[1]][1] == "addLegendEachExp"){
			addLegendEachExp<-as.logical(ta[[1]][2])
		}
		if(ta[[1]][1] == "capLimitPerc"){
			capLimitPerc<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "enrichScoreMax"){
			enrichScoreMax<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "medianMode"){
			medianMode<-as.logical(ta[[1]][2])
		}

	}
}
setwd(wd)

valueLen<-as.integer((scale/bin_size_align)/(step/bin_size_align))*2+1

dataStore<-array(0, dim = c(length(category_names),length(sample_names),length(experiment_names),valueLen))

enrichScoreScaleMax=-Inf
enrichScoreScaleMin=Inf

percentageExpNum<-NULL
enrichExpNum<-NULL
fileNum=1
for(file_name_list in file_name_lists){
	filesInfo<-read.table(file_name_list,sep="\t",header=F)
	for(fileNumInTable in 1:length(filesInfo[,1])){
		categoryNum=as.integer((fileNum-1)/(length(experiment_names)*length(sample_names)))+1
		sampleNum=as.integer((fileNum-1-(categoryNum-1)*length(experiment_names)*length(sample_names))/length(experiment_names))+1
		#fileNumInTable=ifelse(fileNum %% (length(experiment_names)) == 0 , length(experiment_names) ,fileNum %% (length(experiment_names)))
		valueGch<-NULL
		if(is.na(filesInfo[fileNumInTable,1]) | filesInfo[fileNumInTable,1] %in% ""){
			valueGch=rep(-Inf,valueLen) ##when there is no data, show it is fake no enrichment..
			
			
		}else{
			
			gch<-read.table(as.character(filesInfo[fileNumInTable,1]),sep="\t",header=F)		
			axisSeq<-seq(0-scale, scale, by=step)
			dataSeq<-seq(5+((length(gch[1,])-5)/2)-as.integer(scale/bin_size_align), 5+((length(gch[1,])-5)/2)+as.integer(scale/bin_size_align), by=as.integer(step/bin_size_align))	
			if(capLimitPerc>0 & filesInfo[fileNumInTable,2] %in% "enrichement"){
				if(capLimitPerc>=1){
					stop("capLimitPerc should be between 0<capLimitPerc<1")
				}
				breaks=seq(0,1,capLimitPerc)
				capUplimit=quantile(unlist(gch[,5:length(gch[1,])]),probs=breaks,na.rm =T)[length(breaks)-1]
				gch[gch>capUplimit]=capUplimit
			}
			
			for(i in dataSeq){		
					if(floor(i+step/(2*bin_size_align)-1)-ceiling(i-step/(2*bin_size_align)) <= 0){
						if(medianMode & filesInfo[fileNumInTable,2] %in% "enrichement"){
							valueGch<-cbind(valueGch,median(gch[,ceiling(i-step/(2*bin_size_align))], na.rm=T))
						}else{
							valueGch<-cbind(valueGch,mean(gch[,ceiling(i-step/(2*bin_size_align))], na.rm=T))
						}
							
					}else{
						if(medianMode & filesInfo[fileNumInTable,2] %in% "enrichement"){
							valueGch<-cbind(valueGch,median(gch[,ceiling(i-step/(2*bin_size_align)):floor(i+step/(2*bin_size_align)-1)], na.rm=T))	
						}else{
							valueGch<-cbind(valueGch,mean(colMeans(gch[,ceiling(i-step/(2*bin_size_align)):floor(i+step/(2*bin_size_align)-1)], na.rm=T), na.rm=T))	
						}
						
					}			
				
			}
			if(filesInfo[fileNumInTable,2] %in% "enrichment"){ #not methylation percentage type of data
				#valueGch=ifelse(valueGch<=1,log2(1),log2(valueGch))
				#if(logscale){
				#	valueGch=ifelse(valueGch<=1,log2(1),log2(valueGch))
				#}
				#if(!is.null(enrichScoreMax)){
				#	valueGch[valueGch>enrichScoreMax]=enrichScoreMax
				#}
				
				
				if(max(valueGch,na.rm=T)>enrichScoreScaleMax){
					enrichScoreScaleMax=max(valueGch,na.rm=T)
				}
				if(min(valueGch,na.rm=T)<enrichScoreScaleMin){
					enrichScoreScaleMin=min(valueGch,na.rm=T)
				}
				enrichExpNum<-c(enrichExpNum,fileNumInTable)
			}else if(filesInfo[fileNumInTable,2] %in% "percentage"){
				percentageExpNum<-c(percentageExpNum,fileNumInTable)
			}
			
		}
		dataStore[categoryNum,sampleNum, fileNumInTable,]=valueGch
		fileNum=fileNum+1
	}
	
}

percentageExpNum<-percentageExpNum[!duplicated(percentageExpNum)]

enrichExpNum<-enrichExpNum[!duplicated(enrichExpNum)]
enrichExpNum<-sort(enrichExpNum)

if(!is.null(enrichScoreMax)){
	enrichScoreScaleMax=enrichScoreMax
}
if(!is.null(enrichScoreMin)){
	enrichScoreScaleMin=enrichScoreMin
}

#####make the distribution of height and width in each block
#style<-cbind(c(1,2),2+matrix(1:(length(sample_names)*2), 2, length(sample_names), byrow = TRUE))
numEqualFracInRows=length(experiment_names)+3
numEqualFracInCols=1+ifelse(addLegendEachExp, 2*length(sample_names)*length(category_names), 2*(length(sample_names)*length(category_names)-1)+4)
numRows=length(rep_num_experiment)+1
numCols=length(sample_names)*length(category_names)+1

titleFrac=3/numEqualFracInRows
ExpFrac=1/numEqualFracInRows


barFrac=1/numEqualFracInCols
lastSampleFrac=ifelse(addLegendEachExp,2/numEqualFracInCols,3/numEqualFracInCols)
otherSampleFrac=2/numEqualFracInCols

lht<-c(titleFrac,rep_num_experiment*ExpFrac)
lwd<-c(barFrac,rep(otherSampleFrac,(length(sample_names)*length(category_names)-1)),lastSampleFrac)


#####make the block number in canvas
style<-matrix(1:(numRows*numCols),  numRows, numCols, byrow = FALSE)
j=1
style[1:length(style[,1]),1]=j
j=j+1

for(i in 2+(length(style[1,])-1)/length(category_names)*(0:(length(category_names)-1))){	
	

		for(z in i:(i+length(sample_names)-1)){	
			style[1,z]=j
		}
		j=j+1
		
		for(z in i:(i+length(sample_names)-1)){	
			for(r in 2:(length(rep_num_experiment)+1)){							
				style[r,z]=j
				j=j+1
			}
			
		}

}


pdf(paste(prefix,"pdf",sep="."), paper="special", height=0.3*numEqualFracInRows, width=1.5*numCols)
layout(style,heights = lht,widths=lwd)

##plot scale bar 
par(mar=c(3, 2.2, 3, 2.2),oma=c(0,1.5,0,1.5))
image(1, seq(enrichScoreScaleMin,enrichScoreScaleMax,by=(enrichScoreScaleMax-enrichScoreScaleMin)/100), matrix(seq(enrichScoreScaleMin,enrichScoreScaleMax,by=(enrichScoreScaleMax-enrichScoreScaleMin)/100), 1, 101, byrow = TRUE), col = jet.colors(100), axes = FALSE,zlim=c(enrichScoreScaleMin,enrichScoreScaleMax),xlab="",ylab="",font.lab=2, cex.lab=1.0)
axis(2, seq(enrichScoreScaleMin,enrichScoreScaleMax,by=(enrichScoreScaleMax-enrichScoreScaleMin)/5), labels=seq(0,100,by=20), las = 2, tick = T, font=2, cex.axis=1)
axis(4, seq(enrichScoreScaleMin,enrichScoreScaleMax,by=(enrichScoreScaleMax-enrichScoreScaleMin)/5), labels=format(seq(enrichScoreScaleMin,enrichScoreScaleMax,by=(enrichScoreScaleMax-enrichScoreScaleMin)/5),digits=1), las = 2, tick = T, font=2, cex.axis=1, cex.lab=1)
mtext("Methylation/Accessibility",side = 2, adj=0.5, outer=T)
mtext(ifelse(logscale,"log2(Z-score)","Z-score"),side = 4, adj=0.5, outer=F)
#plot.new()
#percentageExpNum<-rev(percentageExpNum)

for(categoryNum in 1:length(category_names)){
	
	plot.new()
	##plot title
	par(mar=c(0.5, 0.5, 0.5, 0.5),usr=c(0,1,0,1))
	halfHeight=(par("usr")[4]-par("usr")[3])/2
	fourthHeight=(par("usr")[4]-par("usr")[3])/4
	halfWidth=(par("usr")[2]-par("usr")[1])/2
	equalWidth=(par("usr")[2]-par("usr")[1])/(length(sample_names)*2)
	if(!addLegendEachExp && categoryNum == length(category_names)){
		equalWidth=equalWidth*4/5
	}
	
	text(halfWidth,fourthHeight*3, category_names[categoryNum], font=2, cex=1.5)
	segments(par("usr")[1],halfHeight,par("usr")[2],halfHeight,lwd=2)
	for(sampleNum in 1:length(sample_names)){
		text(sampleNum*equalWidth*2-equalWidth,fourthHeight, labels=sample_names[sampleNum], font=2, cex=1.5, adj=0.5)
	}
	
	
	##then need to figure out the following..
	 
	
	for(sampleNum in 1:length(sample_names)){
		#par(oma=c(0, 0, 0, 0),mar=c(0, 1, 2, 2))
		annoLegend=FALSE
		if((categoryNum == length(category_names) && sampleNum == length(sample_names)) || addLegendEachExp){
			annoLegend=TRUE
		}
		experimentOrder=0
		for(i in rep_num_experiment){
			experimentSeq = (experimentOrder+1):(experimentOrder+i)
			experimentSeq<-rev(experimentSeq)
			inPercentageGroup= ifelse(length(percentageExpNum[percentageExpNum==experimentSeq[1]])>0,TRUE,FALSE)
			dataTmp=as.matrix(dataStore[categoryNum,sampleNum,experimentSeq,])
			if(dim(dataTmp)[2]>1){
				dataTmp=t(dataTmp)
			}
			if(inPercentageGroup){
				limit = c(0,100)
			}else{
				limit = c(enrichScoreScaleMin,enrichScoreScaleMax)
				if(logscale){
					dataTmp=ifelse(dataTmp<=1,log2(1),log2(dataTmp))
				}
				if(!is.null(enrichScoreMax)){
					dataTmp[dataTmp>enrichScoreMax]=enrichScoreMax
				}
				if(!is.null(enrichScoreMin)){
					dataTmp[(!is.infinite(dataTmp)) & dataTmp<enrichScoreMin]=enrichScoreMin
				}
			}
			par(mar=c(0.5, 0.5, 0, ifelse(annoLegend,8,1)))
			image(1:valueLen, 1:length(experimentSeq), dataTmp, col = jet.colors(100), axes = FALSE,zlim=limit,xlab="",ylab="")
			axis(1, valueLen/2, "", las = 2, tick = TRUE, font=2,cex.axis=2)
			if(annoLegend){
				axis(4, 1:(length(experimentSeq)),experiment_names[experimentSeq], las = 2, tick = FALSE, font=2,cex.axis=1.5,cex.lab=2.5)	
			}
			experimentOrder = experimentOrder+i
		}
		
		

		
	}
	
}
dev.off()


