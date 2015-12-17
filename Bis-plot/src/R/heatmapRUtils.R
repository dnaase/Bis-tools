# heatmapRUtils.R
# 
# Author: yaping
# Oct 1, 2013
# 4:05:31 PM
###############################################################################
library(fastcluster)
library(ape)
##################
#make tree face to the left side need 'plotPhylo' function in "ape" package
library(dendextend)
##################
#rotate tree need 'rotate' function in "dendextend" package, which is downloaded from here: https://github.com/talgalili/dendextend/releases
#download the whole package, and then  R CMD INSTALL dendextend.*.tar.gz
#hc_new<-rotate(hc,c(subClusterOrder[subClusterOrder==1], subClusterOrder[subClusterOrder==2], subClusterOrder[subClusterOrder==3], subClusterOrder[subClusterOrder==4]))
#plot(hc_new)

###function to plot one NOMe-seq clustering heatmap + multiple ChIP-seq heatmap
plotMultiWigToHeatmapPlusAveInSingleLoc<-function(files, sampleNames, prefix="test",capUpLimit=NULL, capDownLimit=NULL,fileNumToOrder=1, regionToCluster=c(-500,500), numCluster=2, bin_size_align=20, notPlotIndexMatrix=FALSE,
		bin_size=20, scale=5000, heatmap_clustering_bin_size_align=1, heatmap_clustering_bin_size=20, heatmap_clustering_scale=1000, capLimit=NULL, logScale=NULL, heatMapCols=NULL, addAverage=addAverage, order_dendgram=NULL){##return limit boundary of the data
	##read index matrix, do clustering
	orderStat<-clusteringMultiMatrixNew(files[fileNumToOrder],regionToCluster=regionToCluster, numCluster=numCluster, bin_size_align=heatmap_clustering_bin_size_align, bin_size=heatmap_clustering_bin_size, scale=heatmap_clustering_scale, order_dendgram=order_dendgram)
	
	##read data, follow the order index matrix give
	
	if(!is.null(logScale)){
		#logScale=logScale[-fileNumToOrder]
	}
	dataSummaryFrame<-readMultiMatrix(files[-fileNumToOrder], sampleNames[-fileNumToOrder],rowOrder=orderStat$rowOrder, capUpLimit=capUpLimit[-fileNumToOrder], capDownLimit=capDownLimit[-fileNumToOrder],numCluster=numCluster,
			subClusterOrder=orderStat$subClusterOrder, doClustering=F, capLimit=capLimit[-fileNumToOrder], logScale=logScale[-fileNumToOrder], bin_size_align=bin_size_align, bin_size=bin_size,scale=scale) 
	##get a data frame, 1st is category name, 2nd is average data in each category, 3rd is data matrix ordered in each category, 4th is hclust tree, 5th is value limit
	##name, average, matrix, tree, limit
	
	
	##layout
	#####make the distribution of height and width in each block
	
	numEqualFracInRows=ifelse(addAverage,9.5,7.5) #1:6:2 or when no average plot, 1:6
	numEqualFracInCols=1+0.1+3*length(fileNumToOrder)+length(files)-length(fileNumToOrder) #tree, sideBar, index matrix, chip-seq sample
	numCols=2+length(files)
	numRows=ifelse(addAverage,3,2)
	
	titleFrac=1.5/numEqualFracInRows
	heatmapFrac=6/numEqualFracInRows
	averageFrac=2/numEqualFracInRows
	
	
	treeFrac=1/numEqualFracInCols
	sideBarFrac=0.1/numEqualFracInCols
	indexFrac=1/numEqualFracInCols
	sampleFrac=1/numEqualFracInCols
	
	
	if(addAverage){
		lht<-c(titleFrac,heatmapFrac,averageFrac)
	}else{
		lht<-c(titleFrac,heatmapFrac)
	}
	lwd<-c(treeFrac,sideBarFrac,rep(indexFrac,length(fileNumToOrder)), rep(sampleFrac,(length(files)-length(fileNumToOrder))))
	
	
	#####make the block number in canvas
	style<-matrix(1:(numRows*numCols),  numRows, numCols, byrow = FALSE)
#if(addAverage){
#	style[3,1]=0
#	style[,2:numCols] <- style[,2:numCols]-1
#}
#if(!doClustering){
#	style[2,1]=0
#	style[,2:numCols] <- style[,2:numCols]-1
#}
	
	#pdf(paste(prefix,"pdf",sep="."), paper="special", height=1*numEqualFracInRows, width=1.7*numEqualFracInCols)
	jpeg(paste(prefix,"jpg",sep="."), height=1*numEqualFracInRows, width=1.7*numEqualFracInCols, units = "in",res=72)
	layout(style,heights = lht,widths=lwd)
	
	############plot key
	plotKey(dataSummaryFrame$limit)
	
	
	###########plot tree
	if(!is.null(orderStat$subClusterOrder)){
		plotTree(as.dendrogram(orderStat$tree))
	}else plot.new()
	
	if(addAverage){
		plot.new()	
	}
	
	###########plot side bar
	plot.new()
	if(!is.null(orderStat$subClusterOrder)){
		plotSideBar(orderStat$subClusterOrder,orderStat$rowOrder)
	}else plot.new()
	
	if(addAverage){
		plot.new()	
	}
	
	###########plot index matrix


		##read index matrix
		#dataSummaryFrameIndexMatrix<-readMultiMatrix(files[fileNumToOrder], sampleNames[fileNumToOrder],rowOrder=orderStat$rowOrder, capUpLimit=capUpLimit, capDownLimit=capDownLimit,numCluster=numCluster,
		#		subClusterOrder=orderStat$subClusterOrder, doClustering=F, capLimit=capLimit[fileNumToOrder], logScale=logScale[fileNumToOrder], bin_size_align=heatmap_clustering_bin_size_align, bin_size=heatmap_clustering_bin_size,scale=heatmap_clustering_scale) 
	for(index in fileNumToOrder){
		##read index matrix
		dataSummaryFrameIndexMatrix<-readMultiMatrix(files[index], sampleNames[index],rowOrder=orderStat$rowOrder, capUpLimit=capUpLimit[index], capDownLimit=capDownLimit[index],numCluster=numCluster,
				subClusterOrder=orderStat$subClusterOrder, doClustering=F, capLimit=capLimit[index], logScale=logScale[index], bin_size_align=heatmap_clustering_bin_size_align, bin_size=heatmap_clustering_bin_size,scale=heatmap_clustering_scale) 
		
		
		##plot title
		plot.new()
		par(mar=c(0.5, 0.5, 0.5, 0.5),usr=c(0,1,0,1))
		halfHeight=(par("usr")[4]-par("usr")[3])/2
		halfWidth=(par("usr")[2]-par("usr")[1])/2
		text(halfWidth,halfHeight, sampleNames[index], font=2, cex=4)
		
		##plot heatmap
		#plotHeatmap(orderStat$matrix, categoryNum=1, addNumber=T, xAxisScale=c(-heatmap_clustering_scale,heatmap_clustering_scale), heatMapCols=heatMapCols[1])
		plotHeatmap(dataSummaryFrameIndexMatrix$matrix, 1, addNumber=F, xAxisScale=c(-heatmap_clustering_scale,heatmap_clustering_scale), heatMapCols=heatMapCols[index])
		
		##plot average plot
		if(addAverage){
			#plotAverage(orderStat$average, 1,  xAxisScale=c(-heatmap_clustering_scale,heatmap_clustering_scale),xStep=heatmap_clustering_scale/2)
			plotAverage(dataSummaryFrameIndexMatrix$average, 1, xAxisScale=c(-heatmap_clustering_scale,heatmap_clustering_scale),xStep=heatmap_clustering_scale/2, yAxisScale=c(capDownLimit[index],capUpLimit[index]),yStep=(capUpLimit[index]-capDownLimit[index])/2)
		}
	}
	
	
	###############plot value
	#heatMapCols<-heatMapCols[-fileNumToOrder]
	for(i in 1:length(dataSummaryFrame$name)){
		tmpColor<-heatMapCols[-fileNumToOrder]
		tmpUpLimit<-capUpLimit[-fileNumToOrder]
		tmpDownLimit<-capDownLimit[-fileNumToOrder]
		##plot title
		plot.new()
		par(mar=c(0.5, 0.5, 0.5, 0.5),usr=c(0,1,0,1))
		halfHeight=(par("usr")[4]-par("usr")[3])/2
		halfWidth=(par("usr")[2]-par("usr")[1])/2
		text(halfWidth,halfHeight, dataSummaryFrame$name[i], font=2, cex=2.5, adj=0.5,srt = 45, xpd = TRUE)
		
		##plot heatmap
		plotHeatmap(dataSummaryFrame$matrix, i, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=tmpColor[i])
		
		##plot average plot
		if(addAverage){
			plotAverage(dataSummaryFrame$average, i, xAxisScale=c(-scale,scale),xStep=scale,yAxisScale=c(tmpDownLimit[i],tmpUpLimit[i]),yStep=(tmpUpLimit[i]-tmpDownLimit[i])/2)
		}
		
	}
	
	dev.off()
	
	
	################make pdf version:
	
	pdf(paste(prefix,"pdf",sep="."), paper="special", height=1*numEqualFracInRows, width=1.7*numEqualFracInCols)
	#jpeg(paste(prefix,"jpg",sep="."), height=1*numEqualFracInRows, width=1.7*numEqualFracInCols, units = "in",res=72)
	layout(style,heights = lht,widths=lwd)
	
	############plot key
	plotKey(dataSummaryFrame$limit)
	
	
	###########plot tree
	if(!is.null(orderStat$subClusterOrder)){
		plotTree(as.dendrogram(orderStat$tree))
	}else plot.new()
	
	if(addAverage){
		plot.new()	
	}
	
	###########plot side bar
	plot.new()
	if(!is.null(orderStat$subClusterOrder)){
		plotSideBar(orderStat$subClusterOrder,orderStat$rowOrder)
	}else plot.new()
	
	if(addAverage){
		plot.new()	
	}
	
	###########plot index matrix
	
	
	##read index matrix
	#dataSummaryFrameIndexMatrix<-readMultiMatrix(files[fileNumToOrder], sampleNames[fileNumToOrder],rowOrder=orderStat$rowOrder, capUpLimit=capUpLimit, capDownLimit=capDownLimit,numCluster=numCluster,
	#		subClusterOrder=orderStat$subClusterOrder, doClustering=F, capLimit=capLimit[fileNumToOrder], logScale=logScale[fileNumToOrder], bin_size_align=heatmap_clustering_bin_size_align, bin_size=heatmap_clustering_bin_size,scale=heatmap_clustering_scale) 
	for(index in fileNumToOrder){
		dataSummaryFrameIndexMatrix<-readMultiMatrix(files[index], sampleNames[index],rowOrder=orderStat$rowOrder, capUpLimit=capUpLimit[index], capDownLimit=capDownLimit[index],numCluster=numCluster,
				subClusterOrder=orderStat$subClusterOrder, doClustering=F, capLimit=capLimit[index], logScale=logScale[index], bin_size_align=heatmap_clustering_bin_size_align, bin_size=heatmap_clustering_bin_size,scale=heatmap_clustering_scale) 
		
		
		##plot title
		plot.new()
		par(mar=c(0.5, 0.5, 0.5, 0.5),usr=c(0,1,0,1))
		halfHeight=(par("usr")[4]-par("usr")[3])/2
		halfWidth=(par("usr")[2]-par("usr")[1])/2
		text(halfWidth,halfHeight, sampleNames[index], font=2, cex=4)
		
		##plot heatmap
		#plotHeatmap(orderStat$matrix, categoryNum=1, addNumber=T, xAxisScale=c(-heatmap_clustering_scale,heatmap_clustering_scale), heatMapCols=heatMapCols[1])
		plotHeatmap(dataSummaryFrameIndexMatrix$matrix, 1, addNumber=F, xAxisScale=c(-heatmap_clustering_scale,heatmap_clustering_scale), heatMapCols=heatMapCols[index])
		
		##plot average plot
		if(addAverage){
			#plotAverage(orderStat$average, 1,  xAxisScale=c(-heatmap_clustering_scale,heatmap_clustering_scale),xStep=heatmap_clustering_scale/2)
			plotAverage(dataSummaryFrameIndexMatrix$average, 1, xAxisScale=c(-heatmap_clustering_scale,heatmap_clustering_scale),xStep=heatmap_clustering_scale/2,yAxisScale=c(capDownLimit[index],capUpLimit[index]),yStep=(capUpLimit[index]-capDownLimit[index])/2)
		}
	}
	
	
	###############plot value
	#heatMapCols<-heatMapCols[-fileNumToOrder]
	for(i in 1:length(dataSummaryFrame$name)){
		tmpColor<-heatMapCols[-fileNumToOrder]
		tmpUpLimit<-capUpLimit[-fileNumToOrder]
		tmpDownLimit<-capDownLimit[-fileNumToOrder]
		##plot title
		plot.new()
		par(mar=c(0.5, 0.5, 0.5, 0.5),usr=c(0,1,0,1))
		halfHeight=(par("usr")[4]-par("usr")[3])/2
		halfWidth=(par("usr")[2]-par("usr")[1])/2
		text(halfWidth,halfHeight, dataSummaryFrame$name[i], font=2, cex=2.5, adj=0.5,srt = 45, xpd = TRUE)
		
		##plot heatmap
		plotHeatmap(dataSummaryFrame$matrix, i, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=tmpColor[i])
		
		##plot average plot
		if(addAverage){
			plotAverage(dataSummaryFrame$average, i, xAxisScale=c(-scale,scale),xStep=scale, yAxisScale=c(tmpDownLimit[i],tmpUpLimit[i]),yStep=(tmpUpLimit[i]-tmpDownLimit[i])/2)
		}
		
	}
	
	dev.off()
	
	
	dataSummaryFrame$limit
}


##function to plot multiple heatmap for NOMe-seq signal
plotMultiWigToHeatmapPlusAveInSingleLocForNOMeSeqPaper<-function(files, sampleNames, prefix="test",limit=NULL, fileNumToOrder=5, regionToCluster=c(-500,500), numCluster=2, bin_size_align=1,
		bin_size=20, scale=1000, heatmap_clustering_bin_size_align=1, heatmap_clustering_bin_size=20, heatmap_clustering_scale=1000, heatMapCols=c(blue2yellow,blue2yellow,blue2yellow,white2darkgreen,white2darkgreen), addAverage=TRUE, 
		ylim=cbind(c(0,100),c(0,100),c(0,100),c(20,80),c(20,80))){##return NULL
	##read index matrix, do clustering
	orderStat<-clusteringSingleMatrix(files[fileNumToOrder],regionToCluster=regionToCluster, numCluster=numCluster, bin_size_align=heatmap_clustering_bin_size_align, bin_size=heatmap_clustering_bin_size, scale=heatmap_clustering_scale)
	

	dataSummaryFrame<-readMultiMatrix(files, sampleNames,rowOrder=orderStat$rowOrder, numCluster=numCluster,
			subClusterOrder=orderStat$subClusterOrder, doClustering=F,  bin_size_align=bin_size_align, bin_size=bin_size,scale=scale) 
	##get a data frame, 1st is category name, 2nd is average data in each category, 3rd is data matrix ordered in each category, 4th is hclust tree, 5th is value limit
	##name, average, matrix, tree, limit
	
	
	##layout
	#####make the distribution of height and width in each block
	
	numEqualFracInRows=ifelse(addAverage,9,6) #6:2 or when no average plot, 6
	numEqualFracInCols=3*length(files)+0.7+0.15+1 #methylation matrix, insualter, accessbility matrix, sidebar, tree
	numCols=length(files)+1+1+1
	numRows=ifelse(addAverage,2,1)
	
	#titleFrac=1/numEqualFracInRows
	heatmapFrac=6/numEqualFracInRows
	averageFrac=3/numEqualFracInRows
	
	
	#
	sideBarFrac=0.15/numEqualFracInCols
	#indexFrac=3/numEqualFracInCols
	insulaterFrac=0.7/numEqualFracInCols
	sampleFrac=3/numEqualFracInCols
	treeFrac=1/numEqualFracInCols
	
	if(addAverage){
		lht<-c(heatmapFrac,averageFrac)
	}else{
		lht<-c(heatmapFrac)
	}
	lwd<-c(rep(sampleFrac,3),insulaterFrac, rep(sampleFrac,2),sideBarFrac,treeFrac)
	
	
	#####make the block number in canvas
	style<-matrix(1:(numRows*numCols),  numRows, numCols, byrow = FALSE)

	
	#pdf(paste(prefix,"pdf",sep="."), paper="special", height=1*numEqualFracInRows, width=1.7*numEqualFracInCols)
	jpeg(paste(prefix,"jpg",sep="."), height=1.0*numEqualFracInRows, width=1.7*numEqualFracInCols, units = "in",res=300)
	layout(style,heights = lht,widths=lwd)

	
	###############plot metyhylation value

	for(i in 1:3){

		##plot heatmap
		plotHeatmap(dataSummaryFrame$matrix, i, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[i], label=FALSE)
		
		##plot average plot
		if(addAverage){
			plotAverage(dataSummaryFrame$average, i, xAxisScale=c(-scale,scale),xStep=scale/2, yAxisScale=ylim[,i],yStep=(ylim[2,i]-ylim[1,i])/2, label=FALSE)
		}
		
	}
	
	###############plot insulater value
	plot.new()	
	if(addAverage){
		plot.new()	
	}
	###############plot accessbility value
	
	for(i in 4:5){
		
		##plot heatmap
		plotHeatmap(dataSummaryFrame$matrix, i, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[i], label=FALSE)
		
		##plot average plot
		if(addAverage){
			plotAverage(dataSummaryFrame$average, i, xAxisScale=c(-scale,scale),xStep=scale/2, yAxisScale=ylim[,i],yStep=(ylim[2,i]-ylim[1,i])/2, label=FALSE)
		}
		
	}
	
	###########plot side bar
	if(!is.null(orderStat$subClusterOrder)){
		plotSideBar(orderStat$subClusterOrder,orderStat$rowOrder)
	}else plot.new()
	
	if(addAverage){
		par(mar=c(0,0,0,0))
		plot.new()	
	}
	
	###########plot tree
	if(!is.null(orderStat$subClusterOrder)){
		plotPhyloTree(as.phylo(orderStat$tree), direction="leftwards")
	}else plot.new()
	
	if(addAverage){
		plot.new()	
	}
	
	dev.off()

	###PRINT PDF VERSION
	#pdf(paste(prefix,"pdf",sep="."), paper="special", height=1.0*numEqualFracInRows, width=1.7*numEqualFracInCols)
	
	
	#dev.off()
	
}


##function to plot multiple heatmap for NOMe-seq signal
plotMultiWigToHeatmapPlusAveInMultipleLocForNOMeSeqPaper<-function(files, sampleNames, prefix="test",fileNumToOrder=5, regionToCluster=c(-500,500), numCluster=c(3,2,2,2), bin_size_align=1, breaks=4, notPlotIndexMatrix=FALSE,
		capLimit=NULL, capUpLimit=NULL, capDownLimit=NULL, logScale=NULL, motif=NULL, adjust_len_proportion=FALSE,
		bin_size=20, scale=1000, heatmap_clustering_bin_size_align=1, heatmap_clustering_bin_size=20, heatmap_clustering_scale=1000, heatMapCols=c(blue2yellow,blue2yellow,blue2yellow,white2darkgreen,white2darkgreen), addAverage=TRUE, 
		ylim=cbind(c(0,100),c(0,100),c(0,100),c(0,100)), heatmap_anno=NULL, heatmap_anno_col=c(white2black), heatmap_anno_name="test", order_dendgram=NULL, heatmap_autoscale=NULL){##return NULL
	numLoc=as.integer(length(files)/breaks)
	##layout
	#####make the distribution of height and width in each block	
	numEqualFracInRows=ifelse(addAverage,9,6)*numLoc #6:2 or when no average plot, 6
	numEqualFracInCols=3*breaks+0.7*1+0.15*ifelse(is.null(heatmap_anno),0,length(heatmap_anno)+1)+1 #methylation matrix, insualter, accessbility matrix, sidebar, tree
	
	#numCols=breaks+1+1+1
	numCols=breaks+1*(breaks-2)/2+1*ifelse(is.null(heatmap_anno),0,1) + 1
	numRows=ifelse(addAverage,2,1)*numLoc
	
	if(notPlotIndexMatrix){
		numEqualFracInCols = numEqualFracInCols - 3*length(fileNumToOrder)
		numCols = numCols - length(fileNumToOrder)
	}
	
	#titleFrac=1/numEqualFracInRows
	heatmapFrac=6/numEqualFracInRows
	averageFrac=3/numEqualFracInRows

	#
	sideBarFrac=ifelse(is.null(heatmap_anno),0,length(heatmap_anno)+1) * 0.15/numEqualFracInCols
	#indexFrac=3/numEqualFracInCols
	insulaterFrac=0.7/numEqualFracInCols
	sampleFrac=3/numEqualFracInCols
	treeFrac=1/numEqualFracInCols
	
	if(adjust_len_proportion){
		proportion<-file_proportion(files,breaks,bin_size_align=bin_size_align,bin_size=bin_size,scale=scale)
		lht<-proportion$lht
		numEqualFracInRows<-proportion$numEqualFracInRows
	}else if(addAverage){
		lht<-rep(c(heatmapFrac,averageFrac),numLoc)
	}else{
		lht<-rep(c(heatmapFrac),numLoc)
	}
	#lwd<-c(rep(sampleFrac,2),insulaterFrac, rep(sampleFrac,2),sideBarFrac,treeFrac)
	#lwd<-c(rep(sampleFrac,2),insulaterFrac, rep(sampleFrac,2),treeFrac)
	lwd<-c(rep(sampleFrac,2),rep(c(insulaterFrac, rep(sampleFrac,2)),(breaks-2)/2),treeFrac)
	if(!is.null(heatmap_anno)){
		lwd<-c(sideBarFrac, rep(sampleFrac,2),rep(c(insulaterFrac, rep(sampleFrac,2)),(breaks-2)/2),treeFrac)
		#lwd<-c(sideBarFrac,rep(sampleFrac,2),insulaterFrac, rep(sampleFrac,2),treeFrac)
	}

	
	#####make the block number in canvas
	style<-matrix(1:(numRows*numCols),  numRows, numCols, byrow = TRUE)
	print(style)
	print(lwd)
	print(lht)
	#pdf(paste(prefix,"pdf",sep="."), paper="special", height=1*numEqualFracInRows, width=1.7*numEqualFracInCols)
	jpeg(paste(prefix,"jpg",sep="."), height=1.0*numEqualFracInRows, width=1.7*numEqualFracInCols, units = "in",res=300)
	layout(style,heights = lht,widths=lwd)
	
	for(i in seq(1,length(files),by=breaks)){
		start<-i
		end<-i+breaks-1
		filesInOneLoc<-files[start:end]
		##read index matrix, do clustering
		#orderStat<-clusteringSingleMatrix(files[(start+fileNumToOrder-1)],regionToCluster=regionToCluster, numCluster=numCluster[start], bin_size_align=heatmap_clustering_bin_size_align, bin_size=heatmap_clustering_bin_size, scale=heatmap_clustering_scale)
		orderStat<-clusteringMultiMatrixNew(files[(start+fileNumToOrder-1)],regionToCluster=regionToCluster, numCluster=numCluster[start], bin_size_align=heatmap_clustering_bin_size_align, bin_size=heatmap_clustering_bin_size, scale=heatmap_clustering_scale, order_dendgram=order_dendgram)
		
		
		##get a data frame, 1st is category name, 2nd is average data in each category, 3rd is data matrix ordered in each category, 4th is hclust tree, 5th is value limit
		##name, average, matrix, tree, limit
		#dataSummaryFrame<-readMultiMatrix(filesInOneLoc, sampleNames[start:end],rowOrder=orderStat$rowOrder, numCluster=numCluster[start],capLimit=c(FALSE,FALSE,FALSE,TRUE,TRUE), capUpLimit=80, capDownLimit=20,
		dataSummaryFrame<-readMultiMatrix(filesInOneLoc, sampleNames[start:end],rowOrder=orderStat$rowOrder, numCluster=numCluster[start],capLimit=capLimit[start:end],	capUpLimit=capUpLimit, capDownLimit=capDownLimit,
				logScale=logScale[start:end], subClusterOrder=orderStat$subClusterOrder, doClustering=F,  bin_size_align=bin_size_align, bin_size=bin_size,scale=scale) 

		##plot sidebar
		if(!is.null(heatmap_anno)){
				plotMultipleSideBar(orderStat$subClusterOrder,orderStat$rowOrder, heatmap_anno=heatmap_anno, heatmap_anno_col=heatmap_anno_col, heatmap_anno_name=heatmap_anno_name)
		}
		
		##plot heatmap
		if(notPlotIndexMatrix){
			#for(j in fileNumToOrder){
			#	plotHeatmap(dataSummaryFrame$matrix, j, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[1], label=FALSE)
			#}
			#if(length(fileNumToOrder)>0){
			#	par(mar=c(0,0,0,0))
			#	plot.new()
			#}
			for(j in seq(length(fileNumToOrder)+1, end-start+1, by=2)){
				plotHeatmap(dataSummaryFrame$matrix, j, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[j], label=FALSE, adjust_len_proportion=adjust_len_proportion)
				if((j+1) > breaks){
					break;	
				}
				plotHeatmap(dataSummaryFrame$matrix, j+1, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[j+1], label=FALSE, adjust_len_proportion=adjust_len_proportion)
				if((j+1) != breaks){
					par(mar=c(0,0,0,0))
					plot.new()
				}	
			}
		}else{
			plotHeatmap(dataSummaryFrame$matrix, 1, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[1], label=FALSE, adjust_len_proportion=adjust_len_proportion,colLim=c(0,100))
			plotHeatmap(dataSummaryFrame$matrix, 2, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[2], label=FALSE, adjust_len_proportion=adjust_len_proportion,colLim=c(0,100))
			par(mar=c(0,0,0,0))
			plot.new()
			for(j in seq(3, end-start+1, by=2)){
				plotHeatmap(dataSummaryFrame$matrix, j, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[j], label=FALSE, adjust_len_proportion=adjust_len_proportion)
				if((j+1) > breaks){
					break;	
				}
				plotHeatmap(dataSummaryFrame$matrix, j+1, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[j+1], label=FALSE, adjust_len_proportion=adjust_len_proportion)
				if((j+1) != breaks){
					par(mar=c(0,0,0,0))
					plot.new()
				}
			}
		}
		
		
		###############plot metyhylation value
	#	for(i in 1:2){			
	#		plotHeatmap(dataSummaryFrame$matrix, i, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[i], label=FALSE)
	#	}
		###############plot insulater value
	#	par(mar=c(0,0,0,0))
	#	plot.new()	
		###############plot accessbility value
	#	for(i in 3:4){
	#		plotHeatmap(dataSummaryFrame$matrix, i, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[i], label=FALSE,colLim=c(0,100))
	#	}
		
		###########plot side bar
		#if(!is.null(orderStat$subClusterOrder)){
			#plotMultipleSideBar(orderStat$subClusterOrder,orderStat$rowOrder, heatmap_anno=heatmap_anno, heatmap_anno_col=heatmap_anno_col, heatmap_anno_name=heatmap_anno_name)
		#}else plot.new()
		
		###########plot tree
		if(!is.null(orderStat$subClusterOrder)){
			plotPhyloTree(as.phylo(orderStat$tree), direction="leftwards")
		}else plot.new()
		
		
		if(addAverage){
			if(!is.null(heatmap_anno)){
				plot.new()
			}
			plotAverage(dataSummaryFrame$average, 1, xAxisScale=c(-scale,scale),xStep=scale/2, yAxisScale=ylim[,1],yStep=(ylim[2,1]-ylim[1,1])/2, label=TRUE, heatmap_autoscale=heatmap_autoscale[1])
			plotAverage(dataSummaryFrame$average, 2, xAxisScale=c(-scale,scale),xStep=scale/2, yAxisScale=ylim[,2],yStep=(ylim[2,2]-ylim[1,2])/2, label=TRUE, heatmap_autoscale=heatmap_autoscale[2])
			for(j in seq(3, end-start+1, by=2)){
				par(mar=c(0,0,0,0))
				plot.new()	
				plotAverage(dataSummaryFrame$average, j, xAxisScale=c(-scale,scale),xStep=scale/2, yAxisScale=ylim[,j],yStep=(ylim[2,j]-ylim[1,j])/2, label=TRUE, heatmap_autoscale=heatmap_autoscale[j])	
				if((j+1) > breaks){
					break;	
				}
				plotAverage(dataSummaryFrame$average, j+1, xAxisScale=c(-scale,scale),xStep=scale/2, yAxisScale=ylim[,j+1],yStep=(ylim[2,j+1]-ylim[1,j+1])/2, label=TRUE, heatmap_autoscale=heatmap_autoscale[j+1])
			}
			
			
			###############plot metyhylation value
		#	for(i in 1:2){	
		#		#print(ylim)
		#		plotAverage(dataSummaryFrame$average, i, xAxisScale=c(-scale,scale),xStep=scale/2, yAxisScale=ylim[,i],yStep=(ylim[2,i]-ylim[1,i])/2, label=FALSE)	
		#	}
			###############plot insulater value
		#	par(mar=c(0,0,0,0))
		#	plot.new()	
			###############plot accessbility value
		#	for(i in 3:4){
		#		#print(ylim)
		#		plotAverage(dataSummaryFrame$average, i, xAxisScale=c(-scale,scale),xStep=scale/2, yAxisScale=ylim[,i],yStep=(ylim[2,i]-ylim[1,i])/2, label=FALSE)
		#	}
			##############rest of vacant value
			par(mar=c(0,0,0,0))
			plot.new()	
		}

	}
	
	dev.off()
	
	###PRINT PDF VERSION
	pdf(paste(prefix,"pdf",sep="."), paper="special", height=1.0*numEqualFracInRows, width=1.7*numEqualFracInCols)
	layout(style,heights = lht,widths=lwd)
	
	for(i in seq(1,length(files),by=breaks)){
		start<-i
		end<-i+breaks-1
		filesInOneLoc<-files[start:end]
		##read index matrix, do clustering
		#orderStat<-clusteringSingleMatrix(files[(start+fileNumToOrder-1)],regionToCluster=regionToCluster, numCluster=numCluster[start], bin_size_align=heatmap_clustering_bin_size_align, bin_size=heatmap_clustering_bin_size, scale=heatmap_clustering_scale)
		orderStat<-clusteringMultiMatrixNew(files[(start+fileNumToOrder-1)],regionToCluster=regionToCluster, numCluster=numCluster[start], bin_size_align=heatmap_clustering_bin_size_align, bin_size=heatmap_clustering_bin_size, scale=heatmap_clustering_scale, order_dendgram=order_dendgram)
		
		
		##get a data frame, 1st is category name, 2nd is average data in each category, 3rd is data matrix ordered in each category, 4th is hclust tree, 5th is value limit
		##name, average, matrix, tree, limit
		#dataSummaryFrame<-readMultiMatrix(filesInOneLoc, sampleNames[start:end],rowOrder=orderStat$rowOrder, numCluster=numCluster[start],capLimit=c(FALSE,FALSE,FALSE,TRUE,TRUE), capUpLimit=80, capDownLimit=20,
		dataSummaryFrame<-readMultiMatrix(filesInOneLoc, sampleNames[start:end],rowOrder=orderStat$rowOrder, numCluster=numCluster[start],capLimit=capLimit[start:end],	capUpLimit=capUpLimit, capDownLimit=capDownLimit,
				logScale=logScale[start:end], subClusterOrder=orderStat$subClusterOrder, doClustering=F,  bin_size_align=bin_size_align, bin_size=bin_size,scale=scale) 
		##plot sidebar
		if(!is.null(heatmap_anno)){
			plotMultipleSideBar(orderStat$subClusterOrder,orderStat$rowOrder, heatmap_anno=heatmap_anno, heatmap_anno_col=heatmap_anno_col, heatmap_anno_name=heatmap_anno_name)
		}
		
		##plot heatmap
		if(notPlotIndexMatrix){
			for(j in seq(length(fileNumToOrder)+1, end-start+1, by=2)){
				plotHeatmap(dataSummaryFrame$matrix, j, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[j], label=FALSE, adjust_len_proportion=adjust_len_proportion)
				if((j+1) > breaks){
					break;	
				}
				plotHeatmap(dataSummaryFrame$matrix, j+1, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[j+1], label=FALSE, adjust_len_proportion=adjust_len_proportion)
				if((j+1) != breaks){
					par(mar=c(0,0,0,0))
					plot.new()
				}	
			}
		}else{
			
			plotHeatmap(dataSummaryFrame$matrix, 1, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[1], label=FALSE, adjust_len_proportion=adjust_len_proportion,colLim=c(0,100))
			plotHeatmap(dataSummaryFrame$matrix, 2, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[2], label=FALSE, adjust_len_proportion=adjust_len_proportion,colLim=c(0,100))
			par(mar=c(0,0,0,0))
			plot.new()
			for(j in seq(3, end-start+1, by=2)){
				plotHeatmap(dataSummaryFrame$matrix, j, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[j], label=FALSE, adjust_len_proportion=adjust_len_proportion)
				if((j+1) > breaks){
					break;	
				}
				plotHeatmap(dataSummaryFrame$matrix, j+1, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[j+1], label=FALSE, adjust_len_proportion=adjust_len_proportion)
				if((j+1) != breaks){
					par(mar=c(0,0,0,0))
					plot.new()
				}
			}
		}
		#if(!(notPlotIndexMatrix  && (1 %in% fileNumToOrder))){
		#	plotHeatmap(dataSummaryFrame$matrix, 1, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[1], label=FALSE)
		#}
		#if(!(notPlotIndexMatrix  && (2 %in% fileNumToOrder))){
		#	plotHeatmap(dataSummaryFrame$matrix, 2, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[2], label=FALSE)
		#}
		#if(!(notPlotIndexMatrix  && (1 %in% fileNumToOrder) && (2 %in% fileNumToOrder))){
		#	par(mar=c(0,0,0,0))
		#	plot.new()
		#}
		#for(j in seq(3, end-start+1, by=2)){
			
		#	if(!(notPlotIndexMatrix  && (j %in% fileNumToOrder))){
		#		plotHeatmap(dataSummaryFrame$matrix, j, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[j], label=FALSE)
		#	}
		#	if(!(notPlotIndexMatrix  && ((j+1) %in% fileNumToOrder))){
		#		plotHeatmap(dataSummaryFrame$matrix, j+1, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[j+1], label=FALSE)
		#	}
		#	if((j+1) != breaks){
		#		par(mar=c(0,0,0,0))
		#		plot.new()
		#	}		
		#}
		
		
		###############plot metyhylation value
		#	for(i in 1:2){			
		#		plotHeatmap(dataSummaryFrame$matrix, i, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[i], label=FALSE)
		#	}
		###############plot insulater value
		#	par(mar=c(0,0,0,0))
		#	plot.new()	
		###############plot accessbility value
		#	for(i in 3:4){
		#		plotHeatmap(dataSummaryFrame$matrix, i, addNumber=F, xAxisScale=c(-scale,scale), heatMapCols=heatMapCols[i], label=FALSE,colLim=c(0,100))
		#	}
		
		###########plot side bar
		#if(!is.null(orderStat$subClusterOrder)){
		#plotMultipleSideBar(orderStat$subClusterOrder,orderStat$rowOrder, heatmap_anno=heatmap_anno, heatmap_anno_col=heatmap_anno_col, heatmap_anno_name=heatmap_anno_name)
		#}else plot.new()
		
		###########plot tree
		if(!is.null(orderStat$subClusterOrder)){
			plotPhyloTree(as.phylo(orderStat$tree), direction="leftwards")
		}else plot.new()
		
		
		if(addAverage){
			if(!is.null(heatmap_anno)){
				plot.new()
			}
			plotAverage(dataSummaryFrame$average, 1, xAxisScale=c(-scale,scale),xStep=scale/2, yAxisScale=ylim[,1],yStep=(ylim[2,1]-ylim[1,1])/2, label=TRUE, heatmap_autoscale=heatmap_autoscale[1])
			plotAverage(dataSummaryFrame$average, 2, xAxisScale=c(-scale,scale),xStep=scale/2, yAxisScale=ylim[,2],yStep=(ylim[2,2]-ylim[1,2])/2, label=TRUE, heatmap_autoscale=heatmap_autoscale[2])
			for(j in seq(3, end-start+1, by=2)){
				par(mar=c(0,0,0,0))
				plot.new()	
				plotAverage(dataSummaryFrame$average, j, xAxisScale=c(-scale,scale),xStep=scale/2, yAxisScale=ylim[,j],yStep=(ylim[2,j]-ylim[1,j])/2, label=TRUE, heatmap_autoscale=heatmap_autoscale[j])
				if((j+1) > breaks){
					break;	
				}
				plotAverage(dataSummaryFrame$average, j+1, xAxisScale=c(-scale,scale),xStep=scale/2, yAxisScale=ylim[,j+1],yStep=(ylim[2,j+1]-ylim[1,j+1])/2, label=TRUE, heatmap_autoscale=heatmap_autoscale[j+1])
			}
			
			
			###############plot metyhylation value
			#	for(i in 1:2){	
			#		#print(ylim)
			#		plotAverage(dataSummaryFrame$average, i, xAxisScale=c(-scale,scale),xStep=scale/2, yAxisScale=ylim[,i],yStep=(ylim[2,i]-ylim[1,i])/2, label=FALSE)	
			#	}
			###############plot insulater value
			#	par(mar=c(0,0,0,0))
			#	plot.new()	
			###############plot accessbility value
			#	for(i in 3:4){
			#		#print(ylim)
			#		plotAverage(dataSummaryFrame$average, i, xAxisScale=c(-scale,scale),xStep=scale/2, yAxisScale=ylim[,i],yStep=(ylim[2,i]-ylim[1,i])/2, label=FALSE)
			#	}
			##############rest of vacant value
			par(mar=c(0,0,0,0))
			plot.new()	
		}
		
	}
	
	dev.off()
	
}



######function for heatmap

plotHeatmap<-function(x, categoryNum,addNumber=FALSE, xAxisScale=c(-5000,5000), na.color=par("bg"),label=TRUE,heatMapCols="white2red",colLim=NULL, adjust_len_proportion=FALSE,...){
	##make heatmap.col
	heatmap.col<-unlist(strsplit(heatMapCols, "2")[1])
	if(is.na(heatmap.col[2])){
		heatmap.col<-jet.rev.colors
	}
	col<-colorRampPalette(heatmap.col)(100)
	
	q<-x[categoryNum,,]
	if(!is.null(colLim)){
		d<-q
	}else{
		d<-1/(1+exp(-(q-25)/3))
	}
	
	#d<-d/mean(d,na.rm=T)
	breaks<-seq(min(d,na.rm=T), max(d,na.rm=T), length = (length(col)+1))
	if(adjust_len_proportion){
		par(mar=c(0, 1.5, 0,2.0),usr=c(0,1,0,1))
	}else{
		par(mar=c(1.5, 1.5, 1.5,2.0),usr=c(0,1,0,1))
	}
	
	#plot(0,0,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="")
	#rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="") ##draw a rectangle around it..
	#par(new=T,mar=c(0.5, 0.5, 0.5, 1.5),usr=c(0,1,0,1))
	nr=dim(d)[1]
	nc=dim(d)[2]
	#missingExist<-(!invalid(na.color) & any(is.na(d)))
	missingExist <- any(is.na(d))
	if (missingExist) {
		mmat <- ifelse(is.na(d), 1, NA)
		image(1:nc, 1:nr, t(mmat), axes = FALSE, xlab = "", ylab = "", col = na.color)
	}
	if(!is.null(colLim)){
		d[d<colLim[1]]=colLim[1]
		d[d>colLim[2]]=colLim[2]
	}
	image(1:nc, 1:nr, t(d), xlim = 0.5 + c(0, nc), ylim = 0.5 + 
					c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
			breaks = breaks, add = ifelse(missingExist,TRUE,FALSE), ...)
	
	s=FALSE
	if(label){
		s=format(c(as.numeric(xAxisScale[1]),0,as.numeric(xAxisScale[2])),digits=2)
	}
	#s<-ifelse(label,format(c(as.numeric(xAxisScale[1]),0,as.numeric(xAxisScale[2])),digits=2),FALSE)
	
	#axis(1, at= (0.5 + c(0,nc/2,nc)), labels=s, las = 1, tick = TRUE, font=1,cex.axis=1,line=0)
	if(addNumber){
		axis(1, at= (0.5 + c(0,nc/2,nc)), labels=s, las = 1, tick = TRUE, font=1,cex.axis=1,line=0)
		axis(4, 0.5 + seq(0,nr,by=nr/5), labels=format(seq(0,nr,by=nr/5),digits=2), las = 2, tick = TRUE, font=1,cex.axis=1,line=0)
	}
	#box(lty = '1111', col = 'black',cex=1)
}

##TODO: need to make plot multiple sidebar possible...
plotSideBar<-function(subClusterOrder,rowOrder, horizontal=TRUE){
	par(mar=c(1.5,0,1.5,0),usr=c(0,1,0,1))
	subClusterColors<-c("red","blue","orange","black","cyan","purple","palegreen3","yellow","lightpink2","ivory4")
	subClusterColBar=subClusterColors[subClusterOrder]
	names(subClusterColBar)=names(subClusterOrder)
	subClusterColBar=subClusterColBar[rowOrder]
	
	
	
	rsc = subClusterColBar
	rsc.colors = matrix()
	rsc.names = names(table(rsc))
	rsc.i = 1
	for (rsc.name in rsc.names) {
		rsc.colors[rsc.i] = rsc.name
		rsc[rsc == rsc.name] = rsc.i
		rsc.i = rsc.i + 1
	}
	rsc = matrix(as.numeric(rsc), nrow = length(rsc))
	image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
	#image(t(subClusterColBar), col = as.vector(subClusterColBar), axes = FALSE)
	axis(3, 1/2, "", las = 2, tick = FALSE, font=2)
}

plotMultipleSideBar<-function(subClusterOrder,rowOrder, horizontal=TRUE, heatmap_anno=NULL, heatmap_anno_col=c("black"), heatmap_anno_name="test",maxGap=1){
	par(mar=c(1.5,0,1.5,0),usr=c(0,1,0,1))
	subClusterColors<-c("red","blue","orange","black","cyan","purple","palegreen3","yellow","lightpink2","ivory4")
	subClusterColBar=subClusterColors[subClusterOrder]
	names(subClusterColBar)=names(subClusterOrder)
	subClusterColBar=subClusterColBar[rowOrder]
	colBar=NULL
	colBar<-cbind(colBar,subClusterColBar)	
	if(!is.null(heatmap_anno)){
		for(fileOrder in c(1:length(heatmap_anno))){
			rowSideData<-read.table(heatmap_anno[fileOrder],sep="\t",header=F)
			rowSideData_loc<-GRanges(seqnames=rowSideData[,1],ranges=IRanges(rowSideData[,2],rowSideData[,3]),strand="*")
			#print(heatmap_anno[fileOrder])	
			#print(dim(rowSideData))
			data_loc_vector<-strsplit(rowOrder,":")
			data_loc<-GRanges(seqnames=sapply(data_loc_vector,function(x) x[1]),ranges=IRanges(sapply(data_loc_vector,function(x) as.integer(as.numeric(x[2])+as.numeric(x[3]))/2),width=1),strand="*")
			#print(length(rowOrder))
			#print(rowOrder[1])
			#print(dim(data_loc_vector))
			#print(data_loc_vector[1])
			#print(data_loc[1])
			#print(fileOrder)
			#print(heatmap_anno_col[fileOrder])
			rowSideData_bar<-ifelse(countOverlaps(data_loc,rowSideData_loc,maxgap=maxGap)>0,heatmap_anno_col[fileOrder],"white")
			colBar<-cbind(colBar,rowSideData_bar)
		}
		
	}
	
	rsc = colBar
	rsc.colors = matrix()
	rsc.names = names(table(rsc))
	rsc.i = 1
	for (rsc.name in rsc.names) {
		rsc.colors[rsc.i] = rsc.name
		rsc[rsc == rsc.name] = rsc.i
		rsc.i = rsc.i + 1
	}
	print(dim(rsc))
	print(dim(rsc.colors))
	if(is.null(dim(rsc))){
		rsc = matrix(as.numeric(rsc), nrow = length(rsc))
	}else{
		rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
	}
	
	image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
	#axis(3, 1/2, "", las = 2, tick = FALSE, font=2)
}



plotTree<-function(tree, horizontal=TRUE){
	par(mar=c(1.5, 0.5, 1.5, 0.5),usr=c(0,1,0,1))
	plot(tree, horiz = horizontal, axes = FALSE, yaxs = "i", leaflab = "none")
	
}

##tree need to be as.phylo(obj)
plotPhyloTree<-function(tree, direction=c("rightwards")){
	par(mar=c(0, 0, 0, 0))
	#if(is.null(order_dendgram)){
		plot(tree, direction = direction, show.tip.label=FALSE, no.margin=TRUE)
	#}else{
	#	subClusterOrder<-cutree(tree, k = length(order_dendgram))
	#	newOrder<-NULL
	#	for(i in order_dendgram){
	#		newOrder<-c(newOrder,subClusterOrder[subClusterOrder==i])
	#	}
	#	hc<-rotate(tree,newOrder)
	#	plot(hc, direction = direction, show.tip.label=FALSE, no.margin=TRUE)
	#}
	
	
}

plotKey<-function(zlim,col=colorRampPalette(white2red)(100), title="Density"){
	scale01 <- function(x, low = min(x), high = max(x)) {
		x <- (x - low)/(high - low)
		x
	}
	par(mar = c(2.5, 1, 2.5, 1))
	z <- seq(zlim[1], zlim[2], length = length(col))
	breaks<-seq(zlim[1], zlim[2], length = (length(col)+1))
	image(z = matrix(z, ncol = 1), col = col, breaks = breaks, 
			xaxt = "n", yaxt = "n")
	par(usr = c(0, 1, 0, 1))
	lv <- seq(zlim[1], zlim[2], length = 3)
	xv <- scale01(as.numeric(lv), zlim[1], zlim[2])
	axis(1, at = xv, labels = format(lv,digits=2))
	mtext(side = 3, title, line = 0.5,cex=1)
	
	
}

plotAverage<-function(x, categoryNum, xAxisScale=c(-5000,5000), xStep=5000, yAxisScale=c(0,100),yStep=50, label=TRUE, heatmap_autoscale=FALSE,
		subClusterColors=c("red","blue","orange","black","cyan","purple","palegreen3","yellow","lightpink2","ivory4")){
	
	axisSeqForPlot<-seq(xAxisScale[1], xAxisScale[2], by=xStep)	
	d<-x[categoryNum,,]
	if(heatmap_autoscale){
		yAxisScale<-c(min(d,na.rm=T),max(d,na.rm=T))
		yStep=(yAxisScale[2]-yAxisScale[1])/2
	}
	
	if(label){
		par(mar = c(5, 1.5, 2, 2.0),usr=c(0,1,0,1))
	}else{
		par(mar = c(0.5, 1.5, 0, 2.0),usr=c(0,1,0,1))
	}
	if(is.null(dim(d))){
		axisSeq<-seq(1, length(d), by=1)
		plot(axisSeq,d,type="l",xlim = 0.5 + c(min(axisSeq), max(axisSeq)), ylim=c(yAxisScale[1],yAxisScale[2]),xaxt="n",ann=F,yaxt="n",xlab="",ylab="",col=subClusterColors[1],lty=1,font=2,lwd=2)
	}else{
		axisSeq<-seq(1, dim(d)[2], by=1)
		for(i in 1:dim(d)[1]){
			plot(axisSeq,d[i,],type="l",xlim = 0.5 + c(min(axisSeq), max(axisSeq)), ylim=c(yAxisScale[1],yAxisScale[2]),xaxt="n",ann=F,yaxt="n",xlab="",ylab="",col=subClusterColors[i],lty=1,font=2,lwd=2)
			par(new=T)
		}
	}
	s=FALSE
	s2=FALSE
	if(label){
		s=format(axisSeqForPlot,digits=2)
		s2=format(seq(yAxisScale[1],yAxisScale[2],by=yStep),digits=2)
	}
	axis(1,at=seq(axisSeq[1],axisSeq[length(axisSeq)],length.out=length(axisSeqForPlot))+0.5,labels=s,lty=1,font=2,cex.axis=1.0,cex.lab=2.0,font.axis=2,font.lab=2,lwd=1)
	axis(2,at=seq(yAxisScale[1],yAxisScale[2],by=yStep),labels=s2, lty=1,font=2,cex.axis=3.0,cex.lab=3.0,font.axis=2,font.lab=2,lwd=1, las=1)
	#box(lty = '1111', col = 'black',cex=1)
	if(label){
		title(xlab="Distance to center (bp)", ylab="Signal intensity")
	}
	
	abline(v=median(axisSeq)+0.5,lwd=1.5)
	
}

readSingleMatrix<-function(fileName,capLimit=T, capUpLimit=NULL, capDownLimit=NULL, logScale=F, psedoCount=0.1, median=F, bin_size_align=20, bin_size=20, scale=5000){
	content<-read.table(fileName,sep="\t",header=F)	
	content<-content[!duplicated(content),]
	rownames(content)<-paste(content[,1],content[,2],content[,3],content[,4],sep=":")
	dataSeq<-seq(5+((length(content[1,])-5)/2)-as.integer(scale/bin_size_align), 5+((length(content[1,])-5)/2)+as.integer(scale/bin_size_align), by=as.integer(bin_size/bin_size_align))
	value<-NULL	
	for(i in dataSeq){		
		if(floor(i+bin_size/(2*bin_size_align)-1)-ceiling(i-bin_size/(2*bin_size_align)) <= 0){
			value<-cbind(value,as.numeric(content[,ceiling(i-bin_size/(2*bin_size_align))]))	
		}else{
			
			value<-cbind(value,as.numeric(rowMeans(content[,ceiling(i-bin_size/(2*bin_size_align)):floor(i+bin_size/(2*bin_size_align)-1)], na.rm=T)))	
		}
		
	}
	rownames(value)<-rownames(content)	
	if(capLimit){
		if(is.null(capUpLimit)){
			capUpLimit=quantile(unlist(value),probs=seq(0,1,0.05),na.rm =T)[20] #cap top 5%
			
		}
		if(is.null(capDownLimit)){
			capDownLimit=quantile(unlist(value),probs=seq(0,1,0.05),na.rm =T)[2] #cap down 5%
			capDownLimit=ifelse(capDownLimit<0,0,capDownLimit) 
		}
		
		value[value<capDownLimit]=capDownLimit
		value[value>capUpLimit]=capUpLimit
	}
	
	if(logScale){
		value<-log2(value+psedoCount)
	}
	value
}

readMultiMatrix<-function(files, sampleNames, rowOrder=NULL, subClusterOrder=NULL, fileNumToOrder=1, doClustering=T, regionToCluster=c(-500,500),numCluster=2,
		capLimit=NULL, capUpLimit=NULL, capDownLimit=NULL, logScale=NULL, psedoCount=0.1, median=F, bin_size_align=20, bin_size=20, scale=5000){ ##by default, capLimit top 5% and bottom 5% for ChIP-seq, require the first 3 columns in each line have: chr, start, end, strand 
	##do clustering or ordering the index matrix
	binToCluster<-(seq(as.integer((scale+regionToCluster[1])/bin_size_align),as.integer((scale+regionToCluster[2])/bin_size_align),by=as.integer(bin_size/bin_size_align)))/(bin_size/bin_size_align)
	tree<-NULL
	if(!is.null(rowOrder)){
		
	}else{
		indexContent<-readSingleMatrix(files[fileNumToOrder],capLimit=F,bin_size_align=bin_size_align,bin_size=bin_size,scale=scale)
		if(doClustering){
			indexContentToDoClustering<-indexContent[rowSums(is.na(indexContent[,binToCluster]))<=length(binToCluster)/2,]
			rownames(indexContentToDoClustering)<-rownames(indexContent[rowSums(is.na(indexContent[,binToCluster]))<=length(binToCluster)/2,])		
			tree<-hclust(dist(indexContentToDoClustering[,binToCluster]),method="ward")
			rowOrder<-rownames(indexContentToDoClustering[tree$order,])
			subClusterOrder<-cutree(tree, k = numCluster)	
			
		}else{
			if(median){
				#rowOrder<-rownames(indexContent[order(rowMedians(indexContent[,binToCluster], na.rm=T)) , ])
			}else{
				rowOrder<-rownames(indexContent[order(rowMeans(indexContent[,binToCluster], na.rm=T)) , ])
			}
			numCluster=1
		}
	} 
	
	
	nc<-length(seq(0-scale, scale, by=bin_size))
##need to check content row number first, to make sure the matrix is the number
	commonNames<-rowOrder


	matAverage<-array(NA, dim = c(length(sampleNames),numCluster, nc))
	for(i in 1:length(files)){
		content<-readSingleMatrix(files[i],capLimit=ifelse(is.null(capLimit),FALSE,capLimit[i]),capUpLimit=capUpLimit[i], capDownLimit=capDownLimit[i], logScale=ifelse(is.null(logScale),FALSE,logScale[i]), psedoCount=psedoCount, median=median,bin_size_align=bin_size_align,bin_size=bin_size,scale=scale)
		commonNames<-intersect(commonNames,rownames(content))

		if(median){
			#matAverage[i,]<-colMedians(content,na.rm=T)
		}else{
			if(is.null(subClusterOrder) || numCluster<2){
				matAverage[i,1,]<-colMeans(content,na.rm=T)
			}else{
				for(num in c(1:numCluster)){
					cluster_name<-intersect(names(subClusterOrder[subClusterOrder==num]),rownames(content))
					matAverage[i,num,]<-colMeans(content[cluster_name,],na.rm=T)
				}
			}
			
		}
		
	}
	nr<-length(commonNames)
	mat<-array(NA, dim = c(length(sampleNames),nr,nc))
	for(i in 1:length(files)){
		content<-readSingleMatrix(files[i],capLimit=ifelse(is.null(capLimit),FALSE,capLimit[i]),capUpLimit=capUpLimit[i], capDownLimit=capDownLimit[i], logScale=ifelse(is.null(logScale),FALSE,logScale[i]), psedoCount=psedoCount, median=median,bin_size_align=bin_size_align,bin_size=bin_size,scale=scale)
		content<-content[commonNames,]
		mat[i,,]=content
	}
	
	
	min.dat<-ifelse(is.null(capDownLimit),min(mat, na.rm=T), min(capDownLimit,na.rm=T))
	max.dat<-ifelse(is.null(capUpLimit),max(mat, na.rm=T), max(capUpLimit,na.rm=T))
	zlim<-c(min.dat,max.dat)
	dataFrame<-list(name=sampleNames,average=matAverage, matrix=mat,tree=tree, limit=zlim)
	dataFrame
}



clusteringSingleMatrix<-function(fileName, regionToCluster=c(-500,500),numCluster=2,bin_size_align=1, bin_size=20, scale=1000){ ##object read from readSingleMatrix()
	binToCluster<-(seq(as.integer((scale+regionToCluster[1])/bin_size_align),as.integer((scale+regionToCluster[2])/bin_size_align),by=as.integer(bin_size/bin_size_align)))/(bin_size/bin_size_align)
	content<-readSingleMatrix(fileName,capLimit=F,bin_size_align=bin_size_align,bin_size=bin_size,scale=scale)
	indexContentToDoClustering<-content[rowSums(is.na(content[,binToCluster]))<=length(binToCluster)/2,]
	rownames(indexContentToDoClustering)<-rownames(content[rowSums(is.na(content[,binToCluster]))<=length(binToCluster)/2,])		
	tree<-hclust(dist(indexContentToDoClustering[,binToCluster]),method="ward")
	rowOrder<-rownames(indexContentToDoClustering[tree$order,])
	subClusterOrder<-cutree(tree, k = numCluster)
	
	nr<-length(rowOrder)
	nc<-length(seq(0-scale, scale, by=bin_size))
	mat<-array(NA, dim = c(1,nr,nc))
	matAverage<-array(NA, dim = c(1,numCluster, nc))
	mat[1,,]=indexContentToDoClustering[rowOrder,]
	for(num in c(1:numCluster)){
		cluster_name<-names(subClusterOrder[subClusterOrder==num])
		matAverage[1,num,]<-colMeans(content[cluster_name,],na.rm=T)
	}
	
	zlim<-c(min(mat, na.rm=T),max(mat, na.rm=T))
	dataFram<-list(rowOrder=rowOrder, subClusterOrder=subClusterOrder,average=matAverage, tree=tree, matrix=mat, limit=zlim)
	dataFram
	#dataList<-list(dataFram=,as.dendrogram(tree))
} 

##special version for old generateHeatmapForNOMeseq.R
clusteringMultiMatrix<-function(fileNames, regionToCluster=c(-500,500),numCluster=2,bin_size_align=1, bin_size=20, scale=1000){ ##object read from readSingleMatrix()
	binToCluster<-(seq(as.integer((scale+regionToCluster[1])/bin_size_align),as.integer((scale+regionToCluster[2])/bin_size_align),by=as.integer(bin_size/bin_size_align)))/(bin_size/bin_size_align)
	print(binToCluster)
	indexContentToDoClustering<-NULL
	for(fileName in fileNames){
		content<-readSingleMatrix(fileName,capLimit=F,bin_size_align=bin_size_align,bin_size=bin_size,scale=scale)
		print(dim(content))
		content<-content[rowSums(is.na(content[,binToCluster]))<=length(binToCluster)/2,]
		if(is.null(indexContentToDoClustering)){
			indexContentToDoClustering<-content[,binToCluster]
			rownames(indexContentToDoClustering)<-rownames(content)
		}else{
			common<-intersect(rownames(indexContentToDoClustering),rownames(content))
			indexContentToDoClustering<-indexContentToDoClustering[common,]
			indexContentToDoClustering<-cbind(indexContentToDoClustering,content[common,][,binToCluster])
		}
		
	}

			
	tree<-hclust(dist(indexContentToDoClustering),method="ward")
	rowOrder<-rownames(indexContentToDoClustering[tree$order,])
	subClusterOrder<-cutree(tree, k = numCluster)
	
	dataFram<-list(rowNames=rowOrder, subClusterOrder=subClusterOrder, capUpLimit=NULL, capDownLimit=NULL)
	
	dataFram
	#dataList<-list(dataFram=,as.dendrogram(tree))
} 


clusteringMultiMatrixNew<-function(fileNames, regionToCluster=c(-500,500),numCluster=2,bin_size_align=1, bin_size=20, scale=1000, order_dendgram=NULL){ ##object read from readSingleMatrix()
	binToCluster<-(seq(as.integer((scale+regionToCluster[1])/bin_size_align),as.integer((scale+regionToCluster[2])/bin_size_align),by=as.integer(bin_size/bin_size_align)))/(bin_size/bin_size_align)
	print(binToCluster)
	indexContentToDoClustering<-NULL
	for(fileName in fileNames){
		content<-readSingleMatrix(fileName,capLimit=F,bin_size_align=bin_size_align,bin_size=bin_size,scale=scale)
		print(dim(content))
		content<-content[rowSums(is.na(content[,binToCluster]))<=length(binToCluster)/2,]
		if(is.null(indexContentToDoClustering)){
			indexContentToDoClustering<-content[,binToCluster]
			rownames(indexContentToDoClustering)<-rownames(content)
		}else{
			common<-intersect(rownames(indexContentToDoClustering),rownames(content))
			indexContentToDoClustering<-indexContentToDoClustering[common,]
			indexContentToDoClustering<-cbind(indexContentToDoClustering,content[common,][,binToCluster])
		}
		
	}
	
	
	tree<-hclust(dist(indexContentToDoClustering),method="ward")
	
	subClusterOrder<-cutree(tree, k = numCluster)
	
	if(is.null(order_dendgram)){
		rowOrder<-rownames(indexContentToDoClustering[tree$order,])
		dataFram<-list(rowOrder=rowOrder, subClusterOrder=subClusterOrder,average=NULL, tree=tree, matrix=NULL, limit=NULL)
	}else{ ##reordering dendgram by provided order of clusters
		if(length(order_dendgram) != numCluster){
			stop("length of provided clusters' order is different from required cluster number")
		}
		newOrder<-NULL
		for(i in order_dendgram){
			newOrder<-c(newOrder,subClusterOrder[subClusterOrder==i])
		}
		hc<-rotate(tree,newOrder)
		rowOrder<-rownames(indexContentToDoClustering[hc$order,])
		subClusterOrder<-cutree(hc, k = numCluster)
		dataFram<-list(rowOrder=rowOrder, subClusterOrder=subClusterOrder,average=NULL, tree=hc, matrix=NULL, limit=NULL)
	}

	dataFram

} 


file_proportion<-function(files,breaks=4,bin_size_align=1, bin_size=20, scale=1000){	
	matrix_len<-NULL
	for(i in seq(1,length(files),by=breaks)){
		start<-i
		dataSummaryFrame<-readSingleMatrix(files[start],capLimit=F,bin_size_align=bin_size_align,bin_size=bin_size,scale=scale)
		matrix_len<-c(matrix_len,dim(dataSummaryFrame)[1])
	}	
	numEqualFracInRows=6*sum(matrix_len/max(matrix_len,na.rm=T),na.rm=T) #6:2 or when no average plot, 6
#	numRows=as.integer(length(files)/breaks)
#	heatmapFrac=6/numEqualFracInRows

	lht<-6*matrix_len/max(matrix_len,na.rm=T)
	propor<-list(lht=lht,numEqualFracInRows=numEqualFracInRows)
	propor
}



###added 2015-11-20
plotHeatmapGeneral<-function(d, addNumber=FALSE,  na.color=par("bg"),label=TRUE,heatMapCols=jet.colors,colLim=NULL, ...){	
	heatMapColsPlot<-colorRampPalette(heatMapCols)(100)
	
	breaks<-seq(min(d,na.rm=T), max(d,na.rm=T), length = (length(heatMapColsPlot)+1))
	par(mar=c(1.5, 1.5, 1.5,2.0),usr=c(0,1,0,1))
	nr=dim(d)[1]
	nc=dim(d)[2]	
	missingExist <- any(is.na(d))
	if (missingExist) {
		mmat <- ifelse(is.na(d), 1, NA)
		image(1:nc, 1:nr, t(mmat), axes = FALSE, xlab = "", ylab = "", col = na.color)
	}
	if(!is.null(colLim)){
		d[d<colLim[1]]=colLim[1]
		d[d>colLim[2]]=colLim[2]
		breaks<-seq(colLim[1], colLim[2], length = (length(heatMapColsPlot)+1))
	}
	image(1:nc, 1:nr, t(d), xlim = 0.5 + c(0, nc), ylim = 0.5 + 
					c(0, nr), axes = FALSE, xlab = "", ylab = "", col = heatMapColsPlot, 
			breaks = breaks, add = ifelse(missingExist,TRUE,FALSE), ...)
	
	if(addNumber){
		axis(4, 0.5 + seq(0,nr,by=nr/5), labels=format(seq(0,nr,by=nr/5),digits=2), las = 2, tick = TRUE, font=1,cex.axis=0.5,line=0)
	}
	box(lty = '1111', col = 'black',cex=1)
}


plotSideBarGeneral<-function(col_list.reorder, horizontal=TRUE){
	if(horizontal){
		par(mar=c(1.5,0,1.5,0),usr=c(0,1,0,1))
	}else{
		par(mar=c(0,1.5,0,2.0),usr=c(0,1,0,1))
	}
	
	if(horizontal){
		image(1, 1:length(col_list.reorder), t(as.matrix(1:length(col_list.reorder))), col = as.vector(col_list.reorder), axes = FALSE)
	}else{
		image(1:length(col_list.reorder), 1, as.matrix(1:length(col_list.reorder)), col = as.vector(col_list.reorder), axes = FALSE)
	}
	axis(3, 1/2, "", las = 2, tick = FALSE, font=2)
}

plotTreeGeneral<-function(tree, horizontal=TRUE){
	if(horizontal){
		par(mar=c(1.5,0.5,1.5,0.5),usr=c(0,1,0,1))
	}else{
		par(mar=c(0.5, 0, 1.5, 0.5),usr=c(0,1,0,1))
	}
	plot(tree, horiz = horizontal, axes = FALSE, yaxs = "i", leaflab = "none")
	
}


plotTree<-function(tree, horizontal=TRUE){
	par(mar=c(0.5, 0, 1.5, 0.5),usr=c(0,1,0,1))
	plot(tree, horiz = horizontal, axes = FALSE, yaxs = "i", leaflab = "none")
	
}








