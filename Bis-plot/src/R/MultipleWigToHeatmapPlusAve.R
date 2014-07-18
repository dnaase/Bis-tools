# MultipleWigToHeatmapPlusAve.R
# 
# Author: yaping
# Sep 25, 2013
# 11:53:05 AM
## 1. get multiple matrix, clustering 1st matrix(hclust/k-means) or order by mean value in 1st matrix
## 2. layout canvas to transfer multiple matrix into multiple heatmap
## 3. layout average plot under each of matrix
## 4. label multiple matrix name, also cluster tree

###############################################################################
source("/home/uec-00/yapingli/code/mytools/R/heatmapRUtils.R")
library(gplots)
library(colorRamps)
library(RColorBrewer)
library(GenomicRanges)
##hclust in this package is much faster, and allow more rows...
library(fastcluster)

#default value
jet.colors <-c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
jet.rev.colors <-c("#7F0000","red","#FF7F00","yellow","#7FFF7F","cyan","#007FFF","blue","#00007F")
white2red<-c("white","red")

files<-NULL
categoryNames<-NULL
sampleNames<-NULL
heatMapCols<-NULL
ylabForAvePlot<-NULL

capLimit<-NULL
capUpLimit<-NULL
capDownLimit<-NULL
limit<-c(NULL,NULL)
logScale<-NULL
numCluster<-NULL

scale=5000
bin_size_align=20
bin_size=20
heatmap_clustering_scale = 1000
heatmap_clustering_bin_size = 20
heatmap_clustering_bin_size_align = 1

fileNumToOrder=1
regionToClusterLow=-500
regionToClusterHigh=500
addAverage=FALSE


breaks<-NULL 
wd="./"
prefix<-NULL
for (e in commandArgs(TRUE)) {
	ta = strsplit(e,"=",fixed=TRUE)
	if(! is.na(ta[[1]][2])) {
		if(ta[[1]][1] == "wd"){
			wd<-ta[[1]][2]
		}
		if(ta[[1]][1] == "inputFn"){
			files<-c(files,ta[[1]][2])
		}
		if(ta[[1]][1] == "prefix"){
			prefix<-c(prefix,ta[[1]][2])
		}
		if(ta[[1]][1] == "clusterNum"){
			numCluster<-c(numCluster,as.numeric(ta[[1]][2]))
		}
		if(ta[[1]][1] == "categoryNames"){
			categoryNames<-c(categoryNames,ta[[1]][2])
		}
		if(ta[[1]][1] == "sampleNames"){
			sampleNames<-c(sampleNames,ta[[1]][2])
		}
		if(ta[[1]][1] == "heatMapCols"){
			heatMapCols<-c(heatMapCols,ta[[1]][2])
		}

		if(ta[[1]][1] == "scale"){
			scale<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "bin_size_align"){
			bin_size_align<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "move_step"){
			bin_size<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "logScale"){
			logScale<-c(logScale, as.logical(ta[[1]][2]))
		}
		if(ta[[1]][1] == "capLimit"){
			capLimit<-c(capLimit, as.logical(ta[[1]][2]))
		}
		if(ta[[1]][1] == "capUpLimit"){
			capUpLimit<- as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "capDownLimit"){
			capDownLimit<- as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "fileNumToPrintSideBar"){
			fileNumToOrder<-as.numeric(ta[[1]][2])
		}

		if(ta[[1]][1] == "addAverage"){
			addAverage<-as.logical(ta[[1]][2])
		}
		if(ta[[1]][1] == "regionToClusterLow"){
			regionToClusterLow<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "regionToClusterHigh"){
			regionToClusterHigh<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "heatmap_clustering_scale"){
			heatmap_clustering_scale<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "heatmap_clustering_bin_size"){
			heatmap_clustering_bin_size<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "heatmap_clustering_bin_size_align"){
			heatmap_clustering_bin_size_align<-as.numeric(ta[[1]][2])
		}
		if(ta[[1]][1] == "breaks"){
			breaks<-as.numeric(ta[[1]][2])
		}
		
	}
}


setwd(wd)

regionToCluster<-c(regionToClusterLow,regionToClusterHigh)
if(!is.null(capUpLimit)){
	limit<-c(capDownLimit,capUpLimit)
}

for(i in seq(1,length(files),by=breaks)){ ##do it in multiple locations
	start<-i
	end<-i+breaks-1
	filesInOneLoc<-files[start:end]
	limit<-plotMultiWigToHeatmapPlusAveInSingleLoc(filesInOneLoc, sampleNames[start:end],prefix=prefix[start], limit=limit, fileNumToOrder=fileNumToOrder, regionToCluster=regionToCluster,numCluster=numCluster[start], 
			bin_size_align=bin_size_align,bin_size=bin_size, scale=scale,heatmap_clustering_bin_size_align=heatmap_clustering_bin_size_align, heatmap_clustering_bin_size=heatmap_clustering_bin_size,heatmap_clustering_scale=heatmap_clustering_scale,
			capLimit=capLimit[start:end],logScale=logScale[start:end],heatMapCols=heatMapCols[start:end], addAverage=addAverage)

}




