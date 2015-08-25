for (e in commandArgs(TRUE)) {
	ta = strsplit(e,"=",fixed=TRUE)
	if(! is.na(ta[[1]][2])) {
		if(ta[[1]][1] == "wd"){
			wd<-ta[[1]][2]  ## directory to make plot file
		}
		if(ta[[1]][1] == "rpkmFile"){
			rpkmFile<-ta[[1]][2]
		}
		if(ta[[1]][1] == "p"){
			p<-ta[[1]][2]
		}
		if(ta[[1]][1] == "geneInfoFile"){
			geneInfoFile<-ta[[1]][2]
		}

	}
}


geneBodyMethyRpkm<-function(x="57epigenomes.RPKM.pc.txt",p,geneInfoFile="Ensembl_v65.Gencode_v10.ENSG.gene_info"){
	geneInfo<-read.delim(geneInfoFile,sep="\t",header=F)
	geneInfo<-geneInfo[!duplicated(geneInfo[,1]),]
	geneInfo<-geneInfo[!is.na(geneInfo[,1]),]
	geneInfo[,5]<-sub("-1","-",geneInfo[,5])
	geneInfo[,5]<-sub("1","+",geneInfo[,5])
	rownames(geneInfo)<-geneInfo[,1]
	
	rpkm<-read.table(x,sep="\t",header=T)
	#rownames(rpkm)<-rpkm[,1]
	colnames(rpkm)<-colnames(rpkm)[2:length(rpkm[1,])]
	rpkm<-rpkm[,1:(length(rpkm[1,])-1)]
	
	common<-intersect(rownames(geneInfo),rownames(rpkm))
	for(i in 2:length(colnames(rpkm))){
		sampleName<-colnames(rpkm)[i]
		print(sampleName)
		if(!file.exists(paste("/home/unix/lyping/compbio/project/meQTL/prelim/methy_enhancer/methy_signal/",sampleName,"_SBS.hg19.bw",sep=""))){
			next
		}
		
		sample.rpkm<-rpkm[common,][,i]
		sample.rpkm.cor<-cbind(geneInfo[common,][,c(2:4,1)],sample.rpkm,geneInfo[common,][,5])
		sample.rpkm.cor[,1]<-paste("chr",sample.rpkm.cor[,1],sep="")
		sample.rpkm.cor.sort<-sample.rpkm.cor[order(sample.rpkm.cor[,5]),]
		

		sample.perc0_20.rpkm<-sample.rpkm.cor.sort[1:(length(sample.rpkm.cor.sort[,1])/5),]
		sample.perc20_40.rpkm<-sample.rpkm.cor.sort[(length(sample.rpkm.cor.sort[,1])/5+1):(length(sample.rpkm.cor.sort[,1])*2/5),]
		sample.perc40_60.rpkm<-sample.rpkm.cor.sort[(length(sample.rpkm.cor.sort[,1])*2/5+1):(length(sample.rpkm.cor.sort[,1])*3/5),]
		sample.perc60_80.rpkm<-sample.rpkm.cor.sort[(length(sample.rpkm.cor.sort[,1])*3/5+1):(length(sample.rpkm.cor.sort[,1])*4/5),]
		sample.perc80_100.rpkm<-sample.rpkm.cor.sort[(length(sample.rpkm.cor.sort[,1])*4/5+1):(length(sample.rpkm.cor.sort[,1])),]
		write.table(sample.perc0_20.rpkm,paste(sampleName,".",p,".rpkm.cor.percent0-20.bed",sep=""),sep="\t",quote=F,row.names =F,col.names =F)
		write.table(sample.perc20_40.rpkm,paste(sampleName,".",p,".rpkm.cor.percent20-40.bed",sep=""),sep="\t",quote=F,row.names =F,col.names =F)
		write.table(sample.perc40_60.rpkm,paste(sampleName,".",p,".rpkm.cor.percent40-60.bed",sep=""),sep="\t",quote=F,row.names =F,col.names =F)
		write.table(sample.perc60_80.rpkm,paste(sampleName,".",p,".rpkm.cor.percent60-80.bed",sep=""),sep="\t",quote=F,row.names =F,col.names =F)
		write.table(sample.perc80_100.rpkm,paste(sampleName,".",p,".rpkm.cor.percent80-100.bed",sep=""),sep="\t",quote=F,row.names =F,col.names =F)
		getMethy(paste(sampleName,".",p,".rpkm.cor.percent0-20.bed",sep=""),paste(sampleName,".",p,".rpkm.cor.percent0-20.",sampleName,"-SBS",sep=""),sampleName)
		getMethy(paste(sampleName,".",p,".rpkm.cor.percent20-40.bed",sep=""),paste(sampleName,".",p,".rpkm.cor.percent20-40.",sampleName,"-SBS",sep=""),sampleName)
		getMethy(paste(sampleName,".",p,".rpkm.cor.percent40-60.bed",sep=""),paste(sampleName,".",p,".rpkm.cor.percent40-60.",sampleName,"-SBS",sep=""),sampleName)
		getMethy(paste(sampleName,".",p,".rpkm.cor.percent60-80.bed",sep=""),paste(sampleName,".",p,".rpkm.cor.percent60-80.",sampleName,"-SBS",sep=""),sampleName)
		getMethy(paste(sampleName,".",p,".rpkm.cor.percent80-100.bed",sep=""),paste(sampleName,".",p,".rpkm.cor.percent80-100.",sampleName,"-SBS",sep=""),sampleName)
		
		plotMultipleWig(sampleName,p,files=c(paste(sampleName,".",p,".rpkm.cor.percent0-20.",sampleName,"-SBS.",sampleName,"_SBS.alignedTo.",sampleName,".6200.1.txt",sep=""),
											paste(sampleName,".",p,".rpkm.cor.percent20-40.",sampleName,"-SBS.",sampleName,"_SBS.alignedTo.",sampleName,".6200.1.txt",sep=""),
											paste(sampleName,".",p,".rpkm.cor.percent40-60.",sampleName,"-SBS.",sampleName,"_SBS.alignedTo.",sampleName,".6200.1.txt",sep=""),
											paste(sampleName,".",p,".rpkm.cor.percent60-80.",sampleName,"-SBS.",sampleName,"_SBS.alignedTo.",sampleName,".6200.1.txt",sep=""),
											paste(sampleName,".",p,".rpkm.cor.percent80-100.",sampleName,"-SBS.",sampleName,"_SBS.alignedTo.",sampleName,".6200.1.txt",sep="")))

	}
}

getMethy<-function(x, prefix,sampleName){
	cmd<-paste("perl /home/unix/lyping/compbio/code/mytools/perl/alignWigToBed.pl --average --alignment_mode 6 --prefix ",prefix," --locs ", x, " /home/unix/lyping/compbio/project/meQTL/prelim/methy_enhancer/methy_signal/",sampleName,"_SBS.hg19.bw --omit_plot_step --data_matrix_scale 6200 --plot_x_axis_scale 3100 --plot_x_axis_step 3000 --colors black --lty 1  --lengends methy --plot_y_axis_max 1.0 --plot_y_axis_step 0.25",sep="")
	system(cmd,intern=T,ignore.stderr=T)
	
}

plotMultipleWig<-function(sampleName, p, files,colors=c("black","blue","green","orange","red")){
	pdf(paste(sampleName,".",p,".gene_body_methy_vs_expr.pdf",sep=""), paper="special", height=5, width=5)
	par(oma=c(0, 0, 0, 0),mar=c(3, 3, 3, 1))
	scale=6100
	axistep=3000
	axisSeqForPlot<-c(-6100,-3100,3000,5900)
	bin_size_align=1
	bin_size=20
	for(j in 1:length(files)){
		gch<-read.table(files[j],sep="\t",header=F)	
		axisSeq<-seq(0-scale, scale, by=1)
		dataSeq<-seq(5+((length(gch[1,])-5)/2)-scale, 5+((length(gch[1,])-5)/2)-scale+scale * 2, by=1)
		valueGch<-array()
		for(i in dataSeq){
			if(floor(i+bin_size/(2*bin_size_align)-1)-ceiling(i-bin_size/(2*bin_size_align)) <= 0){
				valueGch<-cbind(valueGch,mean(gch[,ceiling(i-bin_size/(2*bin_size_align))], na.rm=T))	
			}else{
				valueGch<-cbind(valueGch,mean(colMeans(gch[,ceiling(i-bin_size/(2*bin_size_align)):floor(i+bin_size/(2*bin_size_align)-1)], na.rm=T), na.rm=T))	
			}
			
		}
		valueGch<-valueGch[,2:length(valueGch[1,])]	
		plot(axisSeq,valueGch,type="l",axes=FALSE,xlab="",ylab="",ylim=c(0,1),col=colors[j],lty=1,font=2,lwd=2.5)
		par(new=T)
	}
	axis(2,at=seq(0,1,by=0.25),lty=1,font=2,cex.axis=1.0,cex.lab=1.0,font.lab=2,lwd=2)
	axis(1,at=axisSeqForPlot,labels=c("-3000","TSS","TTS","+3000"),lty=1,font=2,cex.axis=1.2,cex.lab=1.2,font.lab=2,lwd=2)


	title(sampleName, cex.main = 1.2, font.main= 4, col.main= "black",xlab="Distance(bp)",ylab="DNA methylation")
	legend("topright",c("0-20% expr","20-40% expr","40-60% expr","60-80% expr","80-100% expr"), col=colors,lty=1,cex=0.5,lwd=2)
	abline(v=axisSeqForPlot[2],lty=3,lwd=1.5)
	abline(v=axisSeqForPlot[3],lty=3,lwd=1.5)
	dev.off()
}


setwd(wd)

geneBodyMethyRpkm(x=rpkmFile,p=p, geneInfoFile=geneInfoFile)
