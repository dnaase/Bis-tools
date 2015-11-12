for (e in commandArgs(TRUE)) {
	ta = strsplit(e,"=",fixed=TRUE)
	if(! is.na(ta[[1]][2])) {
		if(ta[[1]][1] == "wd"){
			wd<-ta[[1]][2]  ## directory to make plot file
		}
		if(ta[[1]][1] == "markFile"){
			markFile<-ta[[1]][2]  
		}
		if(ta[[1]][1] == "matrix"){
			matrix<-ta[[1]][2]  
		}
		if(ta[[1]][1] == "figure_file"){
			figure_file<-ta[[1]][2]  
		}
		if(ta[[1]][1] == "figure_mark_file"){
			figure_mark_file<-ta[[1]][2]  
		}
		if(ta[[1]][1] == "bin_size"){
			bin_size<-as.numeric(ta[[1]][2])
		}
	}
}

setwd(wd)

d<-read.table(matrix,sep="\t",header=F)

if((dim(d[rowSums(is.na(d))==d[1,],])[1])>=length(d[,1])){ ##all are NA in the matrix
	stop("all value are NA");
}

d<-log10(d+1)

pdf(figure_file, paper="special", height=10, width=10)
	col<-colorRampPalette(c("green","yellow","red"))(100)
	#col<-colorRampPalette(c("blue","yellow"))(100)
	na.color=par("bg")

	breaks<-seq(min(d,na.rm=T), max(d,na.rm=T), length = (length(col)+1))
	par(mar=c(1.5, 1.5, 1.5,2.0),usr=c(0,1,0,1))
	nr=dim(d)[1]
	nc=dim(d)[2]
	missingExist <- any(is.na(d))
	if (missingExist) {
		mmat <- ifelse(is.na(d), 1, NA)
		image(1:nc, 1:nr, t(mmat), axes = FALSE, xlab = "", ylab = "", col = na.color)
	}
	image(1:nc, 1:nr, t(d), xlim = 0.5 + c(0, nc), ylim = 0.5 + 
					c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
			breaks = breaks, add = ifelse(missingExist,TRUE,FALSE))
	
dev.off()

##plot figure with mark

mark_matrix<-read.table(markFile,sep="\t",header=F)

pdf(figure_mark_file, paper="special", height=10, width=10)
	#col<-colorRampPalette(c("green","yellow","red"))(100)
	col<-colorRampPalette(c("white","black"))(100)
	na.color=par("bg")
	#mark.color<-colorRampPalette(c("white","blue"))(2)
	mark.color<-"blue"
	breaks<-seq(min(d,na.rm=T), max(d,na.rm=T), length = (length(col)+1))
	par(mar=c(1.5, 1.5, 1.5,2.0),usr=c(0,1,0,1))
	nr=dim(d)[1]
	nc=dim(d)[2]
	missingExist <- any(is.na(d))
	if (missingExist) {
		mmat <- ifelse(is.na(d), 1, NA)
		image(1:nc, 1:nr, t(mmat), axes = FALSE, xlab = "", ylab = "", col = na.color)
	}
	image(1:nc, 1:nr, t(d), xlim = 0.5 + c(0, nc), ylim = 0.5 + 
					c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
			breaks = breaks, add = ifelse(missingExist,TRUE,FALSE))
	
	par(new=T,mar=c(1.5, 1.5, 1.5,2.0),usr=c(0,1,0,1))
	#image(1:nc, 1:nr, t(mark_matrix), axes = FALSE, xlab = "", ylab = "", col = mark.color,breaks = seq(0,2, length = 3) )
	image(1:nc, 1:nr, t(mark_matrix), axes = FALSE, xlab = "", ylab = "", col = mark.color)

dev.off()

