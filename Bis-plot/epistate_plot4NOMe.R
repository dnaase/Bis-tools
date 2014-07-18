#########################################
#	Epistate plot for Nome-Seq			#
#	Huy Q. Dinh & Ben Berman			#
#	USC Epigenome Center				#
#	July, 2014:huy.q.dinh@gmail.com		#
#########################################

draw_NoMe_read <- function(meths, GCH_idx, HCG_idx, X, Y, label, color, radius=0.5, is_nome = 1) {
	require(plotrix)
	returned_coors <- NULL

	row_space = 2.5
	range = seq(0, 1, 0.05)
	color = "black"
	cols = rev(colorRampPalette(c(color, 'white'))(length(range)))
	for (i in 1:length(HCG_idx)) {
		if (meths[HCG_idx[i]]>=0) {
			rect(X + i*row_space, Y-1.5, X + i*row_space + 2, Y-1.5 + 2.5, col = cols[which(range>=meths[HCG_idx[i]])[1]], border = color)
			returned_coors <- rbind(returned_coors, c(X + i*row_space, Y-1.5))
		}
	}
	X = X + length(HCG_idx)*row_space + 5

	if (is_nome == 1) {
		row_space = 2
		range = seq(0, 1, 0.05)
		color = "darkgreen"
		cols = rev(colorRampPalette(c(color, 'white'))(length(range)))
	
		for (i in 1:length(GCH_idx)) {
			if (meths[GCH_idx[i]]>=0) {
				draw.circle(X + i*row_space+2.5, Y, radius, col = cols[which(range>=meths[GCH_idx[i]])-1], border = color)
				returned_coors <- rbind(returned_coors, c(X + i*row_space+2.5, Y))
			}
		}
	} 
	
	returned_coors	
}

get_coors <- function(chr_name, start, end, GCHs, HCGs) {
	GCH <- vector("list", length(sample_list))
	for (i in 1:length(GCHs))	 {
		temp <- GCHs[[i]]
		temp <- subset(temp, temp[,1] == chr_name)
		temp <- subset(temp, temp[,2] %in% start:end)
		GCH[[i]] <- temp
	}
	HCG <- vector("list", length(sample_list))
	for (i in 1:length(HCGs))	 {
		temp <- HCGs[[i]]
		temp <- subset(temp, temp[,1] == chr_name)
		temp <- subset(temp, temp[,2] %in% start:end)
		HCG[[i]] <- temp
	}
	list(GCH, HCG)
}

get_reads <- function(chr_name, coors, GCHs, HCGs) {
	meth_strs <- vector("list", length(GCHs))
	for (i in 1:length(GCHs)) {
		temp2 <- get_coors(chr_name, coors[1], coors[length(coors)], GCHs, HCGs)
		temp <- rbind(temp2[[1]][[i]], temp2[[2]][[i]])
		reads <- unique(temp[,6])
		strs <- NULL
		for (j in 1:length(reads)) {
			temp2 <- subset(temp, temp[,6] == reads[j])
			meth_str <- rep("-", length(coors))
			meth_str[which(coors %in% temp2[,2])] = as.character(temp2[,3])
			strs <- rbind(strs, meth_str)
		}
		meth_strs[[i]] = strs
	}
	names(meth_strs) = names(GCHs)
	meth_strs
}

get_dist <- function(t2) {
	m <- matrix(0, nrow = nrow(t2), ncol = nrow(t2))
	for (i in 1:(nrow(m)-1)) {
		for (j in (i+1):nrow(m)) {
			m[i,j] = sum(t2[i,]!=t2[j,])
			m[j,i] = m[i,j]
		}
	}
	m
}

BED_file = "test.bed"
GCH_file = "test.gch.reads.txt"
HCG_file = "test.hcg.reads.txt"
wd="./"
Args <- commandArgs()
for (i in 1:length(Args)) {
	if (Args[i] == "--wd") wd = Args[i+1]
	if (Args[i] == "--f") BED_file = Args[i+1]
	if (Args[i] == "--l") sample_file = Args[i+1]
}
setwd(wd)
regs <- read.table(BED_file)
sample_list <- read.table(sample_file)[,1]

GCHs <- vector("list", length(sample_list))
HCGs <- vector("list", length(sample_list))
for (i in 1:length(sample_list)) {
	GCH_file = paste(sample_list[i], "gch.reads.txt", sep = ".")
	HCG_file = paste(sample_list[i], "hcg.reads.txt", sep = ".")
	
	GCH <- read.table(GCH_file)

	#minus_ind = which(GCH[,5]=="-")
	#GCH[minus_ind,2] = GCH[minus_ind,2]-1
	HCG <- read.table(HCG_file)
	#minus_ind = which(HCG[,5]=="-")
	#HCG[minus_ind,2] = HCG[minus_ind,2]-1
	GCHs[[i]] = GCH
	HCGs[[i]] = HCG
}
names(GCHs) = names(HCGs) = sample_list
print(sample_list)
for (i in 1:nrow(regs)) {
	PDF_filename = paste(paste(regs[i,1], paste(regs[i,2], regs[i,3], sep = "-"), sep = ":"), "pdf", sep = ".")
	chr_name = regs[i,1]
	temp <- get_coors(chr_name, regs[i,2], regs[i,3], GCHs, HCGs)
	GCH_tmp <- temp[[1]]
	HCG_tmp <- temp[[2]]
	GCH_coors <- c()
	HCG_coors <- c()
	for (i in 1:length(GCH_tmp)) {
		GCH_coors <- c(GCH_coors, GCH_tmp[[i]][,2])
		HCG_coors <- c(HCG_coors, HCG_tmp[[i]][,2])
	}
	coors <- sort(unique(c(GCH_coors, HCG_coors)))
	reads <- get_reads(chr_name, coors, GCHs, HCGs)

	size = length(coors)*200/40
	
	canvas.X = c(-size,size)
	canvas.Y = c(-size,size)
	
	HCG_idx = which(coors %in% HCG_coors)
	GCH_idx = setdiff(1:length(coors), HCG_idx)

	CpG_window = length(coors)

	pdf(PDF_filename , paper = "special", width=20, height=13)

	plot(0,0,type = "n", xlim = canvas.X, canvas.Y, xlab='', ylab='', xaxt='n',yaxt='n')
	lineX=canvas.X[2]/1.5
	lineY=canvas.Y[2]
	segments(-lineX,lineY,lineX,lineY)
	length = coors[length(coors)]-coors[1]
	for (i in GCH_idx)
			segments(-lineX+(2*lineX/length)*(coors[i]-coors[1]),lineY,-lineX+(2*lineX/length)*(coors[i]-coors[1]),lineY-2, col = "darkgreen")
	for (i in HCG_idx)
			segments(-lineX+(2*lineX/length)*(coors[i]-coors[1]),lineY,-lineX+(2*lineX/length)*(coors[i]-coors[1]),lineY-2, col = "black")

	text_cex=0.5
	text(-lineX, lineY+5, paste(chr_name, coors[1], sep = ":"), cex=text_cex)
	text(-lineX+30, lineY-5, coors[round(length(coors)/3)], cex=text_cex)
	text(0, lineY+5, paste("("," bps)", sep = as.character(length+1)), cex=text_cex)
	text(30, lineY+5, coors[round(2*length(coors)/3)], cex=text_cex)
	text(lineX, lineY+5, coors[length(coors)], cex=text_cex)

	text(-lineX+45, lineY-32, "CG", col = "black", cex=0.75, font = 2)		   
	meth <- rep(0, length(coors))
	returned_coors <- draw_NoMe_read(meth, GCH_idx, HCG_idx, -lineX+50, lineY-30, meth, "black", radius = 1)
	text(returned_coors[nrow(returned_coors),1]+10, lineY-32, "GC", col = "darkgreen", cex=0.75, font = 2)		   

	for (i in 1:length(HCG_idx))
		segments(-lineX+(2*lineX/length)*(coors[HCG_idx[i]]-coors[1]),lineY-2, returned_coors[i,1]+1, returned_coors[i,2]+3, col = "black")

	for (i in 1:length(GCH_idx))
			segments(-lineX+(2*lineX/length)*(coors[GCH_idx[i]]-coors[1]),lineY-2, returned_coors[i+length(HCG_idx),1], returned_coors[i+length(HCG_idx),2]+2, col = "darkgreen")
	
	line_space = 4
	radius = 0.4	
	cov_cutoff = 5
	Y = lineY - 50
	X = -lineX - 70
	is_nome = 1
		
	for (i in 1:length(reads)) {
		text(X+5, Y-1, names(reads)[i], font = 2)
		Y = Y-line_space-7
		t2 <- reads[[i]]
		hc <- hclust(as.dist(get_dist(t2)))
		t2 <- t2[hc$order,]
		for (j in 1:nrow(t2)) {
			meths <- rep(-1, time=length(coors))
			for (k in 1:ncol(t2)) {
				if (t2[j,k] %in% c("m", "M")) meths[k] = 1
				if (t2[j,k] %in% c("u", "U")) meths[k] = 0
			}
			draw_NoMe_read(meths, GCH_idx, HCG_idx, X, Y, meths, color, radius = 1, is_nome)
			Y = Y-line_space
		}
			Y = lineY-50
			X = X+CpG_window*2 + 25
	}
						
	dev.off()
}
		

