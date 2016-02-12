library("DNAcopy")

lowess.gc <- function(jtkx, jtky) {
        jtklow <- lowess(jtkx, log(jtky), f=0.05)
        jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        return(exp(log(jtky) - jtkz$y))
}


cbs.segment01 <- function(indir, outdir, varbin.gc, varbin.bad, varbin.data, sample.name, alt.sample.name, alpha, nperm, undo.SD, min.width) {
	gc <- read.table(varbin.gc, header=T)

	chrom.numeric <- substring(gc$bin.chrom, 4)
	chrom.numeric[which(gc$bin.chrom == "chrX")] <- "23"
	chrom.numeric[which(gc$bin.chrom == "chrY")] <- "24"
	chrom.numeric <- as.numeric(chrom.numeric)

	thisRatio <- read.table(paste(indir, varbin.data, sep="/"), header=F) 
	names(thisRatio) <- c("chrom", "chrompos", "abspos", "bincount", "ratio")
	thisRatio$chrom <- chrom.numeric
       	a <- thisRatio$bincount + 1
	#### a <- thisRatio$bincount
	thisRatio$ratio <- a / mean(a[which(chrom.numeric < 23)])
	thisRatio$gc.content <- gc$gc.content
	thisRatio$lowratio <- lowess.gc(thisRatio$gc.content, thisRatio$ratio)
	#thisRatio$log2lowratio <- log(thisRatio$lowratio, base=2) 

	set.seed(25) 
	CNA.object <- CNA(log(thisRatio$lowratio, base=2), gc$bin.chrom, thisRatio$chrompos, data.type="logratio", sampleid=sample.name) 
	smoothed.CNA.object <- smooth.CNA(CNA.object) 
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=2) 
	thisShort <- segment.smoothed.CNA.object[[2]]

	sortcol <- thisShort$chrom
	sortcol <- gsub("chr", "", sortcol)
	sortcol <- gsub("p", "", sortcol)
	sortcol <- gsub("q", "", sortcol)
	thisShort <- thisShort[order(as.numeric(sortcol)), ]

	m <- matrix(data=0, nrow=nrow(thisRatio), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	thisRatio$seg.mean.LOWESS <- m[, 1]

	thisGrid <- seq(1.5, 5.5, by=0.05)
	thisOuter <- thisRatio$seg.mean.LOWESS %o% thisGrid
	thisOuterRound <- round(thisOuter)
	thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
	thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
	thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
	thisError <- min(thisOuterColsums)
	thisShredded <- length(which(thisRatio$seg.mean.LOWESS[which(chrom.numeric < 23)] < 0.1)) / length(which(chrom.numeric < 23))

	thisRatio$ratio.quantal <- thisRatio$lowratio * thisMultiplier
	thisRatio$seg.quantal <- thisRatio$seg.mean.LOWESS * thisMultiplier
	
	thisQuantalStats <- data.frame("ploidy"=thisMultiplier, "error"=thisError, "shredded"=thisShredded)
	
	chr <- thisRatio$chrom
	chr.shift <- c(chr[-1], chr[length(chr)])
	vlines <- c(1, thisRatio$abspos[which(chr != chr.shift) + 1], thisRatio$abspos[nrow(thisRatio)])
	hlines <- c(0.5, 1.0, 1.5, 2.0)
	chr.text <- c(1:22, "X", "Y")
	vlines.shift <- c(vlines[-1], 4*10^9)
	chr.at <- vlines + (vlines.shift - vlines) / 2
	x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
	x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")

	png(paste(outdir, "/", sample.name, ".png", sep=""), height=800, width=1200)
	plot(x=thisRatio$abspos, y=thisRatio$lowratio, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
	axis(1, at=x.at, labels=x.labels)
	lines(x=thisRatio$abspos, y=thisRatio$lowratio, col="#CCCCCC")
	points(x=thisRatio$abspos, y=thisRatio$seg.mean.LOWESS, col="#0000AA")
	lines(x=thisRatio$abspos, y=thisRatio$seg.mean.LOWESS, col="#0000AA")
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()
	
	hlines <- c(1, 2, 3, 4, 5, 6)

	png(paste(outdir, "/", sample.name, ".quantal.png", sep=""), height=800, width=1200)
	plot(x=thisRatio$abspos, y=thisRatio$ratio.quantal, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
	axis(1, at=x.at, labels=x.labels)
	lines(x=thisRatio$abspos, y=thisRatio$ratio.quantal, col="#CCCCCC")
	points(x=thisRatio$abspos, y=thisRatio$seg.quantal, col="#0000AA")
	lines(x=thisRatio$abspos, y=thisRatio$seg.quantal, col="#0000AA")
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()

	write.table(thisQuantalStats, sep="\t", file=paste(outdir, "/", sample.name, ".varbin.quantal.stats.txt", sep=""), quote=F, row.names=F) 
	write.table(thisRatio, sep="\t", file=paste(outdir, "/", sample.name, ".varbin.data.txt", sep=""), quote=F, row.names=F) 
	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".varbin.short.txt", sep=""), quote=F, row.names=F) 

        # NOBAD SECTION
	bad <- read.table(varbin.bad, header=F, as.is=T, stringsAsFactors=F)

	thisRatioNobig <- thisRatio[-bad[, 1], ]

	set.seed(25) 
	CNA.object <- CNA(log(thisRatioNobig$lowratio, base=2), gc$bin.chrom[-bad[, 1]], thisRatioNobig$chrompos, data.type="logratio", sampleid=sample.name) 
	smoothed.CNA.object <- smooth.CNA(CNA.object) 
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=2) 
	thisShort <- segment.smoothed.CNA.object[[2]]

	sortcol <- thisShort$chrom
	sortcol <- gsub("chr", "", sortcol)
	sortcol <- gsub("p", "", sortcol)
	sortcol <- gsub("q", "", sortcol)
	thisShort <- thisShort[order(as.numeric(sortcol)), ]

	m <- matrix(data=0, nrow=nrow(thisRatioNobig), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	thisRatioNobig$seg.mean.LOWESS <- m[, 1]

	thisGrid <- seq(1.5, 5.5, by=0.05)
	thisOuter <- thisRatioNobig$seg.mean.LOWESS %o% thisGrid
	thisOuterRound <- round(thisOuter)
	thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
	thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
	thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
	thisError <- min(thisOuterColsums)
	thisShredded <- length(which(thisRatioNobig$seg.mean.LOWESS[which(chrom.numeric < 23)] < 0.1)) / length(which(chrom.numeric < 23))

	thisRatioNobig$ratio.quantal <- thisRatioNobig$lowratio * thisMultiplier
	thisRatioNobig$seg.quantal <- thisRatioNobig$seg.mean.LOWESS * thisMultiplier

	
	thisQuantalStatsNobad <- data.frame("ploidy"=thisMultiplier, "error"=thisError, "shredded"=thisShredded)
	
	chr <- thisRatioNobig$chrom
	chr.shift <- c(chr[-1], chr[length(chr)])
	# Use the vlines abspos positions from above.  Because these start at
	# the third bin on the acrocentric chromosomes the vlines end up to
	# the right of the centromere rather than the left which is wrong.
	#vlines <- c(1, thisRatioNobig$abspos[which(chr != chr.shift) + 1], thisRatioNobig$abspos[nrow(thisRatioNobig)])
	hlines <- c(0.5, 1.0, 1.5, 2.0)
	chr.text <- c(1:22, "X", "Y")
	vlines.shift <- c(vlines[-1], 4*10^9)
	chr.at <- vlines + (vlines.shift - vlines) / 2
	x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
	x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")

	png(paste(outdir, "/", sample.name, ".nobad.png", sep=""), height=800, width=1200)
	plot(x=thisRatioNobig$abspos, y=thisRatioNobig$lowratio, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
	axis(1, at=x.at, labels=x.labels)
	lines(x=thisRatioNobig$abspos, y=thisRatioNobig$lowratio, col="#CCCCCC")
	points(x=thisRatioNobig$abspos, y=thisRatioNobig$seg.mean.LOWESS, col="#0000AA")
	lines(x=thisRatioNobig$abspos, y=thisRatioNobig$seg.mean.LOWESS, col="#0000AA")
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()
	
	hlines <- c(1, 2, 3, 4, 5, 6)

	png(paste(outdir, "/", sample.name, ".nobad.quantal.png", sep=""), height=800, width=1200)
	plot(x=thisRatioNobig$abspos, y=thisRatioNobig$ratio.quantal, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
	axis(1, at=x.at, labels=x.labels)
	lines(x=thisRatioNobig$abspos, y=thisRatioNobig$ratio.quantal, col="#CCCCCC")
	points(x=thisRatioNobig$abspos, y=thisRatioNobig$seg.quantal, col="#0000AA")
	lines(x=thisRatioNobig$abspos, y=thisRatioNobig$seg.quantal, col="#0000AA")
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()

	write.table(thisQuantalStatsNobad, sep="\t", file=paste(outdir, "/", sample.name, ".nobad.varbin.quantal.stats.txt", sep=""), quote=F, row.names=F) 
	write.table(thisRatioNobig, sep="\t", file=paste(outdir, "/", sample.name, ".nobad.varbin.data.txt", sep=""), quote=F, row.names=F) 
	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".nobad.varbin.short.txt", sep=""), quote=F, row.names=F) 



}

# Parse command line - no spaces in arguments
args=(commandArgs(TRUE))
for (i in 1:length(args)) {
    eval(parse(text=args[i]))
}


if (!exists("sample")) {
   stop("missing command line argument: sample")
}
if (!exists("gc")) {
   stop("missing command line argument: gc")
}
if (!exists("bad")) {
   stop("missing command line argument: bad")
}

# alpha 0.05 undo.sd 1.0 good for normal profile
# alpha 0.02 undo.sd 0.5 good for cancer
cbs.segment01(indir="./", outdir="./", varbin.gc=gc, varbin.bad=bad, varbin.data="varbin.txt", sample.name=sample, alt.sample.name="", alpha=0.02, nperm=1000, undo.SD=1.0, min.width=3)





