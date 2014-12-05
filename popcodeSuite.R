#!/home/jweile/bin/Rscript

if (!("Biostrings" %in% installed.packages()[,"Package"])) {
	source("http://bioconductor.org/biocLite.R")
	biocLite("Biostrings")
}
if (!("HELP" %in% installed.packages()[,"Package"])) {
	source("http://bioconductor.org/biocLite.R")
	biocLite("HELP")
}

library("Biostrings")
library("HELP")

###
# get the user-supplied argument with the given name
# if not defined return default value
getArg <- function(name, default=NULL, required=FALSE) {

	if (length(commandArgs(TRUE)) == 0) {
		if (required) {
			stop("Required argument:",name)
		} else {
			return(default)
		}
	}

	#tabulate arguments by name
	argTable <- do.call(rbind,strsplit(commandArgs(TRUE),"="))
	#get index of argument with given name
	i <- which(argTable[,1] == name)


	if (length(i) == 0) {
		#return default value if no matching arguments are found.
		if (required) {
			stop("Required argument:",name)
		} else {
			return(default)
		}
	} else if (length(i) > 1) {
		#if multiple matches are found, throw error message.
		stop("Multiple values for", name, "found!")
	} else {
		#if everything checks out, return the argument value
		return(argTable[i,2])
	}
}

###
# get sequence with given name from sequence list 
getSeqByName <- function(seqs, name) {
	seqs[[which(names(seqs)==name)]]
}

#GET INPUT ARGUMENTS

oligo.length <- as.numeric(getArg("length",33))
wiggle <- as.numeric(getArg("wiggle",5))
input.seqs <- readDNAStringSet(getArg("seq",required=TRUE))
outfile <- getArg("outfile","output/tmplDirectedNNK")

p.intercept <- as.numeric(getArg("p.intercept",0.0693407))
p.coefficient <- as.numeric(getArg("p.coefficient",-0.0008216))

max.clones <- as.numeric(getArg("maxClones",10000))

#build sequence construct
orf <- as.character(getSeqByName(input.seqs,"ORF"))
prefix <- as.character(getSeqByName(input.seqs,"prefix"))
suffix <- as.character(getSeqByName(input.seqs,"suffix"))
construct <- paste(prefix,orf,suffix,sep="")

cat("Construct:",construct,"\n\n")

cat("Exploring possible oligos...\n")
codon.starts <- nchar(prefix) + seq(1+3,nchar(orf),3)
oligo.choices <- lapply(codon.starts, function(codon.start) {
	lapply(-wiggle:wiggle, function(offset) {
		start <- codon.start + 1 - floor(oligo.length/2) + offset
		end <- codon.start + 1 + floor(oligo.length/2) + offset
		sequence <- substr(construct,start,end)
		list(
			codon.start.in.construct=codon.start,
			codon.start.in.oligo=codon.start-start+1,
			oligo.start=start,
			oligo.end=end,
			sequence=sequence,
			tm=calcTm(sequence)
		)
	})
})
cat("Optimizing melting temperatures...\n")
median.tm <- median(unlist(lapply(oligo.choices, lapply, function(item)item$tm)))
best.oligos <- do.call(rbind,lapply(oligo.choices, function(options) {
	tms <- sapply(options, function(item) item$tm)
	min.idx <- which.min(abs(tms - median.tm))
	options[[min.idx]]
}))

cat("Plotting offset distribution...\n")
pdf(paste(outfile,".pdf",sep=""),width=6,height=3)
op <- par(mfrow=c(1,2),cex=.9)

barplot(
	table(unlist(best.oligos[,"codon.start.in.oligo"])-floor(oligo.length/2)),
	xlab="Offset",
	ylab="Frequency",
	col="gold2",
	border="gold4"
)
grid(nx=NA,ny=NULL)

hist(
	unlist(best.oligos[,"tm"]),
	breaks=20,
	col="steelblue3",
	border="steelblue4",
	xlab="Melting temperature (C)",
	main=""
)
grid(nx=NA,ny=NULL)
tmm <- median(unlist(best.oligos[,"tm"]))
abline(v=tmm,lty="dashed",col="darkgray")
text(tmm, (par("usr")[4]-par("usr")[3])/2, paste(format(tmm,digits=4),"C"))

par(op)
invisible(dev.off())

cat("inserting NNKs...\n")

oligos <- do.call(c, lapply(1:nrow(best.oligos), function(i) {
	with(best.oligos[i,], {
		codon <- substr(sequence, codon.start.in.oligo, codon.start.in.oligo+2)
		aa <- as.character(translate(DNAString(codon)))
		outseq <- sequence
		substr(outseq, codon.start.in.oligo, codon.start.in.oligo+2) <- "NNK"
		outseqs <- c(sequence,outseq)
		names(outseqs) <- paste(aa,i,c("wt","X"),sep="")
		cat(i,"\t",substr(sequence,1,1),"\t",substr(sequence,nchar(sequence),nchar(sequence)),"\n")
		outseqs
	})
}))


cat("Writing output...\n\n")

#write results to file
writeXStringSet(DNAStringSet(oligos), paste(outfile,".fa",sep=""), format="fasta")

mut.prob <- p.intercept + p.coefficient * do.call(c,best.oligos[,"tm"])
mut.prob <- mut.prob / sum(mut.prob)


cat("Plotting predicted mutation rates...\n")
pdf(paste(outfile,"_mutrate.pdf",sep=""),width=12,height=3)
xs <- barplot(
	mut.prob,
	xlab="Position",ylab="P(mut)",
	col="steelblue2",border="gray",
	main="Predicted mutation probabability"
)
labels <- c(1,seq(from=20,to=length(xs),by=20))
axis(1,at=xs[labels],labels=labels)
invisible(dev.off())

aa.seq <- translate(DNAString(orf))
names(mut.prob) <- 2:(length(mut.prob)+1)

###
# Implementation of the magic table algorithm to sample from custom probability functions.
# n = number of samples to draw
# tab = named numberic vector representing the probability of each possible value.
magic.table <- function(n, tab) {
	if (any(tab < 0)) stop("No negative values allowed!")
	cumul <- sapply(1:length(tab),function(i) sum(tab[1:i]))
	is <- sapply(runif(n,0,sum(tab)),function(i) min(which(cumul >= i)))
	if (is.null(names(tab))) is else names(tab)[is]
}

# #samples the number of mutations as observed in PopCode data
# rnmuts <- function(n) {
# 	sapply(
# 		rpois(n,lambda=as.numeric(magic.table(n,c(`1.1`=17/24,`7`=2/24,`NA`=5/25)))),
# 		function(x) if(is.na(x)) "ns" else x
# 	)
# }

#samples the number of mutations as observed in PopCode data
rnmuts <- function(n, orf.length) {
	rate <- orf.length * 0.004240506
	tab <- c(.52,.48)
	names(tab) <- c(rate,"NA")
	sapply(
		rpois(n,lambda=as.numeric(magic.table(n,tab))),
		function(x) if(is.na(x)) "ns" else x
	)
}


#simulate mutagenesis. Returns a list of variants. The variants being vectors of mutations.
simulate.mutagenesis <- function(aa.seq, n, pos.p) {
	aa.seq <- as.character(aa.seq)
	aas <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
	lapply(rnmuts(n,nchar(aa.seq)*3), function(nmuts) {
		if (nmuts != "ns") {
			nmuts <- as.numeric(nmuts)
			# pos <- sort(sample(1:nchar(aa.seq),nmuts))
			pos <- if (nmuts > 0) as.numeric(magic.table(nmuts,pos.p)) else numeric(0)
			sapply(pos, function(p) {
				from <- substr(aa.seq,p,p)
				to <- sample(setdiff(aas,from),1)
				paste(from,p,to,sep="")
			})
		} else {
			"frameshift"
		}
	})
}

calc.cov <- function(clones, num.aa) {
	change.matrix <- matrix(0,nrow=21,ncol=num.aa,
		dimnames=list(
			c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'),
			1:num.aa
		)
	)
	for (clone in clones) {
		if (length(clone) > 1 && clone != "frameshift") {
			for (mut in clone) {
				pos <- as.numeric(substr(mut,2,nchar(mut)-1))
				to <- substr(mut,nchar(mut),nchar(mut))
				change.matrix[to,pos] <- change.matrix[to,pos]+1
			}
		}
	}
	mean(apply(change.matrix,2,function(col) sum(col > 0))/19)
}

cat("Simulating mutagenesis...\n")

ns <- seq(1000,max.clones,1000)
clones <- NULL
coverages <- sapply(ns,function(n) {
	cat("Simulating",n,"clones\n")
	clones <<- simulate.mutagenesis(aa.seq, n, mut.prob)
	calc.cov(clones,nchar(aa.seq))
})

cat("Plotting predicted coverage...\n")
pdf(paste(outfile,"_clones.pdf",sep=""),width=4,height=4)
plot(
	ns,coverages,type="l",
	ylim=c(0,1),
	xlab="Picked colonies",ylab="Mutant coverage"
)
grid()
dev.off()

data.frame(n=ns,coverage=coverages)


# Plots the mutation coverage for a given change matrix
plotMutCoverage <- function(mutations, sequence, all=FALSE, main="") {

	num.aa <- nchar(sequence)/3

	nucls <- lapply(c("A","C","G","T"),DNAString)
	levenstein1 <- lapply(1:num.aa, function(pos) {
		codon.start <- 1+(pos-1)*3
		codon <- subseq(sequence,codon.start,codon.start+2)
		aa <- as.character(translate(codon))
		aas <- do.call(c,lapply(1:3,function(i){
			sapply(nucls,function(nucl) {
				.codon <- codon
				subseq(.codon,i,i) <- nucl
				as.character(translate(.codon))
			})
		}))
		setdiff(unique(aas),c(aa,"*"))
	})

	#initialize the change matrix
	change.matrix <- matrix(0,nrow=21,ncol=num.aa,
		dimnames=list(
			c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'),
			1:num.aa
		)
	)
	
	#remove any mutations that come with trucations and nonsense
	muts <- unlist(mutations[sapply(mutations, {
		function(x) length(x) > 0 && !any(x == "frameshift")
	})])

	# Mark mutations in the matrix
	for (sac in muts) {

		if (sac == "silent") next

		pos <- as.numeric(substr(sac,2,nchar(sac)-1))
		aa <- substr(sac,nchar(sac),nchar(sac))
		from <- substr(sac,1,1)

		change.matrix[aa,pos] <- change.matrix[aa,pos] + 1
		change.matrix[from,pos] <- -1
	}

	coverage.all <- apply(change.matrix,2,function(x) sum(na.omit(x) > 0)) / 19
	coverage.lev <- sapply(1:num.aa, function(i) {
		sum(change.matrix[levenstein1[[i]],i] > 0)/length(levenstein1[[i]])
	})

	#define drawing layot, set drawing color to gray, adjust margins
	layout(cbind(1:3,c(5,5,4)), heights=c(1,.6,3),widths=c(9,.5))
	op <- par(fg="gray",mar=c(0,4.1,4.1,2.1),xaxs="i")	

	# draw a bar plot for coverage.all
	barplot(coverage.all,
		main=main, 
		xlab="Position",
		ylab="Cover. All", 
		ylim=c(0,1),
		border=NA,
		names.arg=NA,
		col="darkolivegreen3",
		axes=FALSE
	)
	axis(2,at=c(0,.5,1),labels=c(0,.5,1))
	grid(NA,4)

	op <- par(mar=c(0,4.1,1,2.1))	
	barplot(coverage.lev,
		xlab="Position",
		ylab="Cover. Lev1", 
		ylim=c(0,1),
		border=NA,
		names.arg=NA,
		col="darkolivegreen4",
		axes=FALSE
	)
	axis(2,at=c(0,.5,1),labels=c(0,.5,1))
	grid(NA,4)


	# Compute a color gradient to represent the mutation counts
	maxVal <- max(apply(change.matrix,1,function(x) max(na.omit(x))))
	colors <- colorRampPalette(c("white", "orange"))(5)

	### Draw the diagram
	# use horizontal axis labels
	op <- c(op,par(las=1))
	par(mar=c(5.1,4.1,0,2.1),xaxs="i")
	# create an empty plot
	plot(0,
		type='n',
		axes=FALSE,
		xlim=c(0,ncol(change.matrix)), 
		ylim=c(0,21),
		xlab="Position",
		ylab="Amino acid"
	)
	# iterate over each matrix entry and draw the contents on the plot
	for (x in 1:ncol(change.matrix)) {
		for (y in 1:21) {
			if (change.matrix[y,x] > 0) {
				#observed mutations are drawn in a color shade corresponding to their count
				col <- colors[ceiling(4*change.matrix[y,x]/maxVal)+1]
				rect(x-1,22-y,x,21-y,col=col, lty="blank")
			} else if (change.matrix[y,x] < 0) {
				rect(x-1,22-y,x,21-y,col="gray", lty="blank")
			}
		}
	}
	# draw axes
	axis(1, at=c(1,seq(5,ncol(change.matrix),5))-.5, labels=c(1,seq(5,ncol(change.matrix),5)))
	axis(2, at=(1:21)-.5, labels=rev(rownames(change.matrix)))

	par(op)

	op <- par(mar=c(5.1,0,0,4.1),las=1)
	plot(0,type="n",ylim=c(0,6),xlim=c(0,1),axes=FALSE,ylab="",xlab="")
	tops <- 1:4 * maxVal/4
	bottoms <- round(tops - maxVal/4 + 1)
	tops <- round(tops)
	rect(0,1:4,1,2:5,col=colors[2:5],lty="blank")
	rect(0,5,1,6,col="gray",lty="blank")
	axis(4,at=0:5+.5,labels=c("0",paste(bottoms,tops,sep="-"),"wt"),tick=FALSE)
}


cat("Plotting detailed coverage...\n")
pdf(paste(outfile,"_mutcov.pdf",sep=""),width=14,height=4)
plotMutCoverage(clones,DNAString(orf))
cat("Done!\n")
