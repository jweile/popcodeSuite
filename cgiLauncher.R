#!/home/jweile/bin/Rscript

sink(file("msg.log", open="w"),type="message")

makeUUID <- function() {
	baseuuid <- paste(sample(c(letters[1:6],0:9),30,replace=TRUE),collapse="")
	paste(
	    substr(baseuuid,1,8),"-",
	    substr(baseuuid,9,12),"-","4",
	    substr(baseuuid,13,15),"-",
	    sample(c("8","9","a","b"),1),
	    substr(baseuuid,16,18),"-",
	    substr(baseuuid,19,30),
	    sep="", collapse=""
	)
}

interpolate <- function(filename, values) {
	con <- file(filename,open="r")
	string <- paste(readLines(con),collapse="\n")
	close(con)
	for (name in names(values)) {
		string <- gsub(paste("%",name,sep=""),values[[name]],string)
	}
	string
}

showError <- function(err) {
	cat(interpolate("../../html/popcodeSuite/error.html",c(err=err)))
}

showWait <- function(id) {
	outfile <- paste("../../html/popcodeSuite/",id,".out",sep="")
	if (file.exists(outfile)) {
		f <- file(outfile,open="r")
		output <- readLines(f)
		close(f)
	} else {
		output <- "Waiting for daemon..."
	}
	cat(interpolate("../../html/popcodeSuite/wait.html",
		c(id=id, output=paste(output,collapse="\n"))
	))
}

showResult <- function(id) {
	f <- file(paste("../../html/popcodeSuite/",id,".fa",sep=""),open="r")
	oligos <- readLines(f)
	close(f)

	cat(interpolate("../../html/popcodeSuite/result.html",
		c(oligos=paste(oligos,collapse="\n"), id=id)
	))
}

showInput <- function() {
	cat(interpolate("../../html/popcodeSuite/cgiInput.html",character(0)))
}

#read GET data
getTable <- do.call(rbind,strsplit(strsplit(Sys.getenv("QUERY_STRING"),"&")[[1]],"="))
getData <- lapply(getTable[,2],URLdecode)
names(getData) <- getTable[,1]

#read POST data
f <- file("stdin",open="r")
postRaw <- readLines(f)
close(f)
if (!is.null(postRaw) && length(postRaw) > 0 && nchar(postRaw) > 0) {
  postTable <- do.call(rbind,strsplit(strsplit(postRaw,"&")[[1]],"="))
  postData <- lapply(postTable[,2],URLdecode)
  names(postData) <- postTable[,1]
} else {
  postData <- NULL
}

#has an ID been provided?
if ("id" %in% names(getData)) {
	#if an id was quoted, it's not a new job
	id <- gsub(" ","",getData[["id"]])
	#check if output file exists
	if (file.exists(paste("../../html/popcodeSuite/",id,"_mutcov.png",sep=""))) {
		#then the job is done
		showResult(id)
	} else if (!file.exists(paste("../../html/popcodeSuite/",id,"_in.fa",sep=""))) {
		#invalid id
		showError("This job does not exist (anymore?)")
	} else {
		#the job is not yet done.
		#display wait message and refresh
		showWait(id)
	}
} else {
	#no ID was provided -> it's a new job
	if ("seq" %in% names(postData)) {
		id <- paste(format(Sys.time(),"%Y-%m-%d"),makeUUID(),sep="_")
		#get the input data
		seq.data <- postData[["seq"]]
		if (!is.null(seq.data) && length(seq.data) == 1 && nchar(seq.data) > 0) {
			#check validity
			fasta.match <- gregexpr(
				"(>\\w+\\s*\\n([ACGTNRYSWKMacgtnryswkm]+\\s*\\n?)+)+",
				seq.data,perl=TRUE
			)
			if (length(fasta.match) == 1 && 
				fasta.match[[1]] == 1 && 
				attr(fasta.match[[1]],"match.length")==nchar(seq.data)) {

				#write input file
				f <- file(paste("../../html/popcodeSuite/",id,"_in.fa",sep=""),open="w")
				writeLines(seq.data,f)
				close(f)
				#Daemon will pick up the input file and start processing.

				#show wait message and refresh
				showWait(id)

			} else {
				showError("Input is not a valid FASTA document!")
			}

		} else {
			showError("No input provided!")
		}

	} else {
		showInput()
	}

}









