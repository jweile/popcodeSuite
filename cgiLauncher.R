#!/usr/bin/env Rscript

sink(file("msg.log", open="w"),type="message")

docBase <- "/srv/bundler/Staging/popcodeSuite/html/"

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
	cat(interpolate(paste0(docBase,"error.html"),c(err=err)))
}

showWait <- function(id) {
	outfile <- paste(docBase,id,".out",sep="")
	if (file.exists(outfile)) {
		f <- file(outfile,open="r")
		output <- readLines(f)
		close(f)
		if (is.null(output) || length(output) == 0 || nchar(output) == 0) {
			output <- "Waiting for Sun Grid Engine..."
		}
	} else {
		output <- "Scheduling job..."
	}
	cat(interpolate(paste0(docBase,"wait.html"),
		c(id=id, output=paste(output,collapse="\n"))
	))
}

showResult <- function(id) {
	f <- file(paste(docBase,id,".tsv",sep=""),open="r")
	otable <- readLines(f)
	close(f)
	otable <- paste(
		"<table class=\"oligotable\">",
		paste(paste(
			"<tr><td>",
			gsub("\t","</td><td>",otable),
			"</td></tr>",
			sep=""
		),collapse="\n"),
		"</table>",sep="\n"
	)

	if (file.exists(paste(docBase,id,"_mutcov.png",sep=""))) {
		v1mode <- ""
	} else {
		v1mode <- "disabled"
	}

	cat(interpolate(paste0(docBase,"result.html"),
		c(otable=paste(otable,collapse="\n"), id=id, v1mode=v1mode)
	))
}

showInput <- function() {
	cat(interpolate(paste0(docBase,"cgiInput.html"),character(0)))
}


seq.rx <- "([ACGTNRYSWKMacgtnryswkm]+\\n?)+"
valid <- function(x) {
	if (!is.null(x) && length(x) == 1 && nchar(x) > 0) {
		seq.match <- gregexpr(seq.rx, x, perl=TRUE)
		if (length(seq.match) == 1 && 
				seq.match[[1]] == 1 && 
				attr(seq.match[[1]],"match.length")==nchar(x)) {
			return(TRUE)
		}
	}
	return(FALSE)
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
	if (file.exists(paste(docBase,id,".out",sep=""))) {
		out.tail <- system(paste("tail -1 ",docBase,id,".out",sep=""),intern=TRUE)
		if (!is.null(out.tail) && length(out.tail) > 0 && out.tail == "Done!") {
			#then the job is done
			showResult(id)
		} else {
			#still running
			showWait(id)
		}
	} else if (!file.exists(paste(docBase,id,"_in.fa",sep=""))) {
		#invalid id
		showError("This job does not exist (anymore?)")
	} else {
		#the job is not yet done.
		#display wait message and refresh
		showWait(id)
	}
} else {
	#no ID was provided -> it's a new job
	if (any(c("orf","prefix","suffix") %in% names(postData))) {
		#get the input data
		orf <- postData[["orf"]]
		prefix <- postData[["prefix"]]
		suffix <- postData[["suffix"]]
		oligo.length <- postData[["oligo.length"]]
		wiggle <- postData[["wiggle"]]
		ver <- postData[["version"]]

		#check validity	
		if (valid(orf) && valid(prefix) && valid(suffix)) {

			#create ID for this job
			id <- paste(format(Sys.time(),"%Y-%m-%d"),makeUUID(),sep="_")

			#process options and store
			opts <- data.frame()
			if (!is.null(oligo.length) && length(oligo.length)==1 && 
				nchar(oligo.length)>0 && !is.na(as.numeric(oligo.length))) {
				oligo.length <- as.numeric(oligo.length)
				opts[1,"oligo.length"] <- oligo.length
			}
			if (!is.null(wiggle) && length(wiggle)==1 && 
				nchar(wiggle)>0 && !is.na(as.numeric(wiggle))) {
				wiggle <- as.numeric(wiggle)
				opts[1,"wiggle"] <- wiggle
			}
			if (!is.null(ver) && length(ver)==1 && 
				nchar(ver)>0 && !is.na(as.numeric(ver))) {
				ver <- as.numeric(ver)
				opts[1,"version"] <- ver
			}
			
			write.table(
				opts,
				paste(docBase,id,"_opts.csv",sep=""),
				sep=",",quote=FALSE,row.names=FALSE
			)

			#write input file
			f <- file(paste(docBase,id,"_in.fa",sep=""),open="w")
			writeLines(
				c(
					">prefix",
					prefix,
					">ORF",
					orf,
					">suffix",
					suffix
				),
				f
			)
			close(f)
			#Daemon will pick up the input file and start processing.

			#show wait message and refresh
			showWait(id)


		} else {
			showError("The input is not a valid DNA sequence!")
		}

	} else {
		showInput()
	}

}









