#!/home/jweile/bin/Rscript

#keeps track of the previously processed jobs
processed <- data.frame()
if (file.exists("www/html/popcodeSuite/processed")) {
	processed <- read.csv("www/html/popcodeSuite/processed",stringsAsFactors=FALSE)
} 

exportProcessedTable <- function() {
	write.table(
		processed,"www/html/popcodeSuite/processed",
		sep=",",row.names=FALSE,quote=FALSE
	)
}

#a function to delete all files belonging to a given id
purge <- function(id) {
	files <- list.files("www/html/popcodeSuite/",pattern=id,full.names=TRUE)
	success <- sapply(files,file.remove)
	if (all(success)) {
		processed <<- processed[-which(processed$id == id),,drop=FALSE]
	}
	exportProcessedTable()
}

#a function to start processing a given id
process <- function(id,path) {
	#check if an options file exists
	opts.file <- paste("www/html/popcodeSuite/",id,"_opts.csv",sep="")
	if (file.exists(opts.file)) {
		#read the options file
		opts <- read.csv(opts.file,stringsAsFactors=FALSE)
		oligo.length <- if ("oligo.length" %in% colnames(opts)) opts$oligo.length else 33
		wiggle <- if ("wiggle" %in% colnames(opts)) opts$wiggle else 5
	}
	#queue up the process in SGE
	system(
		paste(
			"qsub -V -e /dev/null -o /dev/null -cwd -b y ",
			"/home/jweile/bin/Rscript ",
			"projects/popcodeSuite/popcodeSuite.R ",
			"seq=",path," ",
			"outfile=www/html/popcodeSuite/",id," ",
			"length=",oligo.length," ",
			"wiggle=",wiggle," ",
			">www/html/popcodeSuite/",id,".out",
			sep=""
		)
	)
	processed[nrow(processed)+1,"id"] <<- id
	exportProcessedTable()
}

#infinite loop checks for new ids and processes them
while(TRUE) {


	files <- list.files("www/html/popcodeSuite/",pattern="_in.fa")

	invisible(lapply(files, function(f) {

		path <- paste("www/html/popcodeSuite/",f,sep="")
		id <- substr(f,1,47)
		age <- difftime(Sys.time(), file.info(path)[,"ctime"], units = "weeks")

		if (!(id %in% processed$id)) {
			process(id,path)
		}

		if (age >= 1) {
			purge(id)
		}

	}))

	Sys.sleep(2)
}



