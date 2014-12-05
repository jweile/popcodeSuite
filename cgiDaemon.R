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
		processed[which(processed$id == id),"id"] <<- NULL
	}
	exportProcessedTable()
}

#a function to start processing a given id
process <- function(id,path) {
	system(
		paste(
			"projects/popcodeSuite/popcodeSuite.R ",
			"seq=",path," ",
			"outfile=www/html/popcodeSuite/",id," ",
			"&>www/html/popcodeSuite/",id,".out",
			sep=""
		),
		wait=FALSE
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



