#!/usr/bin/env Rscript

#keeps track of the previously processed jobs
status <- data.frame(id=character(0),status=character(0),stringsAsFactors=FALSE)
if (file.exists("status")) {
	status <- read.csv("status",stringsAsFactors=FALSE)
} 

user <- "jweile"
host <- "guru-int.mshri.on.ca"
workspace <- "/home/rothlab/jweile/pcs_workspace/"
script <- "/home/rothlab/jweile/projects/popcodeSuite/popcodeSuite.R"

#save the status table to disk
save.status <- function() {
	write.table(
		status,"status",
		sep=",",row.names=FALSE,quote=FALSE
	)
}

#set a status for the given id
set.status <- function(id,st) {
	if (id %in% status$id) {
		status[which(status$id==id),"status"] <<- st
	} else {
		status[nrow(status)+1,] <<- list(id,st)
	}
}

#retrieve the status for the given id
get.status <- function(id) {
	if (id %in% status$id) {
		status[which(status$id==id),"status"]
	} else {
		"new"
	}
}


#a function to delete all files belonging to a given id
purge <- function(id) {
	files <- list.files(".",pattern=id,full.names=TRUE)
	success <- sapply(files,file.remove)
	if (all(success)) {
		status <<- status[-which(status$id == id),,drop=FALSE]
	}
	exportProcessedTable()
}

#start a new job for a given id and input file
process.job <- function(id,f) {

	#read the options file
	opts.file <- paste0(id,"_opts.csv")
	if (file.exists(opts.file)) {
		opts <- read.csv(opts.file,stringsAsFactors=FALSE)
		oligo.length <- if ("oligo.length" %in% colnames(opts)) opts$oligo.length else 33
		wiggle <- if ("wiggle" %in% colnames(opts)) opts$wiggle else 5
		version <- if ("version" %in% colnames(opts)) opts$version else 1
	} else {
		oligo.length <- 33
		wiggle <- 5
		version <- 1
	}

	#SCP input file to cluster host
	errCode <- system(paste0(
		"rsync -az ",f," ",user,"@",host,":",workspace
	),intern=FALSE)
	if (errCode) {
		warning("SCP failed!")
		set.status(id,"error")
		return()
	}
	#Use SSH to schedule an SGE job on the data
	errCode <- system(paste0(
		"ssh ",user,"@",host,
		" \"source /etc/profile; qsub -V ",
		"-e ",workspace,id,".err ",
		"-o ",workspace,id,".out ",
		"-wd ",workspace," -b y ",
		"/software/R/bin/Rscript ",
		script," ",
		"seq=",f," ",
		"outfile=",workspace,id," ",
		"length=",oligo.length," ",
		"wiggle=",wiggle," ",
		"version=",version," ",
		"\""
	),intern=FALSE)
	if (errCode) {
		warning("SSH failed!")
		set.status(id,"error")
		return()
	}

	set.status(id,"live")

}

update.job <- function(id) {

	#check if log file exists
	log.exists <- as.logical(system(paste0(
		"ssh ",user,"@",host,
		" \"[[ -f ",workspace,id,".out ]] && echo TRUE || echo FALSE\""
	),intern=TRUE))
	if (log.exists) {
		#SCP the log file over
		errCode <- system(paste0(
			"rsync -az ",user,"@",host,":",workspace,id,".out ."
		),intern=FALSE)
		if (errCode) {
			warning("Unable to transfer remote logfile!")
		} else {
			#check if it's done
			if (system(paste0("tail -1 ",id,".out"),intern=TRUE) == "Done!") {

				#copy over the results
				errCode <- system(paste0(
					"rsync -az --exclude '*_in.fa' ",user,"@",host,":",workspace,id,"* ."
				),intern=FALSE)
				if (errCode) {
					warning("SCP failed!")
					set.status(id,"error")
					return()
				}

				#delete remote copy of results
				errCode <- system(paste0(
					"ssh ",user,"@",host," rm ",workspace,id,"*"
				),intern=FALSE)
				if (errCode) {
					warning("Unable to delete remote copy!")
				}

				#mark as done
				set.status(id,"done")
			}
		}
	}

}




#infinite loop checks for new ids and processes them
while(TRUE) {

	invisible(lapply(list.files(".",pattern="_in.fa"), function(f) {

		id <- substr(f,1,47)
		age <- difftime(Sys.time(), file.info(f)[,"ctime"], units = "weeks")
		st <- get.status(id)

		if (st == "new") {
			process.job(id,f)
		} else if (st == "live") {
			update.job(id)
		} 

		if (age >= 1) {
			purge(id)
		}

	}))

	Sys.sleep(2)
}










