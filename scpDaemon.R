#!/home/jweile/bin/Rscript

#keeps track of the previously processed jobs
status <- data.frame(id=character(0),status=character(0),stringsAsFactors=FALSE)
if (file.exists("status")) {
	status <- read.csv("status",stringsAsFactors=FALSE)
} 

save.status <- function() {
	write.table(
		status,"status",
		sep=",",row.names=FALSE,quote=FALSE
	)
}

set.status <- function(id,st) {
	if (id %in% status$id) {
		status[which(status$id==id),"status"] <<- st
	} else {
		status[nrow(status)+1,] <<- list(id,st)
	}
}

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

process.job <- function(id,f) {

	#read the options file
	opts.file <- paste0("www/html/popcodeSuite/",id,"_opts.csv")
	if (file.exists(opts.file)) {
		opts <- read.csv(opts.file,stringsAsFactors=FALSE)
		oligo.length <- if ("oligo.length" %in% colnames(opts)) opts$oligo.length else 33
		wiggle <- if ("wiggle" %in% colnames(opts)) opts$wiggle else 5
		version <- if ("version" %in% colnames(opts)) opts$version else 1
	}

	#SCP input file to guru
	errCode <- system(paste0("scp ",f," guru:/home/jweile/pcs_workspace/"),intern=FALSE)
	if (errCode) {
		warning("SCP failed!")
		set.status(id,"error")
		return()
	}
	#Use SSH to schedule an SGE job on the data
	errCode <- system(paste0(
		"ssh \"source /etc/profile; qsub -V ",
		"-e /home/jweile/pcs_workspace/",id,".err ",
		"-o /home/jweile/pcs_workspace/",id,".out ",
		"-wd /home/jweile/pcs_workspace/ -b y ",
		"/software/R/bin/Rscript ",
		"/home/jweile/projects/popcodeSuite/popcodeSuite.R ",
		"seq=",f," ",
		"outfile=/home/jweile/pcs_workspace/",id," ",
		"length=",oligo.length," ",
		"wiggle=",wiggle," ",
		"version=",version," ",
		"\"",
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
	log.exists <- as.logical(system(paste0("ssh guru \"[[ -f /home/jweile/pcs_workspace/",id,".out ]] && echo TRUE || echo FALSE\""),intern=TRUE))
	if (log.exists) {
		#SCP the log file over
		system(paste0("scp guru:/home/jweile/pcs_workspace/",id,".out ."))
		#check if it's done
		if (system(paste0("tail -1 ",id,".out"),intern=TRUE) == "Done!") {

			#copy over the results
			system("scp guru:/home/jweile/pcs_workspace/",id,"* .")
			#mark as done
			set.status(id,"done")
		}
	}

}




#infinite loop checks for new ids and processes them
while(TRUE) {

	invisible(lapply(list.files(".",pattern="_in.fa"), function(f) {

		id <- substr(f,1,47)
		age <- difftime(Sys.time(), file.info(f)[,"ctime"], units = "weeks")
		status <- get.status(id)

		if (status == "new") {
			process.job(id,f)
		} else if (status == "live") {
			update.job(id)
		} 

		if (age >= 1) {
			purge(id)
		}

	}))

	Sys.sleep(2)
}










