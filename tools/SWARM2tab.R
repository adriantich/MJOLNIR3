
if (!require("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
library(optparse)

opts <- list(
  make_option(c("-s", "--swarm"), type = "character", default = NULL,
              help = "Path to the SWARM output file (clusters)."),
  make_option(c("-t", "--tab"), type = "character", default = NULL,
              help = "Path to the abundance table file."),
  make_option(c("-o", "--output_ESV"), type = "character", default = NULL,
              help = "Path to the output file."),
  make_option(c("-O", "--output_MOTU"), type = "character", default = NULL,
              help = "Path to the output file."),
  make_option(c("-c", "--cores"), type = "integer", default = 1,
              help = "Number of cores to use for parallel processing.")
)
opt <- parse_args(OptionParser(option_list = opts))


SWARM2tab <- function(fileswarm, filetab, cores = 1) {
    # Check if the input files exist
    if (!file.exists(fileswarm)) {
        stop("The SWARM output file does not exist: ", fileswarm)
    }
    if (!file.exists(filetab)) {
        stop("The abundance table file does not exist: ", filetab)
    }
    library(parallel)
    library(dplyr)
    # Read cluster list database   
    message("1/7. Reading SWARM results...")
    swarm_db <- readLines(fileswarm)
    total_swarms <- length(swarm_db)
    message("2/7. Read ", total_swarms," total MOTUs.")
    
    
    clusters <- strsplit(swarm_db, "; ")
    message("4/7. Keeping only information of the sequences ",
            "that form each cluster.")
    clusters <- mclapply(clusters,function(x) {sub(";.*", "", x)},
                         mc.cores = cores)
    names(clusters) <- mclapply(clusters, function(x) x[[1]], mc.cores = cores)
    
    # Read counts database and keep only the needed clusters
    message("6/7. Reading the abundance database. ",
            "This could take a while.")
    motu_seqs_names <- stack(clusters) %>% rename(ID = values, MOTU = ind)
    
    # db <- read.table(filetab, sep = "\t", head = TRUE)
    db <- read.delim(filetab, header=TRUE, sep="\t",
                    comment.char="", check.names=FALSE,
                    stringsAsFactors=FALSE)
    names(db) <- sub("^#OTU ","", names(db))     # remove leading '#'
    numseqs <- nrow(db)
    samples <- names(db)[!names(db)%in%c("ID","sequence")]
    db <- merge(motu_seqs_names, db, by = "ID")
    numseqs_reduced <- nrow(db)
    samples <- length(samples)
    message("Finished reading the Database, which includes ", 
            numseqs, " total unique sequences and ", samples, " samples.\n",
            "Kept only ", numseqs_reduced, " sequences for calculations.")
    
    message("7/7. Calculating the number of reads in every sample ",
            "for each MOTU.")
    db_total <- split(db[, names(db) %in% samples], db$MOTU)
    db_total <- mclapply(
        db_total, 
        function(x) {
            as.data.frame(
                t(
                    as.matrix(
                        c(
                            COUNT = sum(x),
                            colSums(x),
                            CLUST_WEIGHT = dim(x)[1]
                        ))))},
        mc.cores = cores)
    db_total <- do.call(rbind, db_total)
    db_total <- cbind(data.frame(ID = rownames(db_total)), db_total)
    db_total <- merge(db_total, db[, grepl("ID|sequence", names(db))], by = "ID")
    db$COUNT <- rowSums(db[, names(db) %in% samples])
    # order the columns
    col_order <- c("ID", "COUNT", "MOTU", samples,
                   "sequence")
    db <- db[, col_order]
    message("Finished calculating the number of reads in every sample for each MOTU.")
    return(list(db = db, db_total = db_total))
}


result <- SWARM2tab(fileswarm = opt$swarm, filetab = opt$tab, cores = opt$cores)
if (!is.null(opt$output_ESV)) {
  write.table(result$db_total, file = opt$output_ESV, sep = "\t", row.names = FALSE, quote = FALSE)
}
if (!is.null(opt$output_MOTU)) {
  write.table(result$db, file = opt$output_MOTU, sep = "\t", row.names = FALSE, quote = FALSE)
}