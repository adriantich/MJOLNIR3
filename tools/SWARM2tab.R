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

# fileswarm="swarm_output.txt"
# filetab="seq_table.tsv"
# cores=1

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
    
    swarm_db <- gsub("size=[0-9]+\\s*","",swarm_db)  # remove leading numbers
    
    clusters <- strsplit(swarm_db, ";")
    message("4/7. Keeping only information of the sequences ",
            "that form each cluster.")
    # clusters <- mclapply(clusters,function(x) {sub(";.*", "", x)},
    #                      mc.cores = cores)
    names(clusters) <- mclapply(clusters, function(x) x[[1]], mc.cores = cores)
    
    # Read counts database and keep only the needed clusters
    message("6/7. Reading the abundance database. ",
            "This could take a while.")
    motu_seqs_names <- stack(clusters) %>% rename(ID = values, MOTU = ind)
    
    # db_ESV <- read.table(filetab, sep = "\t", head = TRUE)
    db_ESV <- read.delim(filetab, header=TRUE, sep="\t",
                    comment.char="", check.names=FALSE,
                    stringsAsFactors=FALSE)
    names(db_ESV) <- sub("^#OTU ","", names(db_ESV))     # remove leading '#'
    numseqs <- nrow(db_ESV)
    samples <- names(db_ESV)[!names(db_ESV)%in%c("ID","sequence")]
    db_ESV <- merge(motu_seqs_names, db_ESV, by = "ID")
    numseqs_reduced <- nrow(db_ESV)
    num_samples <- length(samples)
    message("Finished reading the Database, which includes ", 
            numseqs, " total unique sequences and ", num_samples, " samples.\n",
            "Kept only ", numseqs_reduced, " sequences for calculations.")
    
    message("7/7. Calculating the number of reads in every sample ",
            "for each MOTU.")
    db_MOTU <- split(db_ESV[, names(db_ESV) %in% samples], db_ESV$MOTU)
    db_MOTU <- mclapply(
        db_MOTU, 
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
    db_MOTU <- do.call(rbind, db_MOTU)
    db_MOTU <- cbind(data.frame(ID = rownames(db_MOTU)), db_MOTU)
    db_MOTU <- merge(db_MOTU, db_ESV[, grepl("ID|sequence", names(db_ESV)), drop = FALSE], by = "ID")
    db_ESV$COUNT <- rowSums(db_ESV[, names(db_ESV) %in% samples])
    # order the columns
    col_order <- c("ID", "COUNT", "MOTU", samples)
    db_ESV <- db_ESV[, col_order]
    db_MOTU <- db_MOTU[order(db_MOTU$COUNT, decreasing = TRUE), ]
    message("Finished calculating the number of reads in every sample for each MOTU.")
    return(list(db_ESV = db_ESV, db_MOTU = db_MOTU))
}


result <- SWARM2tab(fileswarm = opt$swarm, filetab = opt$tab, cores = opt$cores)
if (!is.null(opt$output_ESV)) {
  write.table(result$db_ESV, file = opt$output_ESV, sep = "\t", row.names = FALSE, quote = FALSE)
}
if (!is.null(opt$output_MOTU)) {
  write.table(result$db_MOTU, file = opt$output_MOTU, sep = "\t", row.names = FALSE, quote = FALSE)
}