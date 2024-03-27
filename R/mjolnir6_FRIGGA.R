#' FRIGGA: Final Recount and Integration of Generated Genealogies and Abundances
#' 
#' FRIGGA recombined the abundances from ODIN with the taxonomic-annotated TSV files
#' from THOR
#' 
#' @details 
#' Input files must be called <EXPX>_THOR_annotated.tsv for the taxonomic-annotated
#' TSV file and <EXPX>_ODIN_counts.tsv for the read counts of MOTUs or <EXPX>_ODIN_counts.tsv
#' for ESVs.
#' Output file is then called <EXPX>_FRIGGA.tsv
#' 
#' @param experiment Character string. Acronym for the experiment. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. However they can be the same.
#' 
#' @export 
#' 
#' @examples 
#' library(mjolnir)
#'
#' # Define input fastq files (only names of R1 files are needed)
#' R1_filenames <- c("ULO1_R1.fastq.gz", "ULO2_R1.fastq.gz", "ULO3_R1.fastq.gz",
#'                   "ULO4_R1.fastq.gz")
#'
#' # Input identifiers for the individual libraries to be used. 
#' # It should be a 4-character name, matching the information in the 
#' # ngsfilter files.
#' lib_prefixes <- c("ULO1", "ULO2", "ULO3", "ULO4")
#'
#' # experiment identifier
#' experiment <- 'ULOY'
#' # Enter number of cores to be used in parallel.
#' cores <- 7
#'
#' mjolnir1_RAN(R1_filenames, lib_prefix = lib_prefixes, experiment = experiment,
#'              cores = cores, R1_motif = "_R1", R2_motif = "_R2")
#'
#' # Run FREYJA
#' mjolnir2_FREYJA(experiment = experiment, cores = cores, Lmin=299, Lmax=320)
#'
#' # Run HELA
#' mjolnir3_HELA(experiment = experiment, cores = cores)
#'
#' # Run ODIN
#' mjolnir4_ODIN(experiment = experiment, cores = cores, d = 13, 
#'               min_reads_MOTU = 2, min_reads_ESV = 2,
#'               min_relative = 1 / 50000, blank_relative = 0.1, 
#'               metadata_table = "", blank_col = "BLANK", blank_tag = TRUE, 
#'               alpha = 4, entropy = c(0.47, 0.23, 1.02, 313), 
#'               algorithm = "DnoisE_SWARM")
#' 
#' # set the directory where the database is stored
#' tax_db <- "~/taxo_NCBI/DUFA_COI"
#' 
#' # Run THOR
#' mjolnir5_THOR(experiment = experiment, cores = cores, 
#'               tax_db = tax_db, run_ecotag = T)
#' 
#' # Run FRIGGA
#' mjolnir6_FRIGGA(experiment = experiment)

mjolnir6_FRIGGA <- function(experiment=NULL, ...){
  
  if (exists("lib") && is.null(experiment)) {
    # Use lib as experiment
    experiment <- lib
    # Print deprecation warning
    warning("The 'lib' argument is deprecated. Please use 'experiment' instead.")
  }
  message("FRYGGA will produce a combined file.")

  infile=paste0(experiment,"_THOR_annotated.tsv")
  abundances=paste0(experiment,"_ODIN_counts.tsv")
  outfile=paste0(experiment,"_FRIGGA.tsv")
  if (!file.exists(abundances)) {
    abundances=paste0(experiment,"_ODIN_ESV.tsv")
  }

  message("FRYGGA is reading the ecotag-annotated database from THOR...")
  ecotag_db <- read.table(infile,sep="\t",head=T,stringsAsFactors=F)
  message("FRYGGA has read the taxonomy database, with ", nrow(ecotag_db)," total MOTUs.")
  # Delete "None" from the taxonomy database
  ecotag_db[ecotag_db=="None"] <- ""

  message("FRYGGA is reading the abundances database from ODIN...")
  abun_db <- read.table(abundances,sep="\t",head=T,stringsAsFactors=F)
  n_samples <- length(grep("sample",names(abun_db)))
  message("FRYGGA has read the abundances database, including ", nrow(abun_db)," total MOTUs and ",n_samples," samples.")

  # Merge databases
  names(abun_db)[names(abun_db)=="ID"] <- "id"
  names(abun_db)[names(abun_db)=="NUC_SEQ"] <- "sequence"

  db <- merge(ecotag_db,abun_db,by="id")

  names(db)[grep("sample",names(db))] <- gsub("MERGED_sample.","",names(db)[grep("sample",names(db))])
  db$COUNT <- rowSums(db[,grep("sample",names(db))])

  write.table(db,outfile,sep="\t",quote=F,row.names=F)
  message("FRYGGA is done. File ", outfile, " written, including ",nrow(db)," MOTUs with ",sum(db$COUNT)," total reads in ",n_samples," samples.")
  message("(",sum(db$COUNT>1)," non-singletons MOTUs).")
}

