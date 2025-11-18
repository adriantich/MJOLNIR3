#' FREYJA: Filtering of Reads, Enrollment, Yoke-reads Joining and Alignment
#'
#' FREYJA will use OBITools3 commands to merge paired-end reads, trim primer
#' sequences, filter by length, split sequences per sample and dereplicate
#' within each sample.
#'
#' @details
#' Input file fastq files are expected to be without primers sequence and all
#' forward sequences in the R1 file and all reverse sequences in the R2 file.
#'
#' @param experiment Character string. Acronym for the experiment. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. However they can be the same.
#'
#' @param cores Numeric. Number of threads for parallel processing.
#'
#' @param Lmin Numeric. Minimum bp length for a sequence to be accepted.
#'
#' @param Lmax Numeric. Maximum bp length for a sequence to be accepted.
#'
#' @param min_overlap Numeric. Minimum overlap length between R1 and R2 reads.
#' 
#' @param maxdiff Numeric. Maximum number of differences allowed in the overlap.
#' If less than 1, it will be considered as a percentage of the overlap length.
#' 
#' @param error_rate Numeric. Maximum expected error rate allowed for the output
#' merged reads. For a given sequence, the expected error is the sum of error
#' probabilities for all the positions in the sequence. 
#' Since error probabilities can be small but not null, the expected error is
#' always greater than zero, and at most equal to the length of the sequence
#' when all positions in the sequence have an error probability of 1.0.
#'
#' @param R1_motif Character string that distinguish the forward line file from
#' the reverse.
#'
#' @param R2_motif Character string that distinguish the reverse line file from
#' the forward.
#'
#' @param commands_file Character string. Name of the file where all commands
#' will be written. If NULL or missing, commands will not be recorded.
#' 
#' @param only_commands Logical. If TRUE, only the commands will be written to
#' the commands_file. If FALSE, the commands will be executed.
#' 
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

mjolnir2_FREYJA <- function(experiment = NULL, cores = 1, Lmin = 299, Lmax = 320,
                            min_overlap = 40,
                            maxdiff = 0,
                            error_rate = 0.05,
                            R1_motif = "_R1", R2_motif = "_R2",
                            commands_file = "commands_runned_FREYJA.txt",
                            only_commands = FALSE,
                            ...) {
  
  if (exists("lib") && is.null(experiment)) {
    # Use lib as experiment
    experiment <- lib
    # Print deprecation warning
    warning("The 'lib' argument is deprecated. Please use 'experiment' instead.")
  }

  # message("FREYJA will do paired-end alignment, demultiplexing and length filter.")
  message("FREYJA will do paired-end alignment and length filter.")
  suppressPackageStartupMessages(library(parallel))

  if (is.null(experiment)) {
    message("experiment can not be NULL, otherwise HELA won't find the files")
    stop()
  }

  # run all commands
  if (!is.logical(only_commands)){
    stop("Error: only_commands must be TRUE or FALSE")
  }
  if (commands_file !=  "" | !is.null(commands_file) | 
      !missing(commands_file)) {
    record_commands <- TRUE
    commands_file <- commands_file
  } else if (only_commands) {
    record_commands <- TRUE
    commands_file <- "commands_runned_FREYJA.txt"
    message(paste0("commands_file was not specified, so commands will be written to ",
                   commands_file))
  } else {
    record_commands <- FALSE
  }

  filtering_commands <- NULL
  message("FREYJA will first clear the battle field.")
  message("Any directory or file containing the word FREYJA will be removed.")
  # system("rm -r *FREYJA*", intern = TRUE, wait = TRUE)

  metadata <- read.table(paste0(experiment, "_metadata.tsv"),
                         sep = "\t", header = TRUE)

  fastqR1_list <- paste0(metadata$original_samples,
                         R1_motif,
                         ".fastq")
  fastqR2_list <- paste0(metadata$original_samples,
                         R2_motif,
                         ".fastq")

  fastqR1_list <- gsub("..fastq", ".fastq", fastqR1_list, fixed = TRUE)
  fastqR2_list <- gsub("..fastq", ".fastq", fastqR2_list, fixed = TRUE)
      
  agnomens <-  metadata$mjolnir_agnomens
  before_FREYJA <- mclapply(fastqR1_list, function(prefix){
    return(data.frame(file=prefix,
                      num_seqs=as.numeric(gsub(' .*','',system(paste0("wc -l ",prefix," "),intern = T,wait = T)))/4))
  }, mc.cores = cores)

  to_retain <- do.call("rbind",before_FREYJA)
  fastqR1_list <- fastqR1_list[to_retain$num_seqs > 0]
  fastqR2_list <- fastqR2_list[to_retain$num_seqs > 0]
  agnomens <- agnomens[to_retain$num_seqs > 0]
  if (maxdiff < 1) {
    maxdiffs_param <- " --fastq_maxdiffpct "
  } else {
    maxdiffs_param <- " --fastq_maxdiffs "
  }
  if (!exists("additional_params_alignment")) {
    additional_params_alignment <- ""
  }
  for (i in seq_along(agnomens)) {
    print(fastqR1_list[i])
    filtering_commands <-
      c(filtering_commands,
        paste0(
          "vsearch ",
          " --fastq_mergepairs ", fastqR1_list[i],
          " --reverse ", fastqR2_list[i],
          " --fastqout ", experiment, "_", agnomens[i], "_FREYJA_aligned.fastq",
          " --fastq_minovlen ", min_overlap,
          maxdiffs_param, maxdiff,
          " --fastq_maxee ", error_rate,
          " --fastq_minmergelen ", Lmin,
          " --fastq_maxmergelen ", Lmax,
          " --fastq_maxns 0 ",
          " ", additional_params_alignment, " ",
          " ; ",
          " vsearch ",
          " --fastx_uniques ", experiment, "_", agnomens[i], "_FREYJA_aligned.fastq",
          " --sizeout ",
          " --fastaout ", experiment, "_", agnomens[i], "_FREYJA_uniq.fasta"
        )
      )
  }
  if (only_commands) {
    writeLines(filtering_commands, con = commands_file)
    message(paste0("Only commands were written to ", commands_file))
    return(invisible(NULL))
  } else {
    if (record_commands) {
      writeLines(filtering_commands, con = commands_file)
      message(paste0("Commands were written to ", commands_file))
    }
  }
  mclapply(filtering_commands, function(x) system(x,intern=T,wait=T), mc.cores = cores)

  # obi uniq vas performed in HELA in previous versions but now is computed here
  files <- list.dirs(recursive = F)
  files <- files[grepl("sample",files)&grepl("FREYJA_aligned.fastq",files)]
  files <- gsub("./","",gsub("_aligned.fastq","",files))

  after_FREYJA <- mclapply(files,function(file){
    file <- "ULOY_ULO1_sample_001_FREYJA"
    output <- system(paste0("wc -l ", file, "_aligned.fastq"),
                     intern = TRUE, wait = TRUE)
    sequences <- as.numeric(gsub(paste0(" ",file,"_aligned.fastq"), "", output)) / 4
    output <- system(paste0("grep '>' ", file, "_uniq.fasta | wc -l"),
                     intern = TRUE, wait = TRUE)
    uniq_seqs <- as.numeric(output)

    return(data.frame(file=file,
                      version=c("filtered sequences","uniq sequences"),
                      num_seqs=c(sequences, uniq_seqs)))
  },mc.cores = cores)

  variables_FREYJA <- data.frame(variable = c("cores", "Lmin",
                                              "Lmax", "experiment",
                                              "score_obialign"),
                                 value = c(cores, Lmin,
                                           Lmax, experiment,
                                           score_obialign))

  save(file = "summary_FREYJA.RData",
       list = c("before_FREYJA", "after_FREYJA", "variables_FREYJA"))

  message("FREYJA is done.")
}
