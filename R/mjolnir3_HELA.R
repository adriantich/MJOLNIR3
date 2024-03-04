#' HELA: Hierarchical Elimination of Lurking Artifacts
#'
#' This function uses the uchime_denovo algorithm implemented in VSEARCH to
#' remove chimaeric sequences from the dataset.
#'
#' @details
#' HELA works in a sample-by-sample basis. HELA will process all
#' individual fasta files in the current folder
#' matching the pattern EXPX_XXXX_sample_XXX.fasta being EXPX the acronym set
#' by experiment parameter. This allows for parallel
#' computing, significantly decreasing calculation times.
#'
#' @param experiment Character string. Acronym for the experiment. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. However they can be the same.
#'
#' @param cores Numeric. Number of threads for parallel processing.
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


mjolnir3_HELA <- function(experiment = NULL, cores = 1, ...){

  suppressPackageStartupMessages(library(parallel))
  if (exists("lib") && is.null(experiment)) {
    # Use lib as experiment
    experiment <- lib
    # Print deprecation warning
    warning("The 'lib' argument is deprecated. Please use 'experiment' instead.")
  }
  sample_list <- gsub("_FREYJA_uniq.fasta", "",
                      list.files(pattern = paste0("^", experiment,
                                                  "_[a-zA-Z0-9]{4}_sample_[a-zA-Z0-9]{3}_FREYJA_uniq.fasta$")))

  message("HELA will remove chimaeras from each sample")
  X <- NULL
  for (i in sample_list) {
    X <- c(X, paste0("vsearch --uchime_denovo ", i, "_FREYJA_uniq.fasta ",
                     "--sizeout --minh 0.90 ",
                     "--nonchimeras ", i, "_HELA.fasta "))
  }
  mclapply(X, function(x) system(x, intern = TRUE, wait = TRUE),
           mc.cores = cores)

  after_HELA <- mclapply(sample_list,function(file){
    output <- system(paste0("grep '>' ",
                            file,"_HELA.fasta | wc -l"),
                     intern = TRUE, wait = TRUE)
    value <- as.numeric(output)
    return(data.frame(file=paste0(file,"_HELA.fasta"),
                      num_seqs=value))
  },mc.cores = cores)
  after_HELA <- do.call("rbind",after_HELA)
  report_HELA <- paste("HELA removed the chimeras with the uchime_denovo ",
                       "algorithm and kept for each sample the following",
                       "number of non-chimeras:\n",
                       "sample\tsequences\n",
                       paste(apply(after_HELA, 1, paste0, collapse = "\t"),
                             collapse = "\n"),
                       "\n")

  save(file = "summary_HELA.RData",list = c("after_HELA", "report_HELA"))

  message("HELA is done.")
}
