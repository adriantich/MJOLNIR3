#' RAN: Reads Allotment in N samples
#'
#' RAN function will prepare the FASTQ raw data for parallel processing.
#'
#' @details
#' Use RAN if your data consist of multiplexed libraries. These will be split
#' samples to be processed by FREYJA.
#' - ngsfilter file
#'      For each library, a ngsfilter file is needed and must be named
#'      ngsfiler_<library identifier>.tsv. This must contain at least three
#'      columns and no header. The first column with the library identifier, the
#'      second with the mjolnir_agnomens and the third with the sample tags
#' - metadata file
#'      a metadata file containing at least two columns required,
#'      'original_samples' and 'mjolnir_agnomensis', and named as
#'      <experiment identifier>_metadata.tsv.
#'
#' @param R1_filenames Character vector with the names of the forward fastq or
#' fastq.gz files.
#'
#' @param lib_prefix Character vector. Acronym for each sequencing library. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. The latter will be required in following steps. However
#' they can be the same.
#'
#' @param experiment Character string. Acronym for the experiment. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. However they can be the same.
#'
#' @param cores Numeric. Number of parts into which the input files will be split
#' for the parallel processing of the FREYJA function.
#'
#' @param R1_motif Character string that distinguishes the forward line file from
#' the reverse.
#'
#' @param R2_motif Character string that distinguishes the reverse line file from
#' the forward.
#'
#' @examples
#' library(mjolnir)
#'
#' # Define input fastq files (only names of R1 files are needed)
#' R1_filenames <-c("ULO1_R1.fastq.gz","ULO2_R1.fastq.gz","ULO3_R1.fastq.gz","ULO4_R1.fastq.gz")
#'
#' # Input identifiers for the individual libraries to be used. It should be a 4-character name, matching the information in the ngsfilter files
#' lib_prefixes <- c("ULO1","ULO2","ULO3","ULO4")
#'
#' # experiment identifier
#' experiment <- 'ULOY'
#' # Enter number of cores to be used in parallel.
#' cores <- 7
#'
#' mjolnir1_RAN(R1_filenames, lib_prefix = lib_prefixes, experiment = experiment,
#'              cores = cores, R1_motif = "_R1", R2_motif = "_R2")

mjolnir1_RAN <- function(R1_filenames, lib_prefix, experiment = NULL, lib = NULL,
                         cores = 1, R1_motif = "_R1", R2_motif = "_R2") {

  if (!is.null(lib) && is.null(experiment)) {
    # Use lib as experiment
    experiment <- lib
    # Print deprecation warning
    warning("The 'lib' argument is deprecated. Please use 'experiment' instead.")
  }

  message("RAN will demultiplex the following initial FASTQ files:")
  message(paste(R1_filenames))
  ngsfilter_files <- list.files('.', 'ngsfilter')
  metadata <- read.table(paste0(experiment, "_metadata.tsv"), sep = "\t", header = TRUE)
  if (length(ngsfilter_files) == 0) {
    stop("No ngsfilter file found.")
  }
  for (i in seq_along(R1_filenames)) {
    lib <- lib_prefix[i]
    R1_file <- R1_filenames[i]
    message(paste0("RAN is processing the ", lib, " library."))
    ngsfile <- read.csv(ngsfilter_files[grep(lib, ngsfilter_files)], sep = "\t", header = FALSE)
    message(paste0("RAN will split initial ", R1_file, " & ", gsub(R1_motif, R2_motif, R1_file),
                   " files in ", dim(ngsfile)[1], " samples."))
    for (sample_num in seq_len(dim(ngsfile)[1])) {
      fwd_tag <- gsub(':.*', '', ngsfile$V3[sample_num])
      rev_tag <- gsub('.*:', '', ngsfile$V3[sample_num])
      fwd_outfile <- paste0(metadata$original_samples[metadata$mjolnir_agnomens == ngsfile$V2[sample_num]], # nolint: line_length_linter.
                            R1_motif,
                            ".fastq")
      cutadapt_command <- paste0("cutadapt -e 0 ", # allow 0 errors
                                 "-O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required
                                 " --no-indels ", # no indels allowed
                                 "-j ", cores, # number of cores allowed
                                 " --discard-untrimmed ", # discard those reads that have not been assigned to the sample
                                 "--max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja
                                 "-g ", fwd_tag, " -G ", rev_tag, # these are the sample_tags
                                 " -o ", fwd_outfile, " -p ", gsub(R1_motif, R2_motif, fwd_outfile), # sample names as original
                                 " ", R1_file, " ", gsub(R1_motif, R2_motif, R1_file)) # input files
      if (fwd_tag != rev_tag) {
        # if this is the case then is a little bit more complex
        # we need to run cutadapt twice changing the order of sample_tags and then concatenate files
        cutadapt_command <- paste0(cutadapt_command, " ; ", # this is the first command
                                   # the second command change the order of sample_tags
                                   paste0("cutadapt -e 0 ", # allow 0 errors
                                          "-O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required
                                          " --no-indels ", # no indels allowed
                                          "-j ", cores, # number of cores allowed
                                          " --discard-untrimmed ", # discard those reads that have not been assigned to the sample
                                          "--max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja
                                          "-g ^$", rev_tag, " -G ^$", fwd_tag, # these are the sample_tags
                                          " -o temp_fileR1.fastq -p temp_fileR2.fastq", # sample names as original
                                          " ", R1_file, " ", gsub(R1_motif, R2_motif, R1_file)), " ; ",
                                   # now concatenate and remove
                                   "cat temp_fileR1.fastq >>", fwd_outfile, " ; rm temp_file1.fastq ; ",
                                   "cat temp_fileR2.fastq >>", gsub(R1_motif,R2_motif,fwd_outfile), " ; rm temp_file2.fastq ; ")

      }
      system(cutadapt_command, wait = TRUE, intern = TRUE)
    }
  }
  message("Demultiplexing done.")
}
