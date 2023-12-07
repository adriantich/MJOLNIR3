#' RAN: Reads Allotment in N samples
#'
#' RAN function will prepare the FASTQ raw data for every sample without sample
#' tags and primers. The R1 output files will contain all forward sequences and
#' R2 the reverse sequences.
#'
#' @details
#' If samples are already demultiplexed primers need to be set or LERAY_XT primers
#' for COI will be used by default. Also RAN will read the
#' names of each individual R1 fastq files from a column in the LIBX_metadata.tsv
#' file, called fastq_name_R1. In the metadata table, each sample in the
#' original_samples column must have a their corresponding fastq_name_R1 and
#' mjolnir_agnomen (LIBX_sample_XXX).
#' If the fastq_name_R1. column is not in the metadata, RAN will use the lib_prefix.
#' Files required:
#' - ngsfilter file, needed only for multiplexed libraries
#'      For each library, a ngsfilter file is needed and must be named
#'      ngsfiler_<library identifier>.tsv. This must contain five tab-separated
#'      columns and no header. The first column with the library identifier, the
#'      second with the mjolnir_agnomens, the third with the sample tags, the
#'      fourth with the forward primers and the fifth with the reverse primers.
#' - metadata file
#'      a metadata file containing at least two columns required,
#'      'original_samples' and 'mjolnir_agnomens' (fastq_name_R1 for demultiplexed
#'       libraries), and named as <experiment identifier>_metadata.tsv.
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
#' @param primer_F Character string of the Forward primer. Necessary when
#' samples are already demultiplexed.
#'
#' @param primer_R Character string of the Reverse primer. Necessary when
#' samples are already demultiplexed.
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
#' @export 
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

mjolnir1_RAN <- function(R1_filenames = "", lib_prefix = "",
                         experiment = NULL, lib = NULL,
                         primer_F="GGWACWRGWTGRACWNTNTAYCCYCC",primer_R="TANACYTCNGGRTGNCCRAARAAYCA",
                         cores = 1, R1_motif = "_R1", R2_motif = "_R2") {

  if (!is.null(lib) && is.null(experiment)) {
    # Use lib as experiment
    experiment <- lib
    # Print deprecation warning
    warning("The 'lib' argument is deprecated. Please use 'experiment' instead.")
  }
  if (is.null(experiment)){
    stop("Error: experiment argument required")
  }
  metadata <- read.table(paste0(experiment, "_metadata.tsv"),
                         sep = "\t", header = TRUE)
  if ('fastq_name_R1' %in% names(metadata)) {
    fastqR1_list <- metadata$fastq_name_R1
    agnomens <-  metadata$original_samples
    for (j in seq_len(length(fastqR1_list))) {
      R1_file <- fastqR1_list[j]
      R2_file <- gsub(R1_motif, R2_motif, R1_file)
      fwd_outfile <- paste0(agnomens[j],
                            R1_motif,
                            ".fastq")
      rev_outfile <- paste0(agnomens[j],
                            R2_motif,
                            ".fastq")
      if (R1_file == fwd_outfile) {
        stop(paste("error: fastq_name_R1 and original_samples are the same. ",
                   "Please consider changing one of them.",
                   "The ouput would be named as the input and in Valhalla we don't do that.",
                   "You can use .fastq.gz files as input instead."))
      }
      cutadapt_command_primers <- paste0(cutadapt_command_primers,
                                        "cutadapt -e 0.1 ", # allow 0.1 errors for primers in ngsfilter was total of 2
                                        "-O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required # nolint: line_length_linter.
                                        " --no-indels ", # no indels allowed # nolint: line_length_linter.
                                        "-j ", cores, # number of cores allowed # nolint: line_length_linter.
                                        " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                                        " --max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja # nolint: line_length_linter.
                                        "-g ^", primer_F, " -G ^", primer_R, # these are the sample_tags # nolint: line_length_linter.
                                        " -o ", fwd_outfile, " -p ", rev_outfile, # sample names as original # nolint: line_length_linter.
                                        " ", R1_file, " ", R2_file, " ; ") # input files # nolint: line_length_linter.
      if (primer_F != primer_R) {
        cutadapt_command_primers <- paste0(cutadapt_command_primers,
                                          "cutadapt -e 0.1 ", # allow 0.1 errors for primers in ngsfilter was total of 2
                                          "-O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required # nolint: line_length_linter.
                                          " --no-indels ", # no indels allowed # nolint: line_length_linter.
                                          "-j ", cores, # number of cores allowed # nolint: line_length_linter.
                                          " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                                          " --max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja # nolint: line_length_linter.
                                          "-g ^", primer_R, " -G ^", primer_F, # these are the sample_tags # nolint: line_length_linter.
                                          " -o temp_fileR1_rev.fastq -p temp_fileR2_fwd.fastq", # sample names as original # nolint: line_length_linter.
                                          " ", R1_file, " ", R2_file, " ; ",
                                          # now concatenate and remove
                                          "cat temp_fileR1_rev.fastq >>", rev_outfile, " ; rm temp_fileR1_rev.fastq ; ",
                                          "cat temp_fileR2_fwd.fastq >>", fwd_outfile, " ; rm temp_fileR2_fwd.fastq ; ")
      }
      system(cutadapt_command_primers, wait = TRUE, intern = TRUE)
    }
  } else {

    message("RAN will demultiplex the following initial FASTQ files:")
    message(paste(R1_filenames))
    ngsfilter_files <- list.files('.', 'ngsfilter')
    if (length(ngsfilter_files) == 0) {
      stop("No ngsfilter file found.")
    }
    for (i in seq_along(R1_filenames)) {
      lib <- lib_prefix[i]
      R1_file <- R1_filenames[i]
      R2_file <- gsub(R1_motif, R2_motif, R1_file)
      loop <- 1
      R1_file_temp <- paste0(lib, 'untrimmed', loop, R1_motif, '.fastq')
      R2_file_temp <- paste0(lib, 'untrimmed', loop, R2_motif, '.fastq')
      message(paste0("RAN is processing the ", lib, " library."))
      ngsfile <- read.csv(ngsfilter_files[grep(lib, ngsfilter_files)],
                          sep = "\t", header = FALSE)
      message(paste0("RAN will split initial ", R1_file, " & ",
                    gsub(R1_motif, R2_motif, R1_file),
                    " files in ", dim(ngsfile)[1], " samples."))
      if (sum(duplicated(paste0(ngsfile$V3,ngsfile$V4,ngsfile$V5))) > 0) {
        stop(paste("There are duplicated sample tags",
                  "with the same primers in the ",
                  ngsfilter_files[grep(lib, ngsfilter_files)], " file."))
      }
      for (tag_combi in unique(ngsfile$V3)) {
        sample_num <- which(ngsfile$V3 == tag_combi)

        # for each loop every tag will be of length one
        fwd_tag <- gsub(":.*", "", tag_combi)
        rev_tag <- gsub(".*:", "", tag_combi)


        # these can be multivalue in case that the sem combi is found multiple times
        fwd_outfile <- paste0(metadata$original_samples[metadata$mjolnir_agnomens %in% ngsfile$V2[sample_num]], # nolint: line_length_linter.
                              R1_motif,
                              ".fastq")
        rev_outfile <- gsub(R1_motif, R2_motif, fwd_outfile)
        fwd_primer <- ngsfile$V4[sample_num]
        rev_primer <- ngsfile$V5[sample_num]
        temp_file_R1 <- paste0('tag1_temp01.fastq')
        temp_file_R2 <- paste0('tag2_temp01.fastq')

        cutadapt_command_tag <- c()
        cutadapt_command_primers <- c()

        cutadapt_command_tag <- paste0(cutadapt_command_tag,
                                        "cutadapt -e 0 ", # allow 0 errors
                                        "-O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required # nolint: line_length_linter.
                                        " --no-indels ", # no indels allowed # nolint: line_length_linter.
                                        "-j ", cores, # number of cores allowed # nolint: line_length_linter.
                                        # " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                                        " --untrimmed-output ", R1_file_temp, # save those reads that have not been assigned to the sample # nolint: line_length_linter.
                                        " --untrimmed-paired-output ", R2_file_temp, # save those reads that have not been assigned to the sample # nolint: line_length_linter.
                                        "--max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja # nolint: line_length_linter.
                                        "-g ", fwd_tag, " -G ", rev_tag, # these are the sample_tags # nolint: line_length_linter.
                                        " -o ", temp_file_R1, " -p ", temp_file_R2, # sample names as original # nolint: line_length_linter.
                                        " ", R1_file, " ", R2_file, " ; ") # input files # nolint: line_length_linter.
        if (loop > 1) {
          cutadapt_command_tag <- paste0(cutadapt_command_tag,
                                          " rm ", R1_file, " ; rm ", R2_file, " ; ")
        }
        R1_file <- R1_file_temp
        R2_file <- R2_file_temp
        loop <- loop + 1
        R1_file_temp <- paste0(lib, 'untrimmed', loop, R1_motif, '.fastq')
        R2_file_temp <- paste0(lib, 'untrimmed', loop, R2_motif, '.fastq')
        if (fwd_tag != rev_tag) {
          # if this is the case then is a little bit more complex
          # we need to run cutadapt twice changing the order of sample_tags and then concatenate files
          cutadapt_command_tag <- paste0(cutadapt_command_tag, # this is the first command
                                          # the second command change the order of sample_tags
                                          "cutadapt -e 0 ", # allow 0 errors
                                          "-O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required
                                          " --no-indels ", # no indels allowed
                                          "-j ", cores, # number of cores allowed
                                          # " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                                          " --untrimmed-output ", R1_file_temp, # save those reads that have not been assigned to the sample # nolint: line_length_linter.
                                          " --untrimmed-paired-output ", R2_file_temp, # save those reads that have not been assigned to the sample # nolint: line_length_linter.
                                          "--max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja
                                          "-g ", rev_tag, " -G ", fwd_tag, # these are the sample_tags
                                          " -o temp_fileR1.fastq -p temp_fileR2.fastq", # sample names as original
                                          " ", R1_file, " ", R2_file, " ; ",
                                          # now concatenate and remove
                                          "cat temp_fileR1.fastq >>", temp_file_R1, " ; rm temp_file1.fastq ; ",
                                          "cat temp_fileR2.fastq >>", temp_file_R2, " ; rm temp_file2.fastq ; ")
          cutadapt_command_tag <- paste0(cutadapt_command_tag,
                                          " rm ", R1_file, " ; rm ", R2_file, " ; ")
          R1_file <- R1_file_temp
          R2_file <- R2_file_temp
          loop <- loop + 1
          R1_file_temp <- paste0(lib, 'untrimmed', loop, R1_motif, '.fastq')
          R2_file_temp <- paste0(lib, 'untrimmed', loop, R2_motif, '.fastq')

        }

        for (j in seq_len(length(sample_num))) {
          cutadapt_command_primers <- paste0(cutadapt_command_primers,
                                            "cutadapt -e 0.1 ", # allow 0.1 errors for primers in ngsfilter was total of 2
                                            "-O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required # nolint: line_length_linter.
                                            " --no-indels ", # no indels allowed # nolint: line_length_linter.
                                            "-j ", cores, # number of cores allowed # nolint: line_length_linter.
                                            " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                                            " --max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja # nolint: line_length_linter.
                                            "-g ^", fwd_primer[j], " -G ^", rev_primer[j], # these are the sample_tags # nolint: line_length_linter.
                                            " -o ", fwd_outfile[j], " -p ", rev_outfile[j], # sample names as original # nolint: line_length_linter.
                                            " ", temp_file_R1, " ", temp_file_R2, " ; ") # input files # nolint: line_length_linter.
          if (fwd_tag != rev_tag) {
            cutadapt_command_primers <- paste0(cutadapt_command_primers,
                                              "cutadapt -e 0.1 ", # allow 0.1 errors for primers in ngsfilter was total of 2
                                              "-O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required # nolint: line_length_linter.
                                              " --no-indels ", # no indels allowed # nolint: line_length_linter.
                                              "-j ", cores, # number of cores allowed # nolint: line_length_linter.
                                              " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                                              " --max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja # nolint: line_length_linter.
                                              "-g ^", rev_primer[j], " -G ^", fwd_primer[j], # these are the sample_tags # nolint: line_length_linter.
                                              " -o temp_fileR1_rev.fastq -p temp_fileR2_fwd.fastq", # sample names as original # nolint: line_length_linter.
                                              " ", temp_file_R1, " ", temp_file_R2, " ; ",
                                              # now concatenate and remove
                                              "cat temp_fileR1_rev.fastq >>", rev_outfile[j], " ; rm temp_fileR1_rev.fastq ; ",
                                              "cat temp_fileR2_fwd.fastq >>", fwd_outfile[j], " ; rm temp_fileR2_fwd.fastq ; ")

          }
        }
        cutadapt_command_primers <- paste0(cutadapt_command_primers,
                                          "rm ", temp_file_R1, " ; rm ", temp_file_R2, " ; ")
        system(cutadapt_command_tag, wait = TRUE, intern = TRUE)
        system(cutadapt_command_primers, wait = TRUE, intern = TRUE)
      }
    }
  }

  message("Demultiplexing done.")
}
