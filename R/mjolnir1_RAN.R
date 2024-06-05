#' RAN: Reads Allotment in N samples
#'
#' RAN function will prepare the FASTQ raw data for each sample with the sample
#' tags and primers removed. The R1 output files will contain all forward sequences and
#' R2 the reverse sequences. Please, read the Details section before running.
#'
#' @details
#' RAN considers the following scenarios:
#' 
#' - \bold{Scenario 1}. The samples are multiplexed in a library/libraries 
#'    (it is possible to run multiple libraries at once from a single experiment).
#'    Each library has only two raw .fastq(.gz), R1 and R2.
#'    For each library there must be only one ngsfile (see below)
#'    and for the whole experiment only one metadata file (see below).
#'    mjolnir_agnomens (the standard name of the sample during the pipeline)
#'    must be unique, it must not be the same in different ngsfiles.
#'    \emph{lib_prefix} and \emph{R1_filenames} must be the same length.
#' 
#' - \bold{Scenario 2}. The samples are multiplexed but for each library we have more
#'    than two raw data files. For example, when your library has been
#'    sequenced across multiple lanes on a Novaseq. In this case, RAN has
#'    to be run separately for each library and the option
#'    \emph{multilane} must be set to TRUE. In this case, \emph{lib_prefix} can only be
#'    of length one (only one library), but the \emph{R1_filenames} can be longer.
#' 
#' - \bold{Scenario 3}. The samples are demultiplexed but the primer remains.
#'    In this case RAN will trimm the primers and separate the sequences
#'    into fwd and rev files.
#' 
#' 
#' \bold{Starting with multiplexed libraries:}
#' 
#' Files required:
#' 
#' - \emph{ngsfilter file}. Needed only for multiplexed libraries
#'      For each library, a ngsfilter file is needed and must be named
#'      \bold{ngsfilter_<library identifier>.tsv}. This must contain five tab-separated
#'      columns and no header. The first column with the library identifier
#'      (four charachter identifier), the second with the mjolnir_agnomens,
#'      the third with the sample tags, the
#'      fourth with the forward primers and the fifth with the reverse primers.
#' 
#' - \emph{metadata file}.
#'      A metadata file containing at least two columns required,
#'      "original_samples" and "mjolnir_agnomens" ("fastq_name_R1" for
#'      demultiplexed libraries, see below), and named as
#'      \bold{<experiment identifier>_metadata.tsv}.
#' 
#' Important: when the same library has different sequencing lanes (NovaSeq),
#' RAN has to be run separately for each library and the option
#' \emph{multilane} must be set to TRUE.
#' 
#' 
#' \bold{Starting with demultiplexed samples:}
#' 
#' If samples are already demultiplexed, primers need to be set
#' (\emph{primer_F} & \emph{primer_R}) or LERAY_XT primers
#' for COI will be used by default.
#' 
#' Files required:
#' 
#' - \emph{metadata file}. RAN will read the
#'      names of each individual R1 fastq files (full name including extension)
#'      from a column in the \bold{<experiment identifier>_metadata.tsv}
#'      file, called "fastq_name_R1". In the metadata table, each sample in the
#'      "original_samples" column must have a their corresponding fastq_name_R1 and
#'      "mjolnir_agnomen" (LIBX_sample_XXX, i.e LIBA_sample_001).
#'      If the "fastq_name_R1" column is not in the metadata, RAN will use the
#'      \emph{lib_prefix}.
#' 
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
#' @param multilane Logical. If FALSE, the function will consider that the each library
#' has only one sequencing lane. If TRUE, only one library is processed at a time.
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

mjolnir1_RAN <- function(R1_filenames = "", lib_prefix = "",
                         experiment = NULL, 
                         primer_F="GGWACWRGWTGRACWNTNTAYCCYCC",primer_R="TANACYTCNGGRTGNCCRAARAAYCA",
                         cores = 1, R1_motif = "_R1", R2_motif = "_R2",
                         multilane = FALSE, ...) {

  if (exists("lib") && is.null(experiment)) {
    # Use lib as experiment
    experiment <- lib
    # Print deprecation warning
    warning("The 'lib' argument is deprecated. Please use 'experiment' instead.")
  }
  if (is.null(experiment)){
    stop("Error: experiment argument required")
  }
  if (exists("lib_prefixes") &&
      length(lib_prefix) == 1 && lib_prefix[1] == "") {
    # Use lib as experiment
    lib_prefix <- lib_prefixes
    # Print deprecation warning
    warning(paste("The 'lib_prefixes' argument is deprecated. Please use 'lib_prefix' instead.\n",
                  "It's okay by now but don't repeat it, okay?"))
  }
  if  (!multilane && length(lib_prefix) != length(R1_filenames)){
    stop(paste("Error: you have set the multilane to False but there",
               "are different number of R1_files and lib_prefix"))
  }
  metadata <- read.table(paste0(experiment, "_metadata.tsv"),
                         sep = "\t", header = TRUE)
  if ('fastq_name_R1' %in% names(metadata)) {
    # in this case the names for the fastq files are in the metadata which mean 
    # that the samples are demultiplexed but the primer needs to be removed
    fastqR1_list <- metadata$fastq_name_R1
    agnomens <-  metadata$original_samples

    cutadapt_command_primers <- c()
    for (j in seq_len(length(fastqR1_list))) {
      R1_file <- list.files(pattern = paste0("^", fastqR1_list[j]))
      R2_file <- gsub(R1_motif, R2_motif, R1_file)
      fwd_outfile <- paste0(agnomens[j],
                            R1_motif,
                            ".fastq")
      rev_outfile <- paste0(agnomens[j],
                            R2_motif,
                            ".fastq")
      fwd_outfile <- gsub("..fastq", ".fastq", fwd_outfile, fixed = TRUE)
      rev_outfile <- gsub("..fastq", ".fastq", rev_outfile, fixed = TRUE)
      if (R1_file == fwd_outfile) {
        stop(paste("error: fastq_name_R1 ->", R1_file, " and ",
                   "original_samples ->", fwd_outfile, " are the same. ",
                   "Please consider changing one of them.",
                   "The ouput would be named as the input and in Valhalla we don't do that.",
                   "You can use .fastq.gz files as input instead."))
      }
      cutadapt_command_primers <- paste0(cutadapt_command_primers,
                                        "cutadapt -e 0.1 ", # allow 0.1 errors for primers in ngsfilter was total of 2
                                        # "-O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required # nolint: line_length_linter.
                                        " --no-indels ", # no indels allowed # nolint: line_length_linter.
                                        "-j ", cores, # number of cores allowed # nolint: line_length_linter.
                                        " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                                        " --max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja # nolint: line_length_linter.
                                        "-g ", primer_F, " -G ", primer_R, # these are the sample_tags # for demultiplexed not anchored # nolint: line_length_linter.
                                        " -o ", fwd_outfile, " -p ", rev_outfile, # sample names as original # nolint: line_length_linter.
                                        " ", R1_file, " ", R2_file, " ; ") # input files # nolint: line_length_linter.
      if (primer_F != primer_R) {
        cutadapt_command_primers <- paste0(cutadapt_command_primers,
                                          "cutadapt -e 0.1 ", # allow 0.1 errors for primers in ngsfilter was total of 2
                                          # "-O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required # nolint: line_length_linter.
                                          " --no-indels ", # no indels allowed # nolint: line_length_linter.
                                          "-j ", cores, # number of cores allowed # nolint: line_length_linter.
                                          " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                                          " --max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja # nolint: line_length_linter.
                                          "-g ", primer_R, " -G ", primer_F, # these are the sample_tags # nolint: line_length_linter.
                                          " -o temp_fileR1_rev.fastq -p temp_fileR2_fwd.fastq", # sample names as original # nolint: line_length_linter.
                                          " ", R1_file, " ", R2_file, " ; ",
                                          # now concatenate and remove
                                          "cat temp_fileR1_rev.fastq >>", rev_outfile, " ; rm temp_fileR1_rev.fastq ; ",
                                          "cat temp_fileR2_fwd.fastq >>", fwd_outfile, " ; rm temp_fileR2_fwd.fastq ; ")
      }
      system(cutadapt_command_primers, wait = TRUE, intern = TRUE)
    }
  } else {
    # in this case the samples are not demultiplexed. 
    # first the samples are demultiplexed with the sample tags 
    # then the primers are removed from the demultiplexed samples.

    message("RAN will demultiplex the following initial FASTQ files:")
    message(paste(R1_filenames))
    ngsfilter_files <- list.files('.', 'ngsfilter')
    if (length(ngsfilter_files) == 0) {
      stop("No ngsfilter file found.")
    }
    if (multilane) {
      if (length(lib_prefix)>1){
        stop(paste("Error: when running with multilane option, only one library",
                   "can be run at a time"))
      }
      multilane_df <- c()
    }
    for (i in seq_along(R1_filenames)) {
      # this process is done for each Library in each loop.
      if (multilane) {
        lib <- lib_prefix
      } else {
        lib <- lib_prefix[i]
      }
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
        # for each tag combination do cutadapt first for the tag combi
        # and then for cutadapt to remove the primers
        sample_num <- which(ngsfile$V3 == tag_combi)

        # for each loop every tag will be of length one
        fwd_tag <- gsub(":.*", "", tag_combi)
        rev_tag <- gsub(".*:", "", tag_combi)


        # these can be multivalue in case that the sem combi is found multiple times
        filename <- metadata$original_samples[metadata$mjolnir_agnomens %in% ngsfile$V2[sample_num]]
        fwd_outfile <- paste0(filename, # nolint: line_length_linter.
                              R1_motif,".fastq")
        rev_outfile <- gsub(R1_motif, R2_motif, fwd_outfile, fixed = TRUE)
        fwd_outfile <- gsub("..fastq", ".fastq", fwd_outfile, fixed = TRUE)
        rev_outfile <- gsub("..fastq", ".fastq", rev_outfile, fixed = TRUE)
        if(multilane) {
          fwd_outfile_lane <- paste0(fwd_outfile, "_lanefwd_", i)
          rev_outfile_lane <- paste0(rev_outfile, "_lanerev_", i)
          multilane_df <- rbind(multilane_df,
                                data.frame(name = c(fwd_outfile, rev_outfile),
                                           lane = c(i, i),
                                           lane_name = c(fwd_outfile_lane,
                                                         rev_outfile_lane)))
          fwd_outfile <- fwd_outfile_lane
          rev_outfile <- rev_outfile_lane
        }
        fwd_primer <- ngsfile$V4[sample_num]
        rev_primer <- ngsfile$V5[sample_num]
        temp_file_R1 <- paste0('tag1_temp01.fastq')
        temp_file_R2 <- paste0('tag2_temp01.fastq')
        
        # here cutadapt commands for tag and primers
        cutadapt_command_tag <- c()
        cutadapt_command_primers <- c()
        
        cutadapt_command_tag <- paste0(cutadapt_command_tag,
                                        "cutadapt -e 0 ", # allow 0 errors
                                        "-O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required # nolint: line_length_linter.
                                        " --no-indels ", # no indels allowed # nolint: line_length_linter.
                                        "-j ", cores, # number of cores allowed # nolint: line_length_linter.
                                        " --action='none' ", # don't remove the tag taken # nolint: line_length_linter.
                                        # " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                                        " --untrimmed-output ", R1_file_temp, # save those reads that have not been assigned to the sample # nolint: line_length_linter.
                                        " --untrimmed-paired-output ", R2_file_temp, # save those reads that have not been assigned to the sample # nolint: line_length_linter.
                                        " --max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja # nolint: line_length_linter.
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
        if(multilane) {
          R1_file_temp <- paste0("lane_", i, "_", R1_file_temp)
          R2_file_temp <- paste0("lane_", i, "_", R2_file_temp)
        }

        if (fwd_tag != rev_tag) {
          # if this is the case then is a little bit more complex
          # we need to run cutadapt twice changing the order of sample_tags and then concatenate files
          cutadapt_command_tag <- paste0(cutadapt_command_tag, # this is the first command
                                          # the second command change the order of sample_tags
                                          "cutadapt -e 0 ", # allow 0 errors
                                          "-O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required
                                          " --no-indels ", # no indels allowed
                                          "-j ", cores, # number of cores allowed
                                          " --action='none' ", # don't remove the tag taken # nolint: line_length_linter.
                                          # " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                                          " --untrimmed-output ", R1_file_temp, # save those reads that have not been assigned to the sample # nolint: line_length_linter.
                                          " --untrimmed-paired-output ", R2_file_temp, # save those reads that have not been assigned to the sample # nolint: line_length_linter.
                                          " --max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja
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
        if(multilane) {
          R1_file_temp <- paste0("lane_", i, "_", R1_file_temp)
          R2_file_temp <- paste0("lane_", i, "_", R2_file_temp)
        }


        }

        for (j in seq_len(length(sample_num))) {
          cutadapt_command_primers <- paste0(cutadapt_command_primers,
                                            "cutadapt -e 0.1 ", # allow 0.1 errors for primers in ngsfilter was total of 2
                                            # "-O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required # nolint: line_length_linter.
                                            " --no-indels ", # no indels allowed # nolint: line_length_linter.
                                            "-j ", cores, # number of cores allowed # nolint: line_length_linter.
                                            " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                                            " --max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja # nolint: line_length_linter.
                                            "-g ", fwd_primer[j], " -G ", rev_primer[j], # these are the primers # nolint: line_length_linter.
                                            " -o ", fwd_outfile[j], " -p ", rev_outfile[j], # sample names as original # nolint: line_length_linter.
                                            " ", temp_file_R1, " ", temp_file_R2, " ; ") # input files # nolint: line_length_linter.
          if (fwd_tag != rev_tag) {
            cutadapt_command_primers <- paste0(cutadapt_command_primers,
                                              "cutadapt -e 0.1 ", # allow 0.1 errors for primers in ngsfilter was total of 2
                                              # "-O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required # nolint: line_length_linter.
                                              " --no-indels ", # no indels allowed # nolint: line_length_linter.
                                              "-j ", cores, # number of cores allowed # nolint: line_length_linter.
                                              " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                                              " --max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja # nolint: line_length_linter.
                                              "-g ", rev_primer[j], " -G ", fwd_primer[j], # these are the primers # nolint: line_length_linter.
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
    if(multilane) {
      for(final_file in unique(multilane_df$name)) {
        system(paste0("cat ",
                      paste0(multilane_df$lane_name[multilane_df$name == final_file],
                            collapse = " "),
                      ">", final_file, 
                      " ; rm ",
                      paste0(multilane_df$lane_name[multilane_df$name == final_file],
                            collapse = " ")))
    }
  }

  message("Demultiplexing done.")
}
