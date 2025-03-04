#' RAN: Reads Allotment in N samples
#'
#' RAN function will prepare the FASTQ raw data for each sample with the sample
#' tags and primers removed.
#' The R1 output files will contain all forward sequences and
#' R2 the reverse sequences. Please, read the Details section before running.
#'
#' @details
#' RAN considers the following scenarios:
#'
#' - \bold{Scenario 1}. The samples are multiplexed in a library/libraries
#'    (it is possible to run multiple libraries
#'    at once from a single experiment).
#'    Each library has only two raw .fastq(.gz), R1 and R2.
#'    For each library there must be only one ngsfile (see below)
#'    and for the whole experiment only one metadata file (see below).
#'    mjolnir_agnomens (the standard name of the sample during the pipeline)
#'    must be unique, it must not be the same in different ngsfiles.
#'    \emph{lib_prefix} and \emph{R1_filenames} must be the same length.
#'
#' - \bold{Scenario 2}. The samples are multiplexed but for
#'    each library we have more
#'    than two raw data files. For example, when your library has been
#'    sequenced across multiple lanes on a Novaseq. In this case, RAN has
#'    to be run separately for each library and the option
#'    \emph{multilane} must be set to TRUE.
#'    In this case, \emph{lib_prefix} can only be
#'    of length one (only one library),
#'    but the \emph{R1_filenames} can be longer.
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
#'      \bold{ngsfilter_<library identifier>.tsv}.
#'      This must contain five tab-separated
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
#'      file, called "fastq_name_R1".
#'      In the metadata table, each sample in the "original_samples"
#'      column must have a their corresponding fastq_name_R1 and
#'      "mjolnir_agnomen" (LIBX_sample_XXX, i.e LIBA_sample_001).
#'      If the "fastq_name_R1" column is not in the metadata, RAN will use the
#'      \emph{lib_prefix}.
#'
#'
#' @param R1_filenames Character vector with the names of the forward fastq or
#' fastq.gz files. Only needed for multiplexed libraries.
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
#' @param cores Numeric. Number of parts into which the
#' input files will be split
#' for the parallel processing of the FREYJA function.
#'
#' @param R1_motif Character string that distinguishes
#' the forward line file from
#' the reverse.
#'
#' @param R2_motif Character string that distinguishes
#' the reverse line file from
#' the forward.
#'
#' @param metadata_table tsv table. if not specified, the file must be named
#' <EXPX>_metadata.tsv. This table must have: a column named "mjolnir_agnomens"
#' with the names given to the samples during the pipeline in FREYJA; a column
#' named "original_samples" with the samples names that will be given to the
#' samples at the end of the pipeline; and a column with the name specified in
#' in the "blank_col" parameter ("BLANK" by default) where blanks, negatives and
#' controls are tagged with a flag specified in the "blank_tag" parameter (T by
#' default).
#'
#' @param multilane Logical. If FALSE, the function will
#' consider that the each library
#' has only one sequencing lane. If TRUE,
#' only one library is processed at a time.
#'
#' @param tag_error Numeric. From 0 (no errors allowed) to 1. will determine the
#' proportion of nucleotides that can be different in the sample tags.
#'
#' @param primer_error Numeric. From 0 (no errors allowed) to 1.
#' will determine the
#' proportion of nucleotides that can be different in the primers.
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
#' lib_prefix <- c("ULO1", "ULO2", "ULO3", "ULO4")
#'
#' # experiment identifier
#' experiment <- 'ULOY'
#' # Enter number of cores to be used in parallel.
#' cores <- 7
#'
#' mjolnir1_RAN(R1_filenames, lib_prefix = lib_prefix, experiment = experiment,
#'              cores = cores, R1_motif = "_R1", R2_motif = "_R2",
#'              tag_error = 0, primer_error = 0.1)

mjolnir1_RAN <- function(R1_filenames = "",
                         lib_prefix = "",
                         experiment = NULL,
                         primer_F = "GGWACWRGWTGRACWNTNTAYCCYCC",
                         primer_R = "TANACYTCNGGRTGNCCRAARAAYCA",
                         cores = 1,
                         R1_motif = "_R1",
                         R2_motif = "_R2",
                         metadata = "",
                         multilane = FALSE,
                         tag_error = 0,
                         primer_error = 0.1,
                         ...) {
  # Rationale:
  # This function will prepare the FASTQ raw data for each sample.
  # R1 file will contain the fwr seqs and the R2 the rev seqs.
  # Scenarios:
  #    A- samples demultiplexed but primers remain
  #       cutadapt creates two pairs of fastq,
  #       one for each direction and then the files are
  #       concatenated to have all the sequences the same direction.
  #    B- samples multiplexed in a library/libraries
  #       B1 - same sample tag
  #            Here first cutadapt creates a pair of files,
  #            the second creates two pairs that need to be r
  #            edirected for the fwd and rev
  #       B2 - different sample tags
  #            Here the first cutadapt creates a pair of files
  #            that redirecting gives all the correct direction,
  #            the second cutadapt will ensure the primers are
  #            there and trimm them.
  #       B3 - different sample tags with samples sharing the invers
  #            Here the first cutadapt creates two pairs of files.
  #            Concatenating them gives one for the tag A and one for
  #            the tag B. in the second cutadapt those in the tag A
  #            file with fwr primers will belong to one sample and
  #            viceversa. renaming of the R1 and R2 of the second
  #            sample needed.
  args <- list(...)
  
  commands_runned <- c()
  commands_file <- "commands_runned_RAN.txt"
  if (("lib" %in% names(args)) && is.null(experiment)) {
    # Use lib as experiment
    experiment <- args$lib
    # Print deprecation warning
    warning("The 'lib' argument is deprecated.",
            "Please use 'experiment' instead.")
  }
  if (is.null(experiment)){
    stop("Error: experiment argument required")
  }
  if (("lib_prefixes" %in% names(args)) &&
      length(lib_prefix) ==  1 && lib_prefix[1] ==  "") {
    # Use lib as experiment
    lib_prefix <- lib_prefixes
    # Print deprecation warning
    warning(paste("The 'lib_prefixes' argument is deprecated.",
                  "Please use 'lib_prefix' instead.\n",
                  "It's okay by now but don't repeat it, okay?"))
  }
  if (("keep_intermediates" %in% names(args))) {
    keep_intermediates <- args$keep_intermediates
  } else {
    keep_intermediates <- FALSE
  }
  if  (!multilane && length(lib_prefix) !=  length(R1_filenames)){
    stop(paste("Error: you have set the multilane to False but there",
               "are different number of R1_files and lib_prefix"))
  }
  metadata <- read.table(paste0(experiment, "_metadata.tsv"),
                         sep = "\t", header = TRUE)
  # check that metadata has the colnames "original_samples" and "mjolnir_agnomens"
  if (!all(c("original_samples", "mjolnir_agnomens") %in% colnames(metadata))) {
    stop(paste("Error: the metadata file must contain the columns",
               "original_samples and mjolnir_agnomens"))
  }
  if (any(grepl(" ", metadata$original_samples))) {
    stop("ERROR: original_samples cannot have spaces.")
  }
  if(!all(!duplicated(metadata$original_samples))){
    stop("ERROR: original_samples cannot have duplicated names.")
  }
  if (!all(!duplicated(metadata$mjolnir_agnomens))){
    stop("ERROR: mjolnir_agnomens cannot have duplicated names.")
  }
  if ('fastq_name_R1' %in% names(metadata)) {
    # in this case the names for the fastq files are in the metadata which mean
    # that the samples are demultiplexed but the primer needs to be removed
    fastqR1_list <- metadata$fastq_name_R1
    agnomens <-  metadata$original_samples
    
    # check if agnomens have any space
    for (j in seq_len(length(fastqR1_list))) {
      R1_file <- list.files(pattern = paste0("^", fastqR1_list[j]))
      R2_file <- list.files(pattern = paste0("^",gsub(R1_motif, R2_motif, R1_file)))
      # check that files exist
      if (length(R1_file) !=  1) {
        stop(paste("Error: R1 file not found for", R1_filenames[i],
                   "or more than one file matching pattern were found.\n",
                   "found: ", paste(R1_file, collapse = ", ")))
      }
      if (length(R2_file) !=  1) {
        stop(paste("Error: R2 file not found for", R1_filenames[i],
                   "or more than one file matching pattern were found.\n",
                   "found: ", paste(R2_file, collapse = ", ")))
      }
      fwd_outfile <- paste0(agnomens[j],
                            R1_motif,
                            ".fastq")
      rev_outfile <- paste0(agnomens[j],
                            R2_motif,
                            ".fastq")
      fwd_outfile <- gsub("..fastq", ".fastq", fwd_outfile, fixed = TRUE)
      rev_outfile <- gsub("..fastq", ".fastq", rev_outfile, fixed = TRUE)
      if (R1_file ==  fwd_outfile) {
        stop(paste("error: fastq_name_R1 ->", R1_file, " and ",
                   "original_samples ->", fwd_outfile, " are the same. ",
                   "Please consider changing one of them.",
                   "The ouput would be named as the input and in Valhalla",
                   "we don't do that.",
                   "You can use .fastq.gz files as input instead."))
      }
      fwd_fasta <- "fwd_primer.fasta"
      rev_fasta <- "rev_primer.fasta"
      if (primer_F !=  primer_R) {
        fwd_fasta_content <- paste0(">", agnomens[j], "_fwd_rev", "\n",
                                    primer_F, "\n",
                                    ">", agnomens[j], "_rev_fwd", "\n",
                                    primer_R)
        rev_fasta_content <- paste0(">", agnomens[j], "_fwd_rev", "\n",
                                    primer_R, "\n",
                                    ">", agnomens[j], "_rev_fwd", "\n",
                                    primer_F)
        concat_command <- paste0("cat ", agnomens[j], "_fwd_rev_R1.fastq ",
                                 agnomens[j], "_fwd_rev_R2.fastq > ",
                                 fwd_outfile, " ; ",
                                 ifelse(keep_intermediates, "",
                                        paste0(" rm ", agnomens[j], "_fwd_rev_R1.fastq ",
                                               agnomens[j], "_fwd_rev_R2.fastq ; ")),
                                 "cat ", agnomens[j], "_rev_fwd_R1.fastq ",
                                 agnomens[j], "_rev_fwd_R2.fastq > ",
                                 rev_outfile, " ; ",
                                 ifelse(keep_intermediates, "",
                                        paste0(" rm ", agnomens[j], "_rev_fwd_R1.fastq ",
                                               agnomens[j], "_rev_fwd_R2.fastq ; ")))
      } else {
        message(paste("!Warning: both primers are the same.",
                      "This is not a problem but it is not common.",
                      "No redirectioning will be done"))
        fwd_fasta_content <- paste0(">", agnomens[j], "\n", primer_F)
        rev_fasta_content <- paste0(">", agnomens[j], "\n", primer_R)
        if(paste0(agnomens[j], "_R1.fastq") ==  fwd_outfile){
          concat_command <- ""
        } else {
          concat_command <- paste0(" mv ", agnomens[j], "_R1.fastq ",
                                   fwd_outfile, " ; ",
                                   " mv ", agnomens[j], "_R2.fastq ",
                                   rev_outfile, " ; ")
        }
      }
      writeLines(fwd_fasta_content, fwd_fasta)
      writeLines(rev_fasta_content, rev_fasta)
      
      cutadapt_command_primers <-
        paste0(" cutadapt -e ", primer_error, # allow 0.1 errors for primers in ngsfilter was total of 2 # nolint: line_length_linter.
               # " -O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required # nolint: line_length_linter.
               " --no-indels ", # no indels allowed # nolint: line_length_linter.
               " -j ", cores, # number of cores allowed # nolint: line_length_linter.
               " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
               " --max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja # nolint: line_length_linter.
               " --pair-adapters ",
               " -g file:", fwd_fasta, " -G file:", rev_fasta, # these are the sample_tags # for demultiplexed not anchored # nolint: line_length_linter.
               " -o {name}_R1.fastq -p {name}_R2.fastq", # sample names as original # nolint: line_length_linter.
               " ", R1_file, " ", R2_file, " ; ",
               concat_command) # input files # nolint: line_length_linter.
      commands_runned <- c(commands_runned, cutadapt_command_primers)
      writeLines(commands_runned, commands_file)
      sapply(cutadapt_command_primers,
             FUN = system, wait = TRUE, intern = TRUE)
    }
  } else {
    # in this case the samples are not demultiplexed.
    # first the samples are demultiplexed with the sample tags
    # then the primers are removed from the demultiplexed samples.
    
    message("RAN will demultiplex the following initial FASTQ files:")
    message(paste(R1_filenames,collapse = ", "))
    ngsfilter_files <- list.files('.', 'ngsfilter')
    if (length(ngsfilter_files) ==  0) {
      stop("No ngsfilter file found.")
    }
    if (multilane) {
      if (length(lib_prefix)>1){
        stop(paste("Error: when running with multilane option,",
                   "only one library can be run at a time"))
      }
      # multilane_df <- c()
    }
    for (i in seq_along(R1_filenames)) {
      # this process is done for each Library in each loop.
      if (!('original_samples' %in% names(metadata)) ||
          !('mjolnir_agnomens' %in% names(metadata))) {
        stop(paste("Error: the metadata file must contain the column",
                   "original_samples and mjolnir_agnomens"))
      }
      if (multilane) {
        lib <- lib_prefix
        multilane_string <- paste0("_multilane_", i)
        
      } else {
        lib <- lib_prefix[i]
        multilane_string <- ""
      }
      
      R1_file <- list.files(pattern = R1_filenames[i])
      R1_file <- unique(gsub("_untrimmed","",R1_file))
      R2_file <- list.files(pattern = gsub(R1_motif, R2_motif, R1_file))
      R2_file <- unique(gsub("_untrimmed","",R2_file))
      # check that files exist
      if (length(R1_file) !=  1) {
        stop(paste("Error: R1 file not found for", R1_filenames[i],
                   "or more than one file matching pattern were found.\n",
                   "found: ", paste(R1_file, collapse = ", ")))
      }
      if (length(R2_file) !=  1) {
        stop(paste("Error: R2 file not found for", R1_filenames[i],
                   "or more than one file matching pattern were found.\n",
                   "found: ", paste(R2_file, collapse = ", ")))
      }
      message(paste0("RAN is processing the ", lib, " library."))
      ngsfile <- read.csv(ngsfilter_files[grep(lib, ngsfilter_files)],
                          sep = "\t", header = FALSE)
      message(paste0("RAN will split initial ", R1_file, " & ",
                     R2_file,
                     " files in ", dim(ngsfile)[1], " samples."))
      if (sum(duplicated(paste0(ngsfile$V3,ngsfile$V4,ngsfile$V5))) > 0) {
        stop(paste("There are duplicated sample tags",
                   "with the same primers in the ",
                   ngsfilter_files[grep(lib, ngsfilter_files)], " file."))
      }
      if (sum(duplicated(ngsfile$V2)) > 0) {
        stop(paste("There are duplicated mjolnir_agnomens",
                   "in the ", ngsfilter_files[grep(lib, ngsfilter_files)],
                   " file."))
      }
      # check for those tags that are inverted in different samples
      ngsfile$fwd_tags <- gsub(":.*", "", ngsfile$V3)
      ngsfile$rev_tags <- gsub(".*:", "", ngsfile$V3)
      ngsfile$same_tag <- ngsfile$fwd_tags ==  ngsfile$rev_tags
      ngsfile$in_two_steps <- FALSE
      ngsfile$in_two_steps[!ngsfile$same_tag] <-
        paste0(ngsfile$rev_tags[!ngsfile$same_tag],
               ":", ngsfile$fwd_tags[!ngsfile$same_tag])  %in% ngsfile$V3
      
      fwd_fasta <- "fwd_tag.fasta"
      rev_fasta <- "rev_tag.fasta"
      fwd_tag_fasta_content <- ""
      rev_tag_fasta_content <- ""
      ngsfile_same_tag <- ngsfile[ngsfile$same_tag,]
      ngsfile_not_same_tag <-
        ngsfile[!ngsfile$same_tag & !ngsfile$in_two_steps,]
      ngsfile_two_steps <- ngsfile[!ngsfile$same_tag & ngsfile$in_two_steps,]
      if(nrow(ngsfile_same_tag) > 0){
        fwd_tag_fasta_content <- paste0(">", ngsfile_same_tag$V2, "_st", "\n", # st for same tag
                                        ngsfile_same_tag$fwd_tags)
        rev_tag_fasta_content <- paste0(">", ngsfile_same_tag$V2, "_st", "\n",
                                        ngsfile_same_tag$rev_tags)
      }
      if(nrow(ngsfile_not_same_tag) > 0){
        fwd_tag_fasta_content <-
          c(fwd_tag_fasta_content,
            paste0(">", ngsfile_not_same_tag$V2, "_fwd_rev", "\n",
                   ngsfile_not_same_tag$fwd_tags),
            paste0(">", ngsfile_not_same_tag$V2, "_rev_fwd", "\n",
                   ngsfile_not_same_tag$rev_tags))
        rev_tag_fasta_content <-
          c(rev_tag_fasta_content,
            paste0(">", ngsfile_not_same_tag$V2, "_fwd_rev", "\n",
                   ngsfile_not_same_tag$rev_tags),
            paste0(">", ngsfile_not_same_tag$V2, "_rev_fwd", "\n",
                   ngsfile_not_same_tag$fwd_tags))
      }
      if(nrow(ngsfile_two_steps) > 0){
        ngsfile_two_steps_bis <- ngsfile_two_steps
        for(j in seq_len(nrow(ngsfile_two_steps_bis)/2)){
          combi_tags_A <-
            paste0(ngsfile_two_steps_bis$fwd_tags[1], ":",
                   ngsfile_two_steps_bis$rev_tags[1])
          combi_tags_B <-
            paste0(ngsfile_two_steps_bis$rev_tags[1], ":",
                   ngsfile_two_steps_bis$fwd_tags[1])
          ngsfile_two_steps_bis$processing <- FALSE
          ngsfile_two_steps_bis$processing[
            ngsfile_two_steps_bis$V3 ==  combi_tags_A] <- TRUE
          ngsfile_two_steps_bis$processing[
            ngsfile_two_steps_bis$V3 ==  combi_tags_B] <- TRUE
          
          fwd_tag_fasta_content <-
            c(fwd_tag_fasta_content,
              paste0(">",
                     paste(
                       ngsfile_two_steps_bis$V2[
                         ngsfile_two_steps_bis$processing],
                       collapse = "_"), "_combinedA", "\n",
                     ngsfile_two_steps_bis$fwd_tags[1]),
              paste0(">",
                     paste(
                       ngsfile_two_steps_bis$V2[
                         ngsfile_two_steps_bis$processing],
                       collapse = "_"), "_combinedB", "\n",
                     ngsfile_two_steps_bis$rev_tags[1]))
          rev_tag_fasta_content <-
            c(rev_tag_fasta_content,
              paste0(">",
                     paste(
                       ngsfile_two_steps_bis$V2[
                         ngsfile_two_steps_bis$processing],
                       collapse = "_"), "_combinedA", "\n",
                     ngsfile_two_steps_bis$rev_tags[1]),
              paste0(">",
                     paste(
                       ngsfile_two_steps_bis$V2[
                         ngsfile_two_steps_bis$processing],
                       collapse = "_"), "_combinedB", "\n",
                     ngsfile_two_steps_bis$fwd_tags[1]))
          ngsfile_two_steps_bis <-
            ngsfile_two_steps_bis[!ngsfile_two_steps_bis$processing,]
        }
      }
      
      writeLines(fwd_tag_fasta_content, fwd_fasta)
      writeLines(rev_tag_fasta_content, rev_fasta)
      
      cutadapt_command_tag <-
        paste0(" cutadapt -e ", tag_error, # allow 0 errors
               " -O ", min(c(nchar(ngsfile$fwd_tags), nchar(ngsfile$rev_tags))), # min overlap required # nolint: line_length_linter.
               " --no-indels ", # no indels allowed # nolint: line_length_linter.
               " -j ", cores, # number of cores allowed # nolint: line_length_linter.
               " --action = 'none' ", # don't remove the tag taken # nolint: line_length_linter.
               # " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
               " --untrimmed-output ", R1_file, "_untrimmed", # save those reads that have not been assigned to the sample # nolint: line_length_linter.
               " --untrimmed-paired-output ", R2_file, "_untrimmed", # save those reads that have not been assigned to the sample # nolint: line_length_linter.
               " --max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja # nolint: line_length_linter.
               " --pair-adapters ",
               " -g file:", fwd_fasta,
               " -G file:", rev_fasta, # these are the sample_tags # nolint: line_length_linter.
               " -o {name}_R1.fastq ",
               " -p {name}_R2.fastq",
               # " -o ", temp_file_R1, " -p ", temp_file_R2, # sample names as original # nolint: line_length_linter.
               " ", R1_file, " ", R2_file, " ; ") # input files # nolint: line_length_linter.
      
      commands_runned <- c(commands_runned, cutadapt_command_tag)
      writeLines(commands_runned, commands_file)
      sapply(cutadapt_command_tag,
             FUN = system, wait = TRUE, intern = TRUE)
      # now the primers are removed
      if (nrow(ngsfile_same_tag) > 0) {
        for (tag_combi in seq_along(ngsfile_same_tag$V3)) {
          # for each tag combination do cutadapt first for the tag combi
          # and then for cutadapt to remove the primers
          sample_name <-
            metadata$original_samples[metadata$mjolnir_agnomens %in%
                                        ngsfile_same_tag$V2[tag_combi]]
          mjolnir_name <- ngsfile_same_tag$V2[tag_combi]
          # check if the file has reads by checking the size of the file
          if (file.size(paste0(mjolnir_name, "_st_R1.fastq")) ==  0) {
            message(paste("Warning: the sample", sample_name, "has no reads."))
            next()
          }
          fwd_outfile <- paste0(sample_name, # nolint: line_length_linter.
                                R1_motif, multilane_string, ".fastq")
          rev_outfile <- paste0(sample_name, # nolint: line_length_linter.
                                R2_motif, multilane_string, ".fastq")
          # for each loop every tag will be of length one
          
          fwd_primer_fasta <- "fwd_primer.fasta"
          rev_primer_fasta <- "rev_primer.fasta"
          fwd_primer_fasta_content <-
            c(paste0(">", sample_name, "_fr\n",# this is the sequences sorted fwd-rev
                     ngsfile_same_tag$V4[tag_combi]),
              paste0(">", sample_name, "_rf\n",
                     ngsfile_two_steps$V5[tag_combi]))
          rev_primer_fasta_content <-
            c(paste0(">", sample_name, "_fr\n",# this is the rev in the R2
                     ngsfile_same_tag$V5[tag_combi]),
              paste0(">", sample_name, "_rf\n",
                     ngsfile_two_steps$V4[tag_combi]))
          
          writeLines(fwd_primer_fasta_content, fwd_primer_fasta)
          writeLines(rev_primer_fasta_content, rev_primer_fasta)
          cutadapt_command_primers <-
            paste0(" cutadapt -e ", primer_error, # allow 0.1 errors for primers in ngsfilter was total of 2 # nolint: line_length_linter.
                   # " -O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required # nolint: line_length_linter.
                   " --no-indels ", # no indels allowed # nolint: line_length_linter.
                   " -j ", cores, # number of cores allowed # nolint: line_length_linter.
                   # " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                   " --untrimmed-output ",
                   fwd_outfile, "_untrimmed", # save those reads that have not been assigned to the sample # nolint: line_length_linter.
                   " --untrimmed-paired-output ",
                   rev_outfile, "_untrimmed", # save those reads that have not been assigned to the sample # nolint: line_length_linter.
                   " --max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja # nolint: line_length_linter.
                   " --pair-adapters ",
                   " -g file:", fwd_primer_fasta,
                   " -G file:", rev_primer_fasta, # these are the primers # nolint: line_length_linter.
                   " -o {name}_R1.fastq ",
                   " -p {name}_R2.fastq", # sample names as original # nolint: line_length_linter.
                   " ", mjolnir_name, "_st_R1.fastq ", mjolnir_name, "_st_R2.fastq ; ", # input files # nolint: line_length_linter.
                   # concatenate the forward files
                   " cat ", sample_name, "_fr_R1.fastq ", # nolint: line_length_linter.
                   sample_name, "_rf_R2.fastq >",
                   fwd_outfile, " ; ", # nolint: line_length_linter.
                   # concatenate the reverse files
                   " cat ", sample_name, "_fr_R2.fastq ", # nolint: line_length_linter.
                   sample_name, "_rf_R1.fastq >",
                   rev_outfile, " ; ", # nolint: line_length_linter.
                   # remove the temp files
                   " rm " ,
                   paste0(sample_name,
                          c("_fr_R1.fastq ", "_rf_R2.fastq", # nolint: line_length_linter.
                            "_fr_R2.fastq ", "_rf_R1.fastq"), # nolint: line_length_linter.
                          collapse = " "), # nolint: line_length_linter.
                   " ; ") # nolint: line_length_linter.
          if (!keep_intermediates) {
            cutadapt_command_primers <-
              paste0(cutadapt_command_primers,
                     " rm ", mjolnir_name, "_st_R1.fastq ",
                     mjolnir_name, "_st_R2.fastq ; ")
          }
          commands_runned <- c(commands_runned, cutadapt_command_primers)
          writeLines(commands_runned, commands_file)
          sapply(cutadapt_command_primers,
                 FUN = system, wait = TRUE, intern = TRUE)
          if(multilane){
            concat_multilane <-
              paste0("cat ",
                     fwd_outfile, " >> ", sample_name, R1_motif, ".fastq ; ",
                     "cat ",
                     rev_outfile, " >> ", sample_name, R2_motif, ".fastq ; ",
                     "rm ", fwd_outfile, " ", rev_outfile, " ; ")
            commands_runned <- c(commands_runned, concat_multilane)
            writeLines(commands_runned, commands_file)
            sapply(concat_multilane,
                   FUN = system, wait = TRUE, intern = TRUE)
          }
        }
      }
      if(nrow(ngsfile_not_same_tag) > 0){
        for (tag_combi in seq_along(ngsfile_not_same_tag$V3)) {
          # for each tag combination do cutadapt first for the tag combi
          # and then for cutadapt to remove the primers
          sample_name <-
            metadata$original_samples[metadata$mjolnir_agnomens %in%
                                        ngsfile_not_same_tag$V2[tag_combi]]
          mjolnir_name <- ngsfile_not_same_tag$V2[tag_combi]
          # check if the file has reads by checking the size of the file
          if (file.size(paste0(mjolnir_name, "_fwd_rev_R1.fastq")) ==  0 &&
              file.size(paste0(mjolnir_name, "_rev_fwd_R1.fastq")) ==  0) {
            message(paste("Warning: the sample", sample_name, "has no reads."))
            next()
          }
          fwd_outfile <- paste0(sample_name, # nolint: line_length_linter.
                                R1_motif, multilane_string, ".fastq")
          rev_outfile <- paste0(sample_name, # nolint: line_length_linter.
                                R2_motif, multilane_string, ".fastq")
          # for each loop every tag will be of length one
          
          fwd_primer_fasta <- "fwd_primer.fasta"
          rev_primer_fasta <- "rev_primer.fasta"
          fwd_primer_fasta_content <-
            c(paste0(">", sample_name, "\n",# this is the sequences sorted fwd-rev
                     ngsfile_not_same_tag$V4[tag_combi]))
          rev_primer_fasta_content <-
            c(paste0(">", sample_name, "\n",# this is the rev in the R2
                     ngsfile_not_same_tag$V5[tag_combi]))
          
          writeLines(fwd_primer_fasta_content, fwd_primer_fasta)
          writeLines(rev_primer_fasta_content, rev_primer_fasta)
          cutadapt_command_primers <-
            paste0(" cat ", mjolnir_name, "_fwd_rev_R1.fastq ",
                   mjolnir_name, "_rev_fwd_R2.fastq > ",
                   mjolnir_name, "_fwd.fastq ; ",
                   " cat ", mjolnir_name, "_fwd_rev_R2.fastq ",
                   mjolnir_name, "_rev_fwd_R1.fastq > ",
                   mjolnir_name, "_rev.fastq ; ",
                   " rm ", paste0(mjolnir_name,
                                  c("_fwd_rev_R1.fastq ", "_rev_fwd_R2.fastq", # nolint: line_length_linter.
                                    "_fwd_rev_R2.fastq ", "_rev_fwd_R1.fastq"), # nolint: line_length_linter.
                                  collapse = " "), " ; ", # nolint: line_length_linter.
                   " cutadapt -e ", primer_error, # allow 0.1 errors for primers in ngsfilter was total of 2
                   # " -O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required # nolint: line_length_linter.
                   " --no-indels ", # no indels allowed # nolint: line_length_linter.
                   " -j ", cores, # number of cores allowed # nolint: line_length_linter.
                   # " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                   " --untrimmed-output ",
                   fwd_outfile, "_untrimmed", # save those reads that have not been assigned to the sample # nolint: line_length_linter.
                   " --untrimmed-paired-output ",
                   rev_outfile, "_untrimmed", # save those reads that have not been assigned to the sample # nolint: line_length_linter.
                   " --max-n=0.5 ", # I allow a max of half of the read being N. this will be solved by freyja # nolint: line_length_linter.
                   " --pair-adapters ",
                   " -g file:", fwd_primer_fasta,
                   " -G file:", rev_primer_fasta, # these are the primers # nolint: line_length_linter.
                   " -o {name}", R1_motif, multilane_string, ".fastq ",
                   " -p {name}", R2_motif, multilane_string, ".fastq", # sample names as original # nolint: line_length_linter.
                   " ", mjolnir_name, "_fwd.fastq ", mjolnir_name, "_rev.fastq ; " # input files # nolint: line_length_linter.
            ) # nolint: line_length_linter.
          if (!keep_intermediates) {
            cutadapt_command_primers <-
              paste0(cutadapt_command_primers,
                     " rm ", mjolnir_name, "_fwd.fastq ", mjolnir_name, "_rev.fastq ; ")
          }
          commands_runned <- c(commands_runned, cutadapt_command_primers)
          writeLines(commands_runned, commands_file)
          sapply(cutadapt_command_primers,
                 FUN = system, wait = TRUE, intern = TRUE)
          if(multilane){
            concat_multilane <-
              paste0("cat ",
                     fwd_outfile, " >> ", sample_name, R1_motif, ".fastq ; ",
                     "cat ",
                     rev_outfile, " >> ", sample_name, R2_motif, ".fastq ; ",
                     "rm ", fwd_outfile, " ", rev_outfile, " ; ")
            commands_runned <- c(commands_runned, concat_multilane)
            writeLines(commands_runned, commands_file)
            sapply(concat_multilane,
                   FUN = system, wait = TRUE, intern = TRUE)
          }
        }
      }
      if(nrow(ngsfile_two_steps) > 0){
        ngsfile_two_steps_bis <- ngsfile_two_steps
        for(j in seq_len(nrow(ngsfile_two_steps_bis)/2)){
          combi_tags_A <-
            paste0(ngsfile_two_steps_bis$fwd_tags[1], ":",
                   ngsfile_two_steps_bis$rev_tags[1])
          combi_tags_B <-
            paste0(ngsfile_two_steps_bis$rev_tags[1], ":",
                   ngsfile_two_steps_bis$fwd_tags[1])
          ngsfile_two_steps_bis$processing <- FALSE
          ngsfile_two_steps_bis$processing[ngsfile_two_steps_bis$V3 == 
                                             combi_tags_A] <- TRUE
          ngsfile_two_steps_bis$processing[ngsfile_two_steps_bis$V3 == 
                                             combi_tags_B] <- TRUE
          
          prev_file_combiA_R1 <-
            paste0(
              paste(ngsfile_two_steps_bis$V2[ngsfile_two_steps_bis$processing],
                    collapse = "_"),
              "_combinedA_R1.fastq")
          prev_file_combiA_R2 <-
            paste0(
              paste(ngsfile_two_steps_bis$V2[ngsfile_two_steps_bis$processing],
                    collapse = "_"),
              "_combinedA_R2.fastq")
          prev_file_combiB_R1 <-
            paste0(
              paste(ngsfile_two_steps_bis$V2[ngsfile_two_steps_bis$processing],
                    collapse = "_"),
              "_combinedB_R1.fastq")
          prev_file_combiB_R2 <-
            paste0(
              paste(ngsfile_two_steps_bis$V2[ngsfile_two_steps_bis$processing],
                    collapse = "_"),
              "_combinedB_R2.fastq")
          mjolnir_name_1st <-
            ngsfile_two_steps_bis$V2[ngsfile_two_steps_bis$processing][1]
          mjolnir_name_2nd <-
            ngsfile_two_steps_bis$V2[ngsfile_two_steps_bis$processing][2]
          sample_name_1st <-
            metadata$original_samples[metadata$mjolnir_agnomens %in%
                                        mjolnir_name_1st]
          sample_name_2nd <-
            metadata$original_samples[metadata$mjolnir_agnomens %in%
                                        mjolnir_name_2nd]
          # check if the file has reads by checking the size of the file
          if (file.size(prev_file_combiA_R1) ==  0 &&
              file.size(prev_file_combiB_R1) ==  0) {
            message(paste("Warning: the samples", sample_name_1st, "&",
                          sample_name_2nd, "have no reads."))
            next()
          }
          cutadapt_command_primers <-
            paste0(" cat ", prev_file_combiA_R1, " ",
                   prev_file_combiB_R2, " > ",
                   "combined_tagA.fastq ; ",
                   " cat ", prev_file_combiA_R2, " ",
                   prev_file_combiB_R1, " > ",
                   "combined_tagB.fastq ; ",
                   ifelse(keep_intermediates, "",
                          paste0(" rm ",
                                 prev_file_combiA_R1, " ", # nolint: line_length_linter.
                                 prev_file_combiB_R2, " ", # nolint: line_length_linter.
                                 prev_file_combiA_R2, " ", # nolint: line_length_linter.
                                 prev_file_combiB_R1, " ; ")) # nolint: line_length_linter.
            ) # nolint: line_length_linter.
          fwd_primer_1st <-
            ngsfile_two_steps_bis$V4[ngsfile_two_steps_bis$processing][1]
          fwd_primer_2nd <-
            ngsfile_two_steps_bis$V4[ngsfile_two_steps_bis$processing][2]
          rev_primer_1st <-
            ngsfile_two_steps_bis$V5[ngsfile_two_steps_bis$processing][1]
          rev_primer_2nd <-
            ngsfile_two_steps_bis$V5[ngsfile_two_steps_bis$processing][2]
          fwd_primer_fasta <- "fwd_primer.fasta"
          rev_primer_fasta <- "rev_primer.fasta"
          fwd_primer_fasta_content <-
            c(paste0(">", sample_name_1st, "\n",# this is the sequences sorted fwd-rev
                     fwd_primer_1st),
              paste0(">", sample_name_2nd, "_rev\n",
                     rev_primer_2nd))
          rev_primer_fasta_content <-
            c(paste0(">", sample_name_1st, "\n",# this is the rev in the R2
                     rev_primer_1st),
              paste0(">", sample_name_2nd, "_rev\n",
                     fwd_primer_2nd))
          writeLines(fwd_primer_fasta_content, fwd_primer_fasta)
          writeLines(rev_primer_fasta_content, rev_primer_fasta)
          cutadapt_command_primers <-
            c(cutadapt_command_primers,
              paste0(" cutadapt -e ", primer_error, # allow 0.1 errors for primers in ngsfilter was total of 2
                     # " -O ", min(c(nchar(fwd_tag), nchar(rev_tag))), # min overlap required # nolint: line_length_linter.
                     " --no-indels ", # no indels allowed # nolint: line_length_linter.
                     " -j ", cores, # number of cores allowed # nolint: line_length_linter.
                     # " --discard-untrimmed ", # discard those reads that have not been assigned to the sample # nolint: line_length_linter.
                     " --untrimmed-output ",
                     sample_name_1st, "_combinedA_R1.fastq_untrimmed", # save those reads that have not been assigned to the sample # nolint: line_length_linter.
                     " --untrimmed-paired-output ",
                     sample_name_1st, "_combinedA_R2.fastq_untrimmed", # save those reads that have not been assigned to the sample # nolint: line_length_linter.
                     " --max-n=0.5 ", # I allow a max of half of the read being N
                     " --pair-adapters ",
                     " -g file:", fwd_primer_fasta,
                     " -G file:", rev_primer_fasta, # these are the primers # nolint: line_length_linter.
                     " -o {name}", R1_motif, multilane_string, ".fastq ",
                     " -p {name}", R2_motif, multilane_string, ".fastq", # sample names as original # nolint: line_length_linter.
                     " combined_tagA.fastq combined_tagB.fastq ; ", # input files # nolint: line_length_linter.
                     " mv ",
                     sample_name_2nd, "_rev", R1_motif, multilane_string, ".fastq ./",
                     sample_name_2nd, R2_motif, multilane_string, ".fastq ; ",
                     " mv ",
                     sample_name_2nd, "_rev", R2_motif, multilane_string, ".fastq ./",
                     sample_name_2nd, R1_motif, multilane_string, ".fastq ; ") # nolint: line_length_linter.
            )
          if (!keep_intermediates) {
            cutadapt_command_primers <-
              paste0(cutadapt_command_primers,
                     " rm combined_tagA.fastq combined_tagB.fastq ; ")
          }
          commands_runned <- c(commands_runned, cutadapt_command_primers)
          writeLines(commands_runned, commands_file)
          sapply(cutadapt_command_primers,
                 FUN = system, wait = TRUE, intern = TRUE)
          if(multilane){
            concat_multilane <-
              paste0("cat ",
                     sample_name_1st, R1_motif, multilane_string, ".fastq >> ",
                     sample_name_1st, R1_motif, ".fastq ; ",
                     "cat ",
                     sample_name_1st, R2_motif, multilane_string, ".fastq >> ",
                     sample_name_1st, R2_motif, ".fastq ; ",
                     "rm ",
                     sample_name_1st, R1_motif, multilane_string, ".fastq ",
                     sample_name_1st, R2_motif, multilane_string, ".fastq ; ",
                     "cat ",
                     sample_name_2nd, R1_motif, multilane_string, ".fastq >> ",
                     sample_name_2nd, R1_motif, ".fastq ; ",
                     "cat ",
                     sample_name_2nd, R2_motif, multilane_string, ".fastq >> ",
                     sample_name_2nd, R2_motif, ".fastq ; ",
                     "rm ",
                     sample_name_2nd, R1_motif, multilane_string, ".fastq ",
                     sample_name_2nd, R2_motif, multilane_string, ".fastq ; ")
            commands_runned <- c(commands_runned, concat_multilane)
            writeLines(commands_runned, commands_file)
            sapply(concat_multilane,
                   FUN = system, wait = TRUE, intern = TRUE)
          }
          ngsfile_two_steps_bis <-
            ngsfile_two_steps_bis[!ngsfile_two_steps_bis$processing,]
        }
      }
      rm(ngsfile, ngsfile_same_tag, ngsfile_not_same_tag, ngsfile_two_steps)
    }
    message("Demultiplexing done.")
    
  }
  message("RAN finished.")
}
