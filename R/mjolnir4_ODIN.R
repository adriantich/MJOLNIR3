#' ODIN: OTU Delimitation Inferred by Networks
#'
#' ODIN performs MOTU clustering and/or denoising. It is one of the main steps of MJOLNIR.
#'
#' @details
#' The function mjolnir4_ODIN() uses the four different strategies to delimit
#' MOTUs and/or ESVs. This strategies are set with the algorithm parameter:
#' 
#' a)"DnoisE_SWARM", b)"SWARM", c)"SWARM_DnoisE" and d)"DnoisE".
#'
#' In short, DnoisE refers to the denoising process with DnoisE to obtain ESV
#' and SWARM to a clustering process with SWARM to obtain MOTUs. DnoisE is a
#' software to merge spurious sequences into their "mothers" (see Antich et al.
#' 2022) to obtain Exact (also Amplicon) Sequence variants. DnoisE is an open
#' source and parallelizable alternative to Unoise that allows to introduce an
#' entropy correction based on the different entropies of each position in the
#' codon of coding genes. This is highly recommended for markers as COI for
#' which this program was intended. However, with the entropy=FALSE parameter,
#' this programs performs the same denoising procedure as described for Unoise3.
#' SWARM is an algorithm to delimit MOTUs, based on linkage-networks created by
#' step-by-step agregation. This clustering algorithm is not based on an
#' arbitrary, constant, absolute identity threshold. Conversely, SWARM is based
#' on an iterative aggregation of sequences that differ less than a given
#' distance d. This strategy results into linkage-networks of different sizes,
#' which can have different effective values for within-MOTU identity threshold,
#' depending on the complexity of the natural variability of the sequences
#' present in the sample. This procedure is very convenient in the case of
#' hypervariable metabarcoding markers, such as COI, which usually feature
#' extremely high levels of natural diversity, in addition to the random sequencing errors.
#' Dereplication step takes place after joining all samples into the same file
#' before SWARM in algorithms a, b and c and after DnoisE in algorithms a and d
#' 
#' NEW: now ODIN performs most of the filters that were applied in RAGNAROC. The
#' idea is to retrieve an output file for each sample that is independent from the
#' rest of the samples (blank and neg corrected; total reads filtered) and with
#' denoised ASV/ESV using the name 
#' 
#' "<mjolnir_agnoment>_ODIN_ESV_<original_name>.fasta.".
#' 
#' To do so these are the the different steps of this function:
#' 
#' <<D -> "DnoisE";
#' 
#' DS -> "DnoisE_SWARM"
#' 
#' S -> "SWARM";
#' 
#' SD -> "SWARM_DnoisE";
#' 
#' SaD (SWARM after DnoisE) --> "DnoisE_SWARM" & run_dnoise = F>>
#' 
#' 0: define variables and load metadata table
#' 
#' 1: D and DS -> denoise the fasta files
#' 
#' 2: D,DS,SD,S,SaD -> cat all fasta
#' 
#' 3: D,DS,SD,S,SaD -> dereplicate
#' 
#' 4: D,DS,SD,S,SaD -> annotate new names
#'
#'
#' filter1A: D,DS -> blank relative abundances filter AS MOTU
#' 
#' filter2A: D,DS -> min relative abundances filter WS ESV
#' 
#' filter3A: D,DS -> remove nletons # improves SWARM for DS
#' 
#' filter1B: SD,S,SaD -> remove nletons # this will improve SWARM performance
#'                                      # apply this to SaD just in case
#'
#' 5: D,DS,SD,S,SaD -> export csv file of sequences (if denoised, they are ESV)
#' 
#' 6: D,DS -> export fasta DnoisEd & blank filt
#' 
#' 7: D -> create fasta files for taxonomic assignment
#' 
#' 8: DS,SD,S,SaD -> do SWARM
#'
#'
#' filter2B: SD,S -> blank relative abundances filter AS MOTU
#' 
#' filter3B: SD,S -> min relative abundances filter WS MOTU
#' 
#' filter4B: SD,S -> remove nletons AS MOTUs # remove artefacts
#'
#'
#' 9: DS,SD,S,SaD -> create csv files of MOTUs and ESV or seqs.
#' 
#' 10: DS,SD,S,SaD -> Also create fasta files for taxonomic assignment
#' 
#' 11: SD -> run DnoisE over the csv files
#' 
#' Blank filter: remove any MOTU for which abundance in the blank or negative 
#' controls is higher than "blank_relative" of its total read abundance
#' and remove blank and NEG samples
#' 
#' Minimum relative abundance filter: Apply a minimum relative abundance 
#' threshold for each sample, setting to zero any abundance below "min_relative" 
#' of the total reads of this sample. It also applies a "min_reads" filter
#'
#' @param experiment Character string. Acronym for the experiment. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. However they can be the same.
#'
#' @param cores Numeric. Number of threads for parallel processing.
#'
#' @param d Numeric value for d parameter of SWARM that refers to the maximum
#' number of differences between two sequences to be linked in the same cluster
#'
#' @param min_reads_MOTU Numeric. Minimum number of reads that a MOTU needs to
#' have to be retained.
#'
#' @param min_reads_ESV Numeric. Minimum number of reads that an ESV needs to
#' have to be retained. Also works works for non-denoised sequences. In the
#' case of algorithms involving SWARM, the removal is performed before SWARM.
#'
#' @param min_relative Number of the minimum relative abundance for a unit in the
#' sample to be retained.
#' 
#' @param blank_relative Relative abundance threshold for a unit to be removed if the 
#' total abundance in the blank/neg/control is higher than this value in terms of
#' relative abundance of the total reads in all samples (see Details).
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
#' @param blank_col Column name of the blank column in the "metadata_table"
#' 
#' @param blank_tag Unique flag to tag the samples that are blank/neg/controls 
#' in the "blank_col"
#' 
#' @param alpha Numeric. Alpha value for DnoisE to run.
#'
#' @param entropy Logical, numeric or character vector specifying whether to run
#' DnoisE with entropy correction and how.
#'
#' a) c(0.47,0.23,1.02,313) - formulation refers to the entropy values for the
#' first, second and third position in the codon and the length of the main
#' sequence length expected. See Antich et al. 2022 for further details.
#'
#' b) FALSE - this will disable the entropy correction. Recommended for non
#' coding markers
#'
#' c) c("auto_sample",313) - this will compute the entropy values for 313
#' (plus or minus multiple of 3) bp within DnoisE and use them to perform the entropy correction
#'
#' d) c("auto_dataset") - this will compute the entropy values for all the
#' dataset and all sequence lengths and use the main sequence length's values for
#' the entropy correction
#'
#' @param algorithm Character. It specifies the algorithm to obtain MOTUs and/or ESVs. Ther are four options:
#'
#' a)"DnoisE_SWARM" - This option will run DnoisE before SWARM. This option is the best
#' choice for highly diverse data sets so it will reduce the computation time of
#' SWARM. It also allows the analysis of metaphylogeographical approaches and
#' retrieves denoised and quality filter fasta files for each sample that can
#' be used for other experiment without previous run of RAN, FREYJA and HELA. It
#' will result on a table of MOTUs and the ESV clustered into them.
#'
#' b)"SWARM" - This option will run only SWARM to obtain a MOTU table.
#'
#' c)"SWARM_DnoisE" - This option will run SWARM before DnoisE. This option
#' responds to a philosophical point of view where the denoising of the
#' sequences have to be performed within MOTUs so closer sequences are compared
#' and to avoid that high abundant sequences from different MOTUs absorb
#' sequences (and thus their reads) from different MOTUs. However, there are no major
#' differences between this option and the "DnoisE_SWARM" option. It allows the
#' analysis of metaphylogeographical and will result on a table of MOTUs and
#' the ESV clustered into them.
#'
#' d)"DnoisE" - This option will run only DnoisE. This option will retrieve
#' quality filter fasta files for each sample that can
#' be used for other experiment without previous run of RAN, FREYJA and HELA.
#'
#' @param run_dnoise Logical. In the case of the algorithm='DnoisE_SWARM' there is
#' the option of not running the DnoisE if it has been already run for a previous
#' experiment. The denoised and filtered fasta files are needed.
#'
#' @param remove_DMS Logical. If TRUE, it will delete all obidms objects that are
#' created during the process. This can save a lot of hard disk space. The FALSE
#' option is useful for developing and debugging.
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


mjolnir4_ODIN <- function(experiment = NULL, cores = 1, d = 13,
                          min_reads_MOTU = 2, min_reads_ESV = 2,
                          min_relative = 1 / 50000,
                          blank_relative = 0.1,
                          metadata_table = "",
                          blank_col = "BLANK", blank_tag = TRUE,
                          alpha = 4,
                          entropy = c(0.47, 0.23, 1.02, 313),
                          algorithm = "DnoisE_SWARM", run_dnoise = TRUE,
                          remove_singletons = NULL, remove_DMS = TRUE, ...) {
  
  
  #####
  # the steps in this function are:
  # 0: define variables and load metadata table
  # 1: D and DS -> denoise the fasta files
  # 2: D,DS,SD,S,SaD -> cat all fasta
  # 3: D,DS,SD,S,SaD -> dereplicate
  # 4: D,DS,SD,S,SaD -> annotate new names
  #
  # filter1A: D,DS -> blank relative abundances filter AS MOTU
  # filter2A: D,DS -> min relative abundances filter WS ESV
  # filter3A: D,DS -> remove nletons # improves SWARM for DS
  # filter1B: SD,S,SaD -> remove nletons # this will improve SWARM performance
  #                                      # apply this to SaD just in case
  #
  # 5: D,DS,SD,S,SaD -> export csv file of sequences (if denoised, they are ESV)
  # 6: D,DS -> export fasta DnoisEd & blank filt
  # 7: D -> create fasta files for taxonomic assignment
  # 8: DS,SD,S,SaD -> do SWARM
  #
  # filter2B: SD,S -> blank relative abundances filter AS MOTU
  # filter3B: SD,S -> min relative abundances filter WS MOTU
  # filter4B: SD,S -> remove nletons AS MOTUs # remove artefacts
  #
  # 9: DS,SD,S,SaD -> create csv files of MOTUs and ESV or seqs.
  # 10: DS,SD,S,SaD -> Also create fasta files for taxonomic assignment
  # 11: SD -> run DnoisE over the csv files
  #
  
  #####
  # this is the table of contents of the function
  # process                           | DS | D | S after D | SD | S
  # DnoisE fasta files                | X  | X |           |    |
  # cat all fasta                     | X  | X |     X     | X  | X
  # dereplicate sequences             | X  | X |     X     | X  | X
  # annotate new names                | X  | X |     X     | X  | X
  # blank relative abundances filter  | X  | X |           |    |
  #     AS ESV
  # min relative abundances filter    | X  | X |           |    |
  #     WS ESV
  # remove nletons                    | X  | X |     X     | X  | X
  #     AS ESV
  # export csv                        | X  | X |     X     | X  | X
  # export fasta DnoisEd & blank filt | X  | X |           |    |
  # generate fasta file               | X  |   |           |    |
  # SWARM                             | X  |   |     X     | X  | X
  # blank relative abundances filter  |    |   |           | X  | X
  #     AS MOTUs
  # min relative abundance            |    |   |           | X  | X
  #     WS MOTUs
  # remove nletons                    |    |   |           | X  | X
  #     AS MOTUs
  # generate fasta file               |    | X |     X     | X  | X
  # DnoisE csv file                   |    |   |           | X  |
  
  options(scipen = 999)
  
  #####
  # 0: define variables
  #####
  algorithm <- tolower(algorithm)
  swarm <- "swarm"
  dnoise <- "dnoise"
  
  suppressPackageStartupMessages(library(parallel))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))
  
  if (exists("lib") && is.null(experiment)) {
    # Use lib as experiment
    experiment <- lib
    # Print deprecation warning
    warning("The 'lib' argument is deprecated. Please use 'experiment' instead.")
  }
  
  
  # deifine output names
  outfile <- paste0(experiment,"_ODIN")# "_part_",sprintf("%02d",part),".fasta"
  if (algorithm == "swarm_dnoise" || algorithm == "swarm" ||
      algorithm == "dnoise_swarm") {
    # file with the reads of each MOTU in each sample
    outfile_MOTU <- paste0(outfile, "_counts.tsv")
    # file with the reads of each ESV in each sample and the MOTU it belongs to
    outfile_ESV <- paste0(outfile, "_ESV.tsv")
    # for algorithm swarm_dnoise this is the csv that will be denoised
    outfile_preDnoisE <- paste0(outfile, "_seqs_to_dnoise.tsv")
    # + _Adcorr_denoised_ratio_d.csv | _Adcorr_denoised_ratio_d.csv
    # intermediate input/ouput from SWARM
    fileswarm <- paste0(outfile, "_SWARM_output")
  }
  
  message("ODIN will first clear the battle field.")
  message("All obidms objects called *ODIN.obidms will be removed.")
  system("rm -r *ODIN.obidms", intern = TRUE, wait = TRUE)
  
  if (!(algorithm=="dnoise_swarm" || algorithm=="dnoise" || algorithm=="swarm_dnoise" || algorithm=="swarm")) {
    message("ERROR: algorithm has to be one of the following:\nDnoisE_SWARM\nSWARM_DnoisE\nSWARM\nDnoisE")
    quit()
  }
  run_entropy <- !is.logical(entropy)
  
  # checkpoint
  report_ODIN <- paste("ODIN was used to obtain meaningful units. In your",
                       "case you chose the ", algorithm, " algorithm.\n")
  
  #####
  # 1: D and DS -> denoise the fasta files
  #####
  
  # list files
  sample_files <- list.files(pattern=paste0("^",experiment,"_.*_HELA.fasta$"))
  sample_list <- gsub("_HELA.fasta", "", sample_files)
  
  if (file.exists("summary_HELA.RData")) {
    load("summary_HELA.RData")
    before_1_ODIN <- after_HELA # nolint: object_usage_linter.
  } else {
    before_1_ODIN <- mclapply(sample_list,function(file){
      output <- system(paste0("grep '>' ",file,"_HELA.fasta | wc -l"),
                       intern = TRUE, wait = TRUE)
      value <- as.numeric(output)
      return(data.frame(file=paste0(file,"_HELA.fasta"),
                        num_seqs=value))
    },mc.cores = cores)
    before_1_ODIN <- do.call(rbind, before_1_ODIN)
  }
  # keep not void files
  sample_files <- sample_files[sample_files %in%
                                 before_1_ODIN$file[before_1_ODIN$num_seqs > 0]]
  sample_list <- gsub("_HELA.fasta", "", sample_files)
  
  if (algorithm == "dnoise_swarm" || algorithm == "dnoise") {
    message("ODIN will denoise each sample file")
    if(run_dnoise || algorithm == "dnoise") {
      dnoise_fasta(before_1_ODIN, entropy, cores, sample_list,
                   alpha, min_reads_ESV, dnoise, run_entropy)
      if (run_entropy) {
        sample_files <- gsub("_HELA.fasta", "_ODIN_Adcorr_denoised_ratio_d.fasta",
                             sample_files)
      } else  {
        sample_files <- gsub("_HELA.fasta", "_ODIN_denoised_ratio_d.fasta",
                             sample_files)
      }
    } else {
      sample_files <- list.files(pattern = paste0("^",experiment,"_.*_ODIN_ESV.*.fasta$"))
      sample_list <- gsub("_ODIN_ESV.*.fasta", "", sample_files)
    }
  }
  
  # remove from sample list those samples with zero reads
  after_1 <- mclapply(sample_files,function(file){
    output <- system(paste0("grep '>' ",file," | wc -l"),
                     intern = TRUE, wait = TRUE)
    value <- as.numeric(output)
    return(data.frame(file=paste0(file),
                      num_seqs=value))
  },mc.cores = cores)
  after_1 <- do.call(rbind, after_1)
  sample_files <- sample_files[sample_files %in%
                                 after_1$file[after_1$num_seqs > 0]]
  sample_list <- gsub("_ODIN.*.fasta", "", sample_files)
  
  #####
  # 2: D,DS,SD,S,SaD -> cat all fasta
  #####
  
  # import fasta files into obidms and annotate with sample name
  samples_2_obidms(sample_list, sample_files, cores)
  
  # cat all samples into one obidms
  cat_samples_obidms(sample_list, experiment)
  
  #####
  # 3: D,DS,SD,S,SaD -> dereplicate
  #####
  
  dereplicate_obidms(experiment)
  
  #####
  # 4: D,DS,SD,S,SaD -> annotate new names
  #####
  
  # annotate
  annotate_obidms(experiment)
  
  ############################################
  annotate_samples(sample_list, sample_files, cores, experiment)
  cat_samples(sample_list, experiment)
  dereplicate_vsearch(experiment)
  rename_sequences(experiment)
  # export to csv and read to apply filters
  if (algorithm == "dnoise"){
    filetab <- paste0(experiment, "_ODIN_ESV.csv")
  } else{
    filetab <- paste0(experiment, "_ODIN_seqs.csv")
  }
  table_creation()
  ############################################

  # checkpoint
  output <- system(paste0("obi ls ", experiment, "_ODIN | grep 'Line count'"), intern = T, wait = T)
  values <- as.numeric(gsub(".*count: ", "", output))
  version <- gsub(".*# ", "", gsub(": Date.*", "", output))
  after_2_ODIN <- data.frame(algorithm = algorithm,
                             version = version,
                             num_seqs = values)
  
  
  # export to csv and read to apply filters
  if (algorithm == "dnoise"){
    filetab <- paste0(experiment, "_ODIN_ESV.csv")
  } else{
    filetab <- paste0(experiment, "_ODIN_seqs.csv")
  }
  
  system(paste0("obi export --tab-output --sep ','  -o ",
                filetab, " ",
                experiment, "_ODIN/seq_id"),
         intern = T, wait = T)
  seqs_abund <- read.csv(filetab, sep = ",", head = TRUE)
  names(seqs_abund) <- gsub("MERGED_sample.", "", names(seqs_abund))
  if(metadata_table == '') {
    metadata_table <- paste0(experiment, "_metadata.tsv")
  }
  metadata <- read.table(metadata_table,
                         sep = "\t", header = TRUE)
  
  #####
  # filter1A: D,DS -> blank relative abundances filter AS ESV
  #####
  #####
  # filter2A: D,DS -> min relative abundances filter WS ESV
  #####
  if((algorithm == "dnoise_swarm" && run_dnoise) || algorithm == "dnoise") {
    # filter1A: D,DS -> blank relative abundances filter AS ESV
    seqs_abund <- blank_correction(seqs_abund,
                                   blank_col, metadata, blank_tag,
                                   blank_relative)
    sample_cols <- grep("_sample_", names(seqs_abund))
    
    # filter2A: D,DS -> min relative abundances filter WS ESV
    seqs_abund <- min_rel_correction(seqs_abund, min_relative)
  }
  sample_cols <- grep("_sample_", names(seqs_abund))
  
  #####
  # filter3A: D,DS -> remove nletons # improves SWARM for DS
  #####
  #####
  # filter1B: SD,S,SaD -> remove nletons # this will improve SWARM performance
  #####
  if (!is.null(remove_singletons)) {
    # Message warning that remove singletons is not in use anymore
    warning(paste("The 'remove_singletons' argument is deprecated.",
                  "Please use 'min_reads_ESV' instead."))
    if(min_reads_ESV < 2) min_reads_ESV <- 2
  }
  if (min_reads_ESV > 1) {
    message(paste("ODIN will remove sequences below", min_reads_ESV, "reads"))
    seqs_abund[rowSums(seqs_abund[, sample_cols]) < min_reads_ESV,
               sample_cols] <- 0
  }
  
  #####
  # 5: D,DS,SD,S,SaD -> export csv file of sequences (if denoised, they are ESV)
  #####
  write.table(seqs_abund, filetab, row.names = FALSE, col.names = TRUE,
              sep = "\t", quote = FALSE)
  
  #####
  # 6: D,DS -> export fasta DnoisEd & blank filt
  #####
  if(algorithm == "dnoise" || (algorithm == "dnoise_swarm" && run_dnoise)) {
    for (sample in names(seqs_abund)[sample_cols]) {
      df_fasta <- seqs_abund[seqs_abund[, sample] > 0,
                             c("ID", sample, "NUC_SEQ")]
      id_fasta <- paste(">", df_fasta$ID, ";size=", df_fasta[, sample],
                        sep = "")
      fasta_to_print <- paste(id_fasta, df_fasta$NUC_SEQ, sep = "\n")
      original_name <- metadata$original_samples[metadata$mjolnir_agnomens ==
                                                   sample]
      writeLines(paste0(fasta_to_print,
                        collapse = "\n"),
                 paste0(sample, "_ODIN_ESV_", original_name, ".fasta"))
    }
    
    #####
    # 7: D -> create fasta files for taxonomic assignment
    #####
    if(algorithm == "dnoise"){
      seqs_abund <- seqs_abund[, c("ID", "NUC_SEQ")]
      seqs_abund <- paste(paste0(">", seqs_abund$ID), seqs_abund$NUC_SEQ, sep = "\n")
      writeLines(paste0(seqs_abund, collapse = "\n"), paste0(outfile, ".fasta"))
    }
  }
  
  
  
  #####
  # 8: DS,SD,S,SaD -> do SWARM and create csv files of MOTUs and ESV or seqs.
  #####
  
  if (algorithm == "swarm_dnoise" || algorithm == "swarm" ||
      algorithm == "dnoise_swarm") { 
    message("ODIN will cluster sequences into MOTUs with SWARM.")
    system(paste0("obi export --fasta-output --only-keys \"COUNT\" ",
                  experiment, "_ODIN/seq_id > ",outfile,".fasta ; ",
                  "sed -i 's/COUNT/size/g' ",outfile,".fasta ; ",
                  "sed -i 's/;//g' ",outfile,".fasta ; ",
                  "sed -E -i 's/(size=[0-9]*).*/\\1;/g' ", outfile, ".fasta ; ",
                  "sed -i 's/ /;/g' ",outfile, ".fasta "))
    system(paste0(swarm, " -d ", d, " -z -t ", cores,
                  " -o ", outfile, "_SWARM_output ",
                  " -s ", outfile, "_SWARM", d, "nc_stats ",
                  " -w ", outfile, "_SWARM_seeds.fasta ", outfile, ".fasta"),
           intern = TRUE, wait = TRUE)
    
    message("ODIN will recount abundances for every MOTU after Swarm.")
    
    # Read cluster list database
    message("1. ODIN is reading SWARM results...")
    swarm_db <- readLines(fileswarm)
    total_swarms <- length(swarm_db)
    message("2. ODIN has read ", total_swarms," total MOTUs.")
    
    
    clusters <- strsplit(swarm_db, "; ")
    message("4. ODIN is keeping only information of the sequences ",
            "that form each cluster.")
    clusters <- mclapply(clusters,function(x) {sub(";.*", "", x)},
                         mc.cores = cores)
    names(clusters) <- mclapply(clusters, function(x) x[[1]], mc.cores = cores)
    
    # Read counts database and keep only the needed clusters
    message("6. ODIN is reading the abundance database. ",
            "This could take Him a while, since He has just ",
            "one eye left, after all.")
    motu_seqs_names <- stack(clusters) %>% rename(ID = values, MOTU = ind)
    
    db <- read.table(filetab, sep = "\t", head = TRUE)
    numseqs <- nrow(db)
    db <- merge(motu_seqs_names, db, by = "ID")
    numseqs_reduced <- nrow(db)
    samples <- sum(grepl("sample", names(db)))
    message("7. ODIN finished reading the Database, which includes ", 
            numseqs, " total unique sequences and ", samples, " samples.\n",
            "ODIN kept only ", numseqs_reduced, " sequences for calculations.")
    
    message("8. ODIN will now calculate the number of reads in every sample ",
            "for each MOTU.")
    db_total <- split(db[, grepl("sample",names(db))], db$MOTU)
    db_total <- mclapply(db_total, function(x) {
      as.data.frame(t(as.matrix(c(COUNT = sum(x),
                                  colSums(x),
                                  CLUST_WEIGHT = dim(x)[1]))))},
      mc.cores = cores)
    db_total <- do.call(rbind, db_total)
    db_total <- cbind(data.frame(ID = rownames(db_total)), db_total)
    db_total <- merge(db_total, db[, grepl("ID|NUC_SEQ", names(db))], by = "ID")
    
    # order the columns
    col_order <- c("ID", "COUNT", "MOTU", names(db)[grepl("sample", names(db))],
                   "NUC_SEQ")
    db <- db[, col_order]
  }
  
  #####
  # filter2B: SD,S -> blank relative abundances filter AS MOTU
  #####
  #####
  # filter3B: SD,S -> min relative abundances filter WS MOTU
  #####
  #####
  # filter4B: SD,S -> remove nletons AS MOTUs # remove artefacts
  #####
  if(algorithm == "swarm_dnoise" || algorithm == "swarm") {
    sample_cols <- grep("_sample_", names(db_total))
    # filter2B: SD,S -> blank relative abundances filter AS MOTU
    db_total <- blank_correction(db_total,
                                 blank_col, metadata, blank_tag,
                                 blank_relative)
    
    # filter3B: SD,S -> min relative abundances filter WS MOTU
    db_total <- min_rel_correction(db_total, min_relative)
    
    sample_cols <- grep("_sample_", names(db_total))
    # filter4B: SD,S -> remove nletons AS MOTUs # remove artefacts
    if(min_reads_MOTU > 1) {
      message(paste("ODIN will remove sequences below", min_reads_MOTU, "reads"))
      db_total[rowSums(db_total[, sample_cols]) < min_reads_MOTU,
               sample_cols] <- 0
    }
  }
  
  
  #####
  # 9: DS,SD,S,SaD -> create csv files of MOTUs and ESV or seqs.
  #####
  #####
  # 10: D,DS,SD,S,SaD -> Also create fasta files for taxonomic assignment
  #####
  if (algorithm == "swarm_dnoise" || algorithm == "swarm" ||
      algorithm == "dnoise_swarm") {
    
    db <- db[db$MOTU %in% db_total$ID,]
    s_opt <- min(grep("sample",names(db)))
    z_opt <- max(grep("sample",names(db)))
    
    # print datasets
    # MOTUs
    write.table(db_total, outfile_MOTU, sep = "\t", quote = FALSE,
                row.names = FALSE)
    # ESV/seqs
    write.table(db, ifelse(algorithm == "swarm_dnoise",
                           outfile_preDnoisE, outfile_ESV),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # 9: DS,SD,S,SaD -> Also create fasta files for taxonomic assignment
    db_total <- db_total[, c("ID", "NUC_SEQ")]
    db_total <- paste(paste0(">", db_total$ID), db_total$NUC_SEQ, sep = "\n")
    writeLines(paste0(db_total, collapse = "\n"), paste0(outfile, ".fasta"))
    
  }
  
  #####
  # 11: SD -> run DnoisE over the csv files
  #####
  if (algorithm=="swarm_dnoise") { # 4ab
    message("ODIN will generate now a list of ESVs for every non-singleton MOTU, using DnoisE.")
    if(entropy[1] == "auto_dataset" || entropy[1] == "auto_sample") {
      message("for the SWARM_DnoisE algorithm entropy equals to entropy_sample or auto_dataset are equivalent")
      system(paste0(dnoise," --csv_input ",outfile_preDnoisE," -g"),intern = T, wait = T)
      entropies  <- read.csv(paste0(outfile_preDnoisE,"_entropy_values.csv"))
      entropy <- c(entropies[1,4:6],entropies[1,1])
    }
    if (is.logical(entropy)) {
      entropy <- ""
    } else {
      entropy <- paste0(" -y -e ",paste0(entropy[1:3],collapse = ",")," -m ",entropy[4])
    }
    
    if (run_entropy) {
      system(paste0(dnoise, " --csv_input ", outfile_preDnoisE, " ",
                    "--csv_output ", outfile_ESV, " ",
                    "-a ", alpha, " -c ", cores, 
                    " -n 'COUNT' -p 1 -q 'NUC_SEQ' ",
                    "-s ", s_opt, " -z ", z_opt, " ",
                    entropy, " -w 'MOTU' ; ",
                    "mv ", outfile_ESV, "_Adcorr_denoised_ratio_d.csv ",
                    outfile_ESV),
             intern = TRUE, wait = TRUE)
      remaining_ESV <- as.numeric(system(paste0("wc -l ", outfile_ESV,
                                                " |  cut -f1 -d ' ' "),
                                         intern = TRUE, wait = TRUE)) - 1
    } else  {
      system(paste0(dnoise," --csv_input ", outfile_preDnoisE, " ",
                    "--csv_output ", outfile_ESV, " ",
                    "-a ", alpha, " -c ", cores,
                    " -n 'COUNT' -p 1 -q 'NUC_SEQ' ",
                    "-s ", s_opt, " -z ", z_opt, " ",
                    "-w 'MOTU' ; ",
                    "mv ", outfile_ESV, "_denoised_ratio_d.csv ",
                    outfile_ESV),
             intern = TRUE, wait = TRUE)
      remaining_ESV <- as.numeric(system(paste0("wc -l ", outfile_ESV,
                                                " |  cut -f1 -d ' ' "),
                                         intern = TRUE, wait = TRUE)) - 1
    }
    after_4a_ODIN <- rbind(after_4a_ODIN,data.frame(version = c("remaining_ESV"),
                                                    value = c(remaining_ESV)))
    message("")
  }
  
  save(file = "summary_ODIN.RData",list = c(c("alpha","min_reads_MOTU","min_reads_ESV","algorithm","run_entropy","entropy","after_2_ODIN"),
                                            c("before_1_ODIN")[exists("before_1_ODIN")],c("after_4a_ODIN")[exists("after_4a_ODIN")]))
  if (remove_DMS) {
    system(paste0("rm -r *ODIN.obidms "), intern = TRUE, wait = TRUE)
  }
  message("ODIN is done.")
}

blank_correction <- function(seqs_abund, 
                             blank_col, metadata, blank_tag,
                             blank_relative){
  if(blank_col %in% names(metadata) &&
     sum(grepl(blank_tag, metadata[, blank_col])) > 0) {
    blank_columns <- grep(paste0(metadata$mjolnir_agnomens[grep(blank_tag,
                                                                metadata[,blank_col])], 
                                 collapse = "|"), names(seqs_abund))
    if(length(blank_columns) != 0) {
      neg_samples <- seqs_abund[, blank_columns]
      seqs_abund <- seqs_abund[, -blank_columns]
      message("The samples used in this steps as blanks are: ", 
              paste(metadata$mjolnir_agnomens[grep(blank_tag,metadata[,blank_col])],
                    collapse = ', '))
      if (length(blank_columns) == 1) {
        neg_reads <- neg_samples
      } else {
        neg_reads <- rowSums(neg_samples)
      }
      if (!is.data.frame(seqs_abund)) { 
        message("You have the honor to have found an error predicted by Adria. Tell him to solve it.")
      }
      sample_cols <- grep("_sample_", names(seqs_abund))
      sample_reads <- rowSums(seqs_abund[,sample_cols])
      # data_neg_filt_deleted <- seqs_abund[neg_reads/(sample_reads+neg_reads) > blank_relative,]
      seqs_abund <- seqs_abund[!neg_reads/(sample_reads+neg_reads) > blank_relative,]
      # sample_cols <- grep("_sample_", names(seqs_abund))
    }
  }
  return(seqs_abund)
}

min_rel_correction <- function(seqs_abund, min_relative){
  if(min_relative > 0) {
    sample_cols <- grep("_sample_", names(seqs_abund))
    message("ODIN will remove sequences with relative abundances below ", 
            min_relative)
    num_seqs <- colSums(seqs_abund[,sample_cols]>0)
    reads_seqs <- colSums(seqs_abund[,sample_cols])
    
    change_matrix <- apply(seqs_abund[,sample_cols], 2,
                           relabund, min_relative = min_relative)
    
    seqs_abund[, sample_cols][change_matrix] <- 0
    seqs_abund <- seqs_abund[rowSums(seqs_abund[,sample_cols]) > 0,]
    # num_seqs_after <- colSums(seqs_abund[,sample_cols]>0)
    # reads_seqs_after <- colSums(seqs_abund[,sample_cols]) 
    # after_3_ODIN <- data.frame(version = c("min_relative","num_seqs","num_seqs_reduced"),
    #                           value = c(min_relative,nrow(seqs_abund),nrow(seqs_abund)))
    # } else {
    # after_3_ODIN <- data.frame(version = c("min_relative","num_seqs","num_seqs_reduced"),
    #                           value = c(min_relative, nrow(seqs_abund), NA))
  }
  return(seqs_abund)
}

dnoise_fasta <- function(before_1_ODIN, entropy, cores, sample_list, alpha, min_reads_ESV, dnoise, run_entropy) {
  if(entropy[1] == "auto_dataset") {
    system(paste0("cat *_HELA.fasta >file_for_entropy.fasta ; ",
                  dnoise, " --fasta_input file_for_entropy.fasta -g"),
           intern = TRUE, wait = TRUE)
    entropies  <- read.csv("file_for_entropy.fasta_entropy_values.csv")
    entropy <- c(entropies[1, 4:6],entropies[1, 1])
  }
  for (file in sample_list) {
    seq_counts <- before_1_ODIN$num_seqs[before_1_ODIN$file == paste0(file, "_HELA.fasta")]
    if (seq_counts==1) {
      if (run_entropy) {
        system(paste0("cp ", file, "_HELA.fasta ",
                      file, "_ODIN_Adcorr_denoised_ratio_d.fasta "),
               intern = TRUE, wait = TRUE)
      } else  {
        system(paste0("cp ", file, "_HELA.fasta ",
                      file, "_ODIN_denoised_ratio_d.fasta "),
               intern = TRUE, wait = TRUE)
      }
    } else {
      if (!run_entropy) {
        entropy_file <- ""
      } else if(entropy[1]=="auto_sample") {
        entropy_file <- paste0(" -y -m ", entropy[2])
      } else {
        entropy_file <- paste0(" -y -e ", paste0(entropy[1:3],
                                                 collapse = ","),
                               " -m ", entropy[4])
      }
      
      message(paste("ODIN will denoise", file))
      command <- paste0(dnoise,
                        " --fasta_input ", file, "_HELA.fasta ",
                        "--fasta_output ", file, "_ODIN ",
                        "-a ", alpha, " -c ", cores,
                        entropy_file, " -r ", min_reads_ESV)
      system(command,
             intern = TRUE, wait = TRUE)
    }
  }
}


samples_2_obidms <- function(sample_list, sample_files, cores) {
  X <- NULL
  for (file in sample_list) {
    input_file <- sample_files[grep(file, sample_files)]
    X <- c(X,
           paste0("tempfile=$(mktemp); ",
                  "sed 's/;size/; COUNT/g' ",input_file, " > $tempfile ; ",
                  "obi import --fasta-input $tempfile ", file,"_ODIN/sample ; ",
                  "obi annotate -S sample:\"", gsub("^[a-zA-Z0-9]{4}_", "", file), "\" ", 
                  file,"_ODIN/sample  ", file,"_ODIN/sample_name ;
                  rm $tempfile"))
    
  }
  mclapply(X, function(x) system(x,intern = TRUE, wait = TRUE), mc.cores = cores)
}

cat_samples_obidms <- function(sample_list, experiment) {
  for (i in seq_along(sample_list)) {
    file <- sample_list[i]
    if (i == 1) {
      system(paste0("obi cat -c ",
                    file,"_ODIN/sample_name",
                    " ", experiment,"_ODIN/version",i),
             intern = TRUE, wait = TRUE)
    } else if (i==length(sample_list)) {
      system(paste0("obi cat -c ",
                    file, "_ODIN/sample_name",
                    " -c ", experiment, "_ODIN/version", c(i - 1),
                    " ", experiment, "_ODIN/samples ; ",
                    "obi rm ", experiment, "_ODIN/version", c(i - 1)),
             intern = TRUE, wait = TRUE)
    } else {
      system(paste0("obi cat -c ",
                    file, "_ODIN/sample_name",
                    " -c ", experiment, "_ODIN/version", c(i - 1),
                    " ", experiment, "_ODIN/version", i, " ; ",
                    "obi rm ", experiment, "_ODIN/version", c(i - 1)),
             intern = TRUE, wait = TRUE)
    }
  }
}

dereplicate_obidms <- function(experiment) {
  system(paste0("obi uniq --merge 'sample' ",
                experiment, "_ODIN/samples ",
                experiment, "_ODIN/samples_uniq"),
         intern = TRUE, wait = TRUE)
}

annotate_obidms <- function(experiment) {
  system(paste0("obi annotate --seq-rank ",
                experiment, "_ODIN/samples_uniq ",
                experiment, "_ODIN/seq_rank"),
         intern = TRUE, wait = TRUE)
  system(paste0("obi annotate --set-identifier ",
                "\'\"\'", experiment, "\'_%09d\" % sequence[\"seq_rank\"]\' ",
                experiment, "_ODIN/seq_rank ",
                experiment, "_ODIN/seq_id"),
         intern = TRUE, wait = TRUE)
}

############################################
# new part
annotate_samples <- function(sample_list, sample_files, cores, experiment) {
  X <- NULL
  for (file in sample_list) {
    input_file <- sample_files[grep(file, sample_files)]
    X <- c(X,
          paste0("sed 's/;$/;sample=",gsub(paste0("^", experiment, "_"), "", file),"/g' ",input_file,
                  " > ", file, "_ODIN_sample_annotated.fasta ; "))
    
  }
  mclapply(X, function(x) system(x,intern = TRUE, wait = TRUE), mc.cores = cores)
}
cat_samples <- function(sample_list, experiment) {
  system(paste0("cat *_ODIN_sample_annotated.fasta > ", experiment, "_ODIN_samples_joined.fasta; rm *_ODIN_sample_annotated.fasta"),
        intern = TRUE, wait = TRUE)
}

dereplicate_vsearch <- function(experiment){
  system(paste0("vsearch ",
                "--fastx_uniques ", experiment, "_ODIN_samples_joined.fasta ",
                "--sizeout ",
                "--fastaout ", experiment, "_ODIN_derep_seqs.fasta"),
        intern = TRUE, wait = TRUE)
}

Rcpp::sourceCpp(system.file("src/rename_fasta.cpp", package = "mjolnir"))
Rcpp::sourceCpp(system.file("src/seq2tab.cpp", package = "mjolnir"))

rename_fasta <- function(input_file, output_file, experiment) {
  .Call('_mjolnir_rename_fasta', PACKAGE = 'mjolnir', input_file, output_file, experiment)
}
seq2tab <- function(input_table_file, fasta_file, id_column) {
  .Call('_mjolnir_seq2tab', PACKAGE = 'mjolnir', input_table_file, fasta_file, id_column)
}

rename_sequences <- function(experiment) {
  input_file <- paste0(experiment, "_ODIN_derep_seqs.fasta")
  output_file <- paste0(experiment, "_ODIN_derep_seqs_renamed.fasta")
  rename_fasta(input_file, output_file, experiment)
}
table_creation <- function(experiment, filetab){
  system(paste0("vsearch ",
                "--search_exact ", experiment, "_ODIN_samples_joined.fasta ",
                "--db ", experiment, "_ODIN_derep_seqs_renamed.fasta ",
                "--otutab ", filetab),
        intern = TRUE, wait = TRUE)
  # Example usage of the C++ function
  fasta_file <- paste0(experiment, "_ODIN_derep_seqs_renamed.fasta")
  seq2tab(filetab, fasta_file, "ID")
}

relabund <- function(x,min_relative) if (sum(x)>0) x/sum(x) < min_relative else FALSE