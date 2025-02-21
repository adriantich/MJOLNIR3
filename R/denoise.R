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

fasta_denoiser <- function(experiment = NULL,
                           sample_df = NULL,
                           cores = 1,
                           alpha = 4,
                           entropy = c(0.47, 0.23, 1.02, 313),
                           min_reads_ESV = 2,
                           min_relative = 1 / 50000,
                           blank_relative = 0.1,
                           ...) {
  experiment <- 'ULOY'
  setwd("~/example_MJOLNIR_multiplexed_data_metafast/")
  sample_df <- read.table("matadata_new.tsv",
                          sep = "\t", header = TRUE)
  entropy = c(0.47, 0.23, 1.02, 313)
  cores = 1
  alpha = 4
  min_reads_ESV = 2
  min_relative = 1 / 50000
  blank_relative = 0.1
  
  options(scipen = 999)
  
  suppressPackageStartupMessages(library(parallel))
  #   suppressPackageStartupMessages(library(dplyr))
  #   suppressPackageStartupMessages(library(tidyr))
  
  #####
  # 0: define variables
  #####
  dnoise <- "dnoise"
  # deifine output names
  outfile <- paste0(experiment,"_denoised_cleaned")
  
  run_entropy <- !is.logical(entropy)
  
  #####
  # 1: denoise the fasta files
  #####
  
  # list files
  # check that fasta_file column exist
  if(!"fasta_file" %in% colnames(sample_df)) {
    stop("fasta_file names not found.")
  }
  # check that metadata has the column original_samples
  if(!"original_samples" %in% colnames(sample_df)) {
    message("original_samples names not found.")
    sample_df$original_samples <- sample_df$fasta_file
  }
  
  # check that all samples have sequences
  before_dnoising <- mclapply(sample_df$fasta_file,function(file){
    output <- system(paste0("grep '>' ",file," | wc -l"),
                     intern = TRUE, wait = TRUE)
    value <- as.numeric(output)
    return(data.frame(file=file,
                      num_seqs=value))
  },mc.cores = cores)
  before_denoising <- do.call(rbind, before_dnoising)
  
  # keep not void files
  sample_df <- sample_df[sample_df$fasta_file %in%
                           before_denoising$file[before_denoising$num_seqs > 0],]
  
  message("starting to denoise")
  sample_df$dnoise_output <- dnoise_function(sample_df, entropy, cores, sample_list,
                                             alpha, min_reads_ESV, dnoise, run_entropy)
  
  after_dnoising <- mclapply(sample_df$dnoise_output,function(file){
    output <- system(paste0("grep '>' ",file," | wc -l"),
                     intern = TRUE, wait = TRUE)
    value <- as.numeric(output)
    return(data.frame(file=file,
                      num_seqs=value))
  },mc.cores = cores)
  after_denoising <- do.call(rbind, after_dnoising)
  
  sample_df <- sample_df[sample_df$dnoise_output %in%
                           after_denoising$file[after_denoising$num_seqs > 0],]
  
  #####
  # 2: write sample label
  #####
  sample_df$annotated_output <- annotate_samples_generic(sample_df$original_samples,
                                                 sample_df$dnoise_output,
                                                 cores, experiment)
  
  #####
  # 3: cat all samples into one fasta
  #####
  cat_samples_generic(sample_df$annotated_output, experiment)
  
  #####
  # 4: dereplicate
  #####
  dereplicate_vsearch_generic(experiment)
  
  #####
  # 5: rename sequences
  #####
  rename_sequences_generic(experiment)
  
  #####
  # 6: export to csv
  #####
  filetab <- paste0(experiment, "_ESV.tsv")
  table_creation_generic(experiment, filetab)
  
  seqs_abund <- read.csv(filetab, sep = "\t", head = TRUE)
  
  #####
  # filter1A: blank relative abundances filter AS ESV
  #####
  seqs_abund <- blank_correction_generic(seqs_abund,
                                 blank_col = "BLANK", metadata = sample_df, blank_tag = TRUE,
                                 blank_relative = blank_relative)
  sample_cols <- c(names(seqs_abund) %in% sample_df$original_samples)
  #####
  # filter2A: min relative abundances filter WS ESV
  #####
  seqs_abund <- min_rel_correction_generic(seqs_abund, sample_cols, min_relative)
  
  
  #####
  # filter3A: remove nletons
  #####
  if (min_reads_ESV > 1) {
    message(paste("ODIN will remove sequences below", min_reads_ESV, "reads"))
    seqs_abund[rowSums(seqs_abund[, sample_cols]) < min_reads_ESV,
               sample_cols] <- 0
  }
  
  #####
  # 7: export csv file of sequences (if denoised, they are ESV)
  #####
  write.table(seqs_abund, filetab, row.names = FALSE, col.names = TRUE,
              sep = "\t", quote = FALSE)
  
  #####
  # 8: export fasta DnoisEd & blank filt
  #####
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
               paste0(original_name, "_denoised_cleaned.fasta"))
  }
  
  #####
  # 9: create fasta files for taxonomic assignment
  #####
  seqs_abund <- seqs_abund[, c("ID", "NUC_SEQ")]
  seqs_abund <- paste(paste0(">", seqs_abund$ID), seqs_abund$NUC_SEQ, sep = "\n")
  writeLines(paste0(seqs_abund, collapse = "\n"), paste0(outfile, ".fasta"))
  
  message("denoising is done.")
}

dnoise_function <- function(sample_df, entropy, cores, alpha, min_reads_ESV, dnoise, run_entropy) {
  if(entropy[1] == "auto_dataset") {
    system(paste0("cat ",paste(sample_df$fasta_file,collapse = " "),
                  " >file_for_entropy.fasta ; ",
                  dnoise, " --fasta_input file_for_entropy.fasta -g"),
           intern = TRUE, wait = TRUE)
    entropies  <- read.csv("file_for_entropy.fasta_entropy_values.csv")
    entropy <- c(entropies[1, 4:6],entropies[1, 1])
  }
  if (!run_entropy) {
    entropy_file <- ""
    output_files <- paste0(sample_df$original_samples,
                           "_denoised_ratio_d.fasta")
  } else if(entropy[1]=="auto_sample") {
    entropy_file <- paste0(" -y -m ", entropy[2])
    output_files <- paste0(sample_df$original_samples,
                           "_Adcorr_denoised_ratio_d.fasta")
  } else {
    entropy_file <- paste0(" -y -e ", paste0(entropy[1:3],
                                             collapse = ","),
                           " -m ", entropy[4])
    output_files <- paste0(sample_df$original_samples,
                           "_Adcorr_denoised_ratio_d.fasta")
  }
  for (file in seq_along(nrow(sample_df))) {   
    message(paste("ODIN will denoise", sample_df$fasta_file[file]))
    command <- paste0(dnoise,
                      " --fasta_input ", sample_df$fasta_file[file],
                      " --fasta_output ", sample_df$original_samples[file],
                      " -a ", alpha, " -c ", cores,
                      entropy_file, " -r ", min_reads_ESV)
    system(command,
           intern = TRUE, wait = TRUE)
  }
  return(output_files)
}

blank_correction_generic <- function(seqs_abund, 
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
  } else {
    message("No blanks found")
  }
  return(seqs_abund)
}

min_rel_correction_generic <- function(seqs_abund, sample_cols, min_relative){
  if(min_relative > 0) {
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

annotate_samples_generic <- function(sample_list, sample_files, cores, experiment) {
  X <- NULL
  for (file in sample_list) {
    input_file <- sample_files[grep(file, sample_files)]
    X <- c(X,
           paste0("sed 's/;$/;sample=",gsub(paste0("^", experiment, "_"), "", file),"/g' ",input_file,
                  " > ", file, "_sample_annotated.fasta ; "))
    
  }
  mclapply(X, function(x) system(x,intern = TRUE, wait = TRUE), mc.cores = cores)
  annotated_files <- paste0(sample_list, "_sample_annotated.fasta")
  return(annotated_files)
}
cat_samples_generic <- function(sample_list, experiment) {
  system(paste0("cat ",paste(sample_list, collapse = " ")," > ", experiment, "_samples_joined.fasta; rm ",paste(sample_list, collapse = " ")),
         intern = TRUE, wait = TRUE)
}

dereplicate_vsearch_generic <- function(experiment){
  system(paste0("vsearch ",
                "--fastx_uniques ", experiment, "_samples_joined.fasta ",
                "--sizeout ",
                "--fastaout ", experiment, "_derep_seqs.fasta ; ",
                " rm ", experiment, "_samples_joined.fasta"),
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

rename_sequences_generic <- function(experiment) {
  input_file <- paste0(experiment, "_derep_seqs.fasta")
  output_file <- paste0(experiment, "_derep_seqs_renamed.fasta")
  rename_fasta(input_file, output_file, experiment)
  system(paste0("rm ", input_file), intern = TRUE, wait = TRUE)
}
table_creation_generic <- function(experiment, filetab){
  system(paste0("vsearch ",
                "--search_exact ", experiment, "_samples_joined.fasta ",
                "--db ", experiment, "_derep_seqs_renamed.fasta ",
                "--otutab ", filetab),
         intern = TRUE, wait = TRUE)
  # Example usage of the C++ function
  fasta_file <- paste0(experiment, "_derep_seqs_renamed.fasta")
  seq2tab(filetab, fasta_file, "ID")
}

relabund <- function(x,min_relative) if (sum(x)>0) x/sum(x) < min_relative else FALSE