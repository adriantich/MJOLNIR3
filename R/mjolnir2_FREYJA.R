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
#' @param score_obilign Numeric. Minimum quality threshold to retain a sequence
#' in the qualiy filtering after pairalignment.
#'
#' @param R1_motif Character string that distinguish the forward line file from
#' the reverse.
#'
#' @param R2_motif Character string that distinguish the reverse line file from
#' the forward.
#'
#' @param remove_DMS Logical. If TRUE, it will delete all obidms objects that are
#' created during the process. This can save a lot of hard disk space. The FALSE
#' option is useful for developing and debugging.
#'
#' @param run_on_tmp Logical. If TRUE, the obidms objects will be created in
#' the /tmp location. This increases the speed as the communication within the
#' processor and the object that is being edited all the time is faster. However,
#' this method will consume much of the /tmp memory and it is recommended to have
#' three to four times the memory available in the /tmp directory than the original
#' forward files and remove_DMS=T
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
                            score_obialign = 40,
                            R1_motif = "_R1", R2_motif = "_R2", remove_DMS = T, 
                            run_on_tmp = F, ...) {
  
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

  if (run_on_tmp) {
    tmp <- "/tmp/"
  } else {
    tmp <- ""
  }

  filtering_commands <- NULL
  message("FREYJA will first clear the battle field.")
  message("Any directory or file containing the word FREYJA will be removed.")
  system("rm -r *FREYJA*", intern = TRUE, wait = TRUE)

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
  agnomens <- agnomens[to_retain$num_seqs > 0]
  # Create obitool commands
  for (i in seq_along(agnomens)) {
    print(fastqR1_list[i])
    filtering_commands <- c(filtering_commands,paste0(
      "obi import --fastq-input ",fastqR2_list[i], " ", tmp, experiment,"_",agnomens[i],"_FREYJA/reads2 ; ",
      "obi import --fastq-input ",fastqR1_list[i], " ", tmp, experiment,"_",agnomens[i],"_FREYJA/reads1 ; ",
      "obi alignpairedend -R ", tmp, experiment,"_",agnomens[i],"_FREYJA/reads2 ", tmp, experiment,"_",agnomens[i],"_FREYJA/reads1 ", tmp, experiment,"_",agnomens[i],"_FREYJA/aligned_seqs ; ",
      ifelse(remove_DMS,paste0("obi rm ", tmp, experiment,"_",agnomens[i],"_FREYJA/reads1 ; obi rm ", tmp, experiment,"_",agnomens[i],"_FREYJA/reads2 ; "),""),
      "obi grep -p \"sequence[\'score\'] > ",score_obialign,"\" ", tmp, experiment,"_",agnomens[i],"_FREYJA/aligned_seqs ", tmp, experiment,"_",agnomens[i],"_FREYJA/good_seqs ; ",
      ifelse(remove_DMS,paste0("obi rm ", tmp, experiment,"_",agnomens[i],"_FREYJA/aligned_seqs ; "),""),
      "obi grep -p \"len(sequence)>",Lmin," and len(sequence)<",Lmax,"\" -S \"^[ACGT]+$\" ", tmp, experiment,"_",agnomens[i],"_FREYJA/good_seqs ", tmp, experiment,"_",agnomens[i],"_FREYJA/filtered_seqs ; "))
  }

  # run all commands
  mclapply(filtering_commands, function(x) system(x,intern=T,wait=T), mc.cores = cores)

  # obi uniq vas performed in HELA in previous versions but now is computed here
  files <- list.dirs(recursive = F)
  files <- files[grepl("sample",files)&grepl("FREYJA.obidms",files)]
  files <- gsub("./","",gsub(".obidms","",files))
  X <- NULL

  for (file in files) {
    X <- c(X,paste0("obi uniq ", tmp, file,"/filtered_seqs ",file,"/uniq ; ",
                    "obi export --fasta-output --only-keys \"COUNT\" ", tmp, file,"/uniq > ",file,"_uniq.fasta ; ",
                    "sed -i 's/COUNT/size/g' ",file,"_uniq.fasta ; ",
                    "sed -i 's/;//g' ",file,"_uniq.fasta ; ",
                    "sed -E -i 's/(size=[0-9]*).*/\\1;/g' ",file,"_uniq.fasta ; ",
                    "sed -i 's/ /;/g' ",file,"_uniq.fasta "))
  }

  mclapply(X, function(x) system(x,
                                 intern = TRUE, wait = TRUE),
           mc.cores = cores)

  after_FREYJA <- mclapply(files,function(file){
    output <- system(paste0("obi ls ", tmp, file, " | grep 'filtered_seqs\\|uniq'"),
                     intern = TRUE, wait = TRUE)
    values <- as.numeric(gsub(".*count: ", "", output))
    if (remove_DMS) {
      system(paste0("rm -r ",tmp, file,".obidms "),
             intern = TRUE, wait = TRUE)
    } else if (run_on_tmp) {
      system(paste0("mv ",tmp, file,".obidms ."),
             intern = TRUE, wait = TRUE)
    }

    return(data.frame(file=file,
                      version=c("filtered sequences","uniq sequences"),
                      num_seqs=values))
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
