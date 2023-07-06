#' HELA: Hierarchical Elimination of Lurking Artifacts
#' 
#' This function uses the uchime_denovo algorithm implemented in VSEARCH to 
#' remove chimaeric sequences from the dataset. 
#' 
#' @details 
#' HELA works in a sample-by-sample basis. HELA will process all 
#' individual fasta files in the current folder 
#' matching the pattern EXPX_XXXX_sample_XXX.fasta being EXPX the acronym set 
#' by lib parameter. This allows for parallel 
#' computing, significantly decreasing calculation times. 
#' 
#' @param lib Character string. Acronym for the experiment. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. However they can be the same.
#' 
#' @param cores Numeric. Number of threads for parallel processing.
#' 
#' @param vsearchpath Character string specifying the PATH to the vsearch package.
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
#' # Enter number of cores to be used in parallel. 
#' cores <- 7
#' 
#' # Run RAN
#' mjolnir1_RAN(R1_filenames,cores,lib_prefixes,R1_motif="_R1",R2_motif="_R2")
#' 
#' # set experiment acronym
#' lib <- "ULOY"
#' 
#' # Run FREYJA
#' mjolnir2_FREYJA(lib_prefix = lib_prefixes,lib = lib,cores = cores,Lmin=299,Lmax=320)
#' 
#' # set the maximum number of cores possible
#' cores <- 16
#' 
#' # Run HELA
#' mjolnir3_HELA(lib,cores)


mjolnir3_HELA <- function(lib, cores, vsearchpath = NULL){

  suppressPackageStartupMessages(library(parallel))
  old_path <- Sys.getenv("PATH")
  if (is.null(vsearchpath)) vsearchpath <- "~/vsearch-2.22.1/bin/"
  Sys.setenv(PATH = paste(old_path, path.expand(vsearchpath), sep = ":"))
  sample_list <- gsub("_FREYJA_uniq.fasta","",list.files(pattern="^[a-zA-Z0-9]{4}_[a-zA-Z0-9]{4}_sample_[a-zA-Z0-9]{3}_FREYJA_uniq.fasta$"))

  message("HELA will remove chimaeras from each sample")
  X <- NULL
  for (i in sample_list) {
    X <- c(X,paste0("vsearch --uchime_denovo ",i,"_FREYJA_uniq.fasta ",
                    "--sizeout --minh 0.90 ",
                    "--nonchimeras ",i,"_HELA_nonchimeras.fasta ",
                    "--chimeras ",i,"_HELA_chimeras.fasta ",
                    "--uchimeout ",i,"_HELA_uchimeout.log"))
  }
  no_cores <- cores
  clust <- makeCluster(no_cores)
  clusterExport(clust, list("X","old_path","vsearchpath"),envir = environment())
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  stopCluster(clust)

  after_HELA <- mclapply(sample_list,function(file){
    output <- system(paste0("grep '>' ",file,"_HELA_nonchimeras.fasta | wc -l"),intern = T,wait = T)
    value <- as.numeric(output)
    return(data.frame(file=paste0(file,"_HELA_nonchimeras.fasta"),
                      num_seqs=value))
  },mc.cores = cores)
  after_HELA <- do.call("rbind",after_HELA)

  save(file = "summary_HELA.RData",list = c("after_HELA"))

  message("HELA is done.")
}


