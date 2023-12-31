#' RAN: Reads Allotment in N portions
#' 
#' RAN function will prepare the FASTQ raw data for parallel processing.
#' 
#' @details 
#' RAN will be run only if your data consist of multiplexed libraries, 
#' which will be split into aliquote parts, to be processed by FREYJA. Until 
#' further optimization, please do not set more than 7 cores for the parallel
#' processing of the next step, FREYJA.
#' 
#' @param R1_filenames Character vector with the names of the forward fastq or 
#' fastq.gz files.
#' 
#' @param cores Numeric. Number of parts into which the input files will be split
#' for the parallel processing of the FREYJA function.
#' 
#' @param lib_prefixes Character vector. Acronym for each sequencing library. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. The latter will be required in following steps. However 
#' they can be the same.
#' 
#' @param R1_motif Character string that distinguish the forward line file from
#' the reverse.
#' 
#' @param R2_motif Character string that distinguish the reverse line file from
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
#' # Enter number of cores to be used in parallel. 
#' cores <- 7
#' 
#' mjolnir1_RAN(R1_filenames,cores,lib_prefixes,R1_motif="_R1",R2_motif="_R2")

mjolnir1_RAN <- function(R1_filenames="",cores=1,lib_prefixes="",R1_motif="_R1",R2_motif="_R2"){

  suppressPackageStartupMessages((library(parallel)))
  suppressPackageStartupMessages((library(stringr)))
  message(paste0("RAN will split initial FASTQ files in ",cores," fragments each."))
  filelist <- NULL
  outfilelist <- NULL
  old_path <- Sys.getenv("PATH")

  # here the name of R1 and R2 of all files
  for (file in R1_filenames) filelist <- c(filelist,file,gsub(R1_motif,R2_motif,file))

  # ouput file names using lib_prefixes
  for (prefix in lib_prefixes) outfilelist <- c(outfilelist,paste0(prefix,"_R1_part"),paste0(prefix,"_R2_part"))

  no_cores <- length(filelist)

  # unzip files if needed
  filelist <- unlist(mclapply(filelist,unzip_files,old_path = old_path,mc.cores = cores))

  # distribute into <cores> files with new names with prefixes
  mclapply(1:length(filelist),split_lines,old_path = old_path,filelist = filelist,outfilelist = outfilelist,cores=cores,mc.cores = cores)

  message("Splitting done.")
}

unzip_files <- function(file,old_path=""){
  if (grepl(".gz",file)) {
    if (file.exists(file)) {
      # run commands to unzip from the system
      Sys.setenv(PATH = old_path)
      system(paste0("gzip -dc ",file, " > ",str_remove(file,".gz")),intern=T,wait=T)
      file <- str_remove(file,".gz")
    } else {
      file <- str_remove(file,".gz")
    }
  }
  return(file)
}

vector_to_split <- function(core_num,total_cores,num_seqs){
  j <- core_num
  cores <- total_cores
  # get the lines on infile for each outfile
  lines_vector <- sort(c(seq(1+(j-1)*4,num_seqs*4,((cores-1)*4+4)),
                         seq(2+(j-1)*4,num_seqs*4,((cores-1)*4+4)),
                         seq(3+(j-1)*4,num_seqs*4,((cores-1)*4+4)),
                         seq(4+(j-1)*4,num_seqs*4,((cores-1)*4+4))))
  return(lines_vector)
}

split_lines <- function(num,old_path = old_path,filelist = filelist,outfilelist = outfilelist,cores=cores){
  file <- filelist[num]
  outfile <- outfilelist[num]
  Sys.setenv(PATH = old_path)

  # get the number of lines for each file
  num_lines <- system(paste0("wc -l ",file," |  cut -f1 -d ' ' "),intern=T,wait=T)
  num_lines <- as.numeric(num_lines)
  num_seqs <- num_lines/4

  # get vectors of lines that will end in each file
  lines_to_split <- lapply(1:cores, vector_to_split, total_cores = cores, num_seqs = num_seqs)

  # read file
  file_read <- readLines(file)

  # write each file
  for (i in 1:cores) {
    writeLines(file_read[lines_to_split[[i]]],
               paste0(outfile,"_",sprintf("%02d",i),".fastq"))
  }
}
