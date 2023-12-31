#' FREYJA: Filtering of Reads, Enrollment, Yoke-reads Joining and Alignment
#' 
#' FREYJA will use OBITools3 commands to merge paired-end reads, demultiplex 
#' libraries into samples (if needed),trim primer sequences, filter by length, 
#' split sequences per sample and dereplicate within each sample.
#' 
#' @details 
#' In case the data is already demultiplexed and consist of individual fastq 
#' files for each sample, use the option demultiplexed=TRUE. When demultiplexed=TRUE, 
#' FREYJA will read the names of each individual R1 fastq files from a column 
#' in the LIBX_metadata.tsv file, called fastq_name_R1.
#' In the metadata table, each sample in the original_samples column must have 
#' a matching fastq_name_R1 and a matching mjolnir_agnomen (LIBX_sample_XXX).
#' When demultiplexed=TRUE, you must also specify the primer_F and primer_R 
#' sequences in the options input to FREYJA. COI Leray-XT primers are specified 
#' by default. Otherwise, when demultiplexed=FALSE, the primers information 
#' must be already written in the LIBX_ngsfilter.tsv files.
#' Until further optimization, please do not set more than 7 cores for the parallel
#' processing of the next step, FREYJA.
#' 
#' @param lib_prefixes Character vector. Acronym for each sequencing library. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. The latter will be required in following steps. However 
#' they can be the same.
#' 
#' @param cores Numeric. Number of threads for parallel processing.
#' 
#' @param Lmin Numeric. Minimum bp length for a sequence to be accepted.
#' 
#' @param Lmax Numeric. Maximum bp length for a sequence to be accepted.
#' 
#' @param lib Character string. Acronym for the experiment. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. However they can be the same.
#' 
#' @param fastq_ouput Logical. If TRUE, a fastq file for each sample with the 
#' demultiplexed and quality filtered sequences will be retrieved. This is useful
#'  if you wish to publish your data in public repositories.
#'  
#' @param score_obilign Numeric. Minimum quality threshold to retain a sequence 
#' in the qualiy filtering after pairalignment.
#' 
#' @param demultiplexed Logical. If TRUE, the data has been already demultiplexed. 
#' otherwise this has to be set to FALSE. See details for further information.
#' 
#' @param primer_F Character string of the Forward primer. Necessary when 
#' demultiplexed=TRUE.
#' 
#' @param primer_R Character string of the Reverse primer. Necessary when 
#' demultiplexed=TRUE.
#' 
#' @param R1_motif Character string that distinguish the forward line file from
#' the reverse.
#' 
#' @param R2_motif Character string that distinguish the reverse line file from
#' the forward.
#' 
#' @param obipath Character string specifying the PATH to the obi binary.
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

mjolnir2_FREYJA <- function(lib_prefix="",cores=1,Lmin=299,Lmax=320,lib="EXPX", #fasta_output=F,
                            fastq_output=T,score_obialign=40,
                            demultiplexed=F,primer_F="GGWACWRGWTGRACWNTNTAYCCYCC",primer_R="TANACYTCNGGRTGNCCRAARAAYCA",
                            R1_motif="_R1",R2_motif="_R2",obipath=NULL,remove_DMS=T, run_on_tmp=F){

  message("FREYJA will do paired-end alignment, demultiplexing and length filter.")
  suppressPackageStartupMessages(library(parallel))
  no_cores <- cores*length(lib_prefix)
  old_path <- Sys.getenv("PATH")
  if (is.null(obipath)) obipath <- "~/obi3-env/bin/"
  obipath <- path.expand(obipath)
  Sys.setenv(PATH = paste(old_path, path.expand(obipath), sep = ":"))
  if (is.null(lib)) {
    message("lib can not be NULL, otherwise HELA won't find the files")
    exit()
  }

  if (run_on_tmp) {
    tmp = "/tmp/"
  } else {
    tmp = ""
  }
  
  import_fastqs <- NULL
  count_seqs <- NULL
  rest_of_commands <- NULL
  libslist <- NULL
  message("FREYJA will first clear the battle field.")
  message("Any directory or file containing the word FREYJA will be removed.")
  system("rm -r *FREYJA*",intern = T, wait = T)
  if (!demultiplexed){
    
    for (i in 1:cores) for (j in 1:length(lib_prefix)) {
      formatted_i <- sprintf("%02d",i)
      import_fastqs <- c(import_fastqs,paste0(
        "obi import --fastq-input ",lib_prefix[j],"_R2_part_",formatted_i,".fastq ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/reads2 ; ",
        "obi import --fastq-input ",lib_prefix[j],"_R1_part_",formatted_i,".fastq ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/reads1 ; "))
      count_seqs <- c(count_seqs, paste0(tmp, lib_prefix[j],"_",formatted_i,"_FREYJA"))
      rest_of_commands <- c(rest_of_commands,paste0(
        "obi import --ngsfilter-input ngsfilter_",lib_prefix[j],".tsv ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/ngsfile ; ",
        "obi alignpairedend -R ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/reads2 ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/reads1 ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/aligned_seqs ; ",
        ifelse(remove_DMS,paste0("obi rm ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/reads1 ; obi rm ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/reads2 ; "),""),
        "obi grep -p \"sequence[\'score\'] > ",score_obialign,"\" ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/aligned_seqs ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/good_seqs ; ",
        ifelse(remove_DMS,paste0("obi rm ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/aligned_seqs ; "),""),
        "obi ngsfilter -t ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/ngsfile -u ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/unidentified_seqs ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/good_seqs ",tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/identified_seqs ; ",
        ifelse(remove_DMS,paste0("obi rm ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/good_seqs ; "),""),
        "obi grep -p \"len(sequence)>",Lmin," and len(sequence)<",Lmax," and sequence[\'forward_tag\']!=None and sequence[\'reverse_tag\']!=None\" -S \"^[ACGT]+$\" ",tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/identified_seqs ",tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/filtered_seqs"))
      libslist <- paste0(libslist,paste0("-c ", tmp, lib_prefix[j],"_",formatted_i,"_FREYJA/filtered_seqs "))
    }
  } else {
    metadata <- read.table(paste0(lib,"_metadata.tsv"),sep="\t",header=T)
    fastqR1_list <- metadata$fastq_name_R1
    agnomens <-  metadata$mjolnir_agnomens
    before_FREYJA <- mclapply(fastqR1_list, function(prefix){
        return(data.frame(file=prefix,
                          num_seqs=as.numeric(system(paste0("grep '>' ",prefix," | wc -l"),intern = T,wait = T))))
      },mc.cores = cores)
    # Create ngsfilter files
    for (ag in agnomens) writeLines(paste(lib,ag,"None:None",primer_F,primer_R,sep="\t"),paste0("ngsfilter_",ag,".tsv"))
    # Create obitool commands
    for (i in 1:length(agnomens)) {
      print(fastqR1_list[i])
      import_fastqs <- c(import_fastqs,paste0(
        "obi import --fastq-input ",gsub(R1_motif,R2_motif,fastqR1_list[i]), " ", tmp, lib,"_",agnomens[i],"_FREYJA/reads2 ; ",
        "obi import --fastq-input ",fastqR1_list[i], " ", tmp, lib,"_",agnomens[i],"_FREYJA/reads1 ; "))
      count_seqs <- c(count_seqs, paste0(tmp, lib,"_",agnomens[i],"_FREYJA"))
      rest_of_commands <- c(rest_of_commands,paste0(
        "obi import --ngsfilter-input ngsfilter_",agnomens[i],".tsv ", tmp, lib,"_",agnomens[i],"_FREYJA/ngsfile ; ",
        "obi alignpairedend -R ", tmp, lib,"_",agnomens[i],"_FREYJA/reads2 ", tmp, lib,"_",agnomens[i],"_FREYJA/reads1 ", tmp, lib,"_",agnomens[i],"_FREYJA/aligned_seqs ; ",
        ifelse(remove_DMS,paste0("obi rm ", tmp, lib,"_",agnomens[i],"_FREYJA/reads1 ; obi rm ", tmp, lib,"_",agnomens[i],"_FREYJA/reads2 ; "),""),
        "obi grep -p \"sequence[\'score\'] > ",score_obialign,"\" ", tmp, lib,"_",agnomens[i],"_FREYJA/aligned_seqs ", tmp, lib,"_",agnomens[i],"_FREYJA/good_seqs ; ",
        ifelse(remove_DMS,paste0("obi rm ", tmp, lib,"_",agnomens[i],"_FREYJA/aligned_seqs ; "),""),
        "obi ngsfilter --no-tags -t ", tmp, lib,"_",agnomens[i],"_FREYJA/ngsfile -u ", tmp, lib,"_",agnomens[i],"_FREYJA/unidentified_seqs ", tmp, lib,"_",agnomens[i],"_FREYJA/good_seqs ", tmp, lib,"_",agnomens[i],"_FREYJA/identified_seqs ; ",
        ifelse(remove_DMS,paste0("obi rm ", tmp, lib,"_",agnomens[i],"_FREYJA/good_seqs ; "),""),
        "obi grep -p \"len(sequence)>",Lmin," and len(sequence)<",Lmax,"\" -S \"^[ACGT]+$\" ", tmp, lib,"_",agnomens[i],"_FREYJA/identified_seqs ", tmp, lib,"_",agnomens[i],"_FREYJA/filtered_seqs ; "))
    }
  }
  # clust <- makeCluster(no_cores)
  # clusterExport(clust, list("X","old_path","obipath"),envir = environment())
  # clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))})
  
  # run the commands in sets so in case of demultiplexed, if you have more samples
  # than cores, you can allow that they don't get saturated. However, I think with 
  # mclapplay this should not happen. # at the end I decide to do it without separating it.
  # nu_sets <- ceiling(length(import_fastqs)/no_cores)
  # for (set_i in 1:nu_sets) {
  #   from <- (set_i-1)*no_cores + 1
  #   if (set_i == nu_sets) {
  #     to <- length(import_fastqs)
  #   } else {
  #     to <- (set_i)*no_cores
  #   }
  #   parLapply(clust,import_fastqs[from:to], function(x) system(x,intern=T,wait=T)) # commented to run with mclapply
  # }
  mclapply(import_fastqs, function(x) system(x,intern=T,wait=T), mc.cores = cores)
  # here the number of lines correspond to the number of sequences. This way is faster
  # than doing this with grep > file | wc -l
  before_FREYJA <- mclapply(count_seqs,function(i){
      return(data.frame(file=i,
                        num_seqs=as.numeric(gsub(".*count: ","",system(paste0("obi ls ", i, " | grep reads1"),intern = T,wait = T), perl = T))))
    },mc.cores = cores)
  # now I run the rest of the commands.
  mclapply(rest_of_commands, function(x) system(x,intern=T,wait=T), mc.cores = cores)
  
  # stopCluster(clust)
  # If not demultiplexed, then join all parts into a joined file and then split it into samples
  if (!demultiplexed){
    message("FREYJA is joining filtered reads into a single file.")
    # this step is necessary to join all sets of sequences to split them into different files
    system(paste0("obi cat ",libslist," ", tmp, lib,"_FREYJA/concatenated"),intern=T,wait=T)
    if (remove_DMS) {
      system(paste0("rm -r ",gsub("-c ","",gsub("/filtered_seqs",".obidms",libslist))),intern=T,wait=T)
    }
    message("Sequence sets concatenated")
    message("FREYJA will create individual files for each sample in the dms directory, this could take a while..")
    # here the concatenated file is split into different samples
    system(paste0("obi split -t \"sample\" -p ",lib,"_ ", tmp, lib,"_FREYJA/concatenated"),intern=T,wait=T)
    # in order to make the next steps parallelizable it is nesary to export each file into a new dms
    files <- system(paste0("obi ls ", tmp, lib,"_FREYJA | cut -f1 -d ':' | cut -f4 -d ' ' | grep 'sample' "),intern=T,wait=T)
    files <- files[grep(lib,files)]
    for (file in files) {
      system(paste0("obi cat -c ",
                    tmp, lib,"_FREYJA/",file, " ",
                    tmp, file,"_FREYJA/filtered_seqs "),intern=T,wait=T)
    }
    if (remove_DMS) {
      system(paste0("rm -r ",tmp, lib,"_FREYJA.obidms "),intern=T,wait=T)
    } else if (run_on_tmp) {
      system(paste0("mv ",tmp, lib,"_FREYJA.obidms ."),intern=T,wait=T)
    }
  }
  # obi uniq vas performed in HELA in previous versions but now is computed here
  files <- list.dirs(recursive = F)
  files <- files[grepl("sample",files)&grepl("FREYJA.obidms",files)]
  files <- gsub("./","",gsub(".obidms","",files))
  X <- NULL
  if (fastq_output) {
    print("Exporting fastq files and dereplicating, this could take a while")
    
  }
  for (file in files) {
    X <- c(X,paste0(ifelse(fastq_output,paste0("obi export --fastq-output -o ",
                                               file,".fastq ",
                                               tmp, file,"/filtered_seqs ; "), ""),
                    "obi uniq ", tmp, file,"/filtered_seqs ",file,"/uniq ; ",
                    "obi export --fasta-output --only-keys \"COUNT\" ", tmp, file,"/uniq > ",file,"_uniq.fasta ; ",
                    "sed -i 's/COUNT/size/g' ",file,"_uniq.fasta ; ",
                    "sed -i 's/;//g' ",file,"_uniq.fasta ; ",
                    "sed -E -i 's/(size=[0-9]*).*/\\1;/g' ",file,"_uniq.fasta ; ",
                    "sed -i 's/ /;/g' ",file,"_uniq.fasta "))
  }
  clust <- makeCluster(no_cores)
  clusterExport(clust, list("X","old_path","obipath"),envir = environment())
  clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))})
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  stopCluster(clust)


  after_FREYJA <- mclapply(files,function(file){
    output <- system(paste0("obi ls ",tmp, file," | grep 'filtered_seqs\\|uniq'"),intern = T,wait = T)
    values <- as.numeric(gsub(".*count: ","",output))
    if (remove_DMS) {
      system(paste0("rm -r ",tmp, file,".obidms "),intern=T,wait=T)
    } else if (run_on_tmp) {
      system(paste0("mv ",tmp, file,".obidms ."),intern=T,wait=T)
    }

    return(data.frame(file=file,
                      version=c("filtered sequences","uniq sequences"),
                      num_seqs=values))
  },mc.cores = cores)

  variables_FREYJA <- data.frame(variable=c("cores","Lmin","Lmax","lib",#"fasta_output",
                                            "fastq_output","score_obialign","demultiplexed","primer_F","primer_R"),
                                 value=c(cores,Lmin,Lmax,lib,#fasta_output,
                                         fastq_output,score_obialign,demultiplexed,primer_F,primer_R))

  save(file = "summary_FREYJA.RData",list = c("before_FREYJA","after_FREYJA","variables_FREYJA"))

  # if required, export to fasta or fastq. These options is to return the files without joining amplicons with obi uniq
  # if (fasta_output) {
  #   files <- list.dirs(recursive = F)
  #   files <- files[grepl("sample",files)&grepl("FREYJA.obidms",files)]
  #   files <- gsub("./","",gsub(".obidms","",files))
  #   X <- NULL
  #   print("Exporting fasta files, this could take a while.")
  #   for (file in files) {
  #     X <- c(X,paste0("obi export --fasta-output ",
  #                     file,"/filtered_seqs > ",
  #                     file,".fasta "))
  #   }
  #   
  #   clust <- makeCluster(no_cores)
  #   clusterExport(clust, list("X","old_path","obipath"),envir = environment())
  #   clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))})
  #   parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  #   stopCluster(clust)
  # 
  # }


  message("FREYJA is done.")
}
