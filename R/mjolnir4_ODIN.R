#' ODIN: OTU Delimitation Inferred by Networks
#' 
#' ODIN performs MOTU clustering and/or denoising. It is one of the main steps of MJOLNIR.
#' 
#' @details 
#' The function mjolnir4_ODIN() uses the four different strategies to delimit 
#' MOTUs and/or ESVs. This strategies are set with the algorithm parameter: 
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
#' @param lib Character string. Acronym for the experiment. This
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
#' have to be retained.
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
#' @param obipath Character string specifying the PATH to the obi binary.
#' 
#' @param python_packages  Character string specifying the PATH to where the 
#' Python3 packages are stored.
#' 
#' @param swarmpath Character string specifying the PATH to the SWARM program.
#' 
#' @param dnoise_path Character string specifying the PATH to the DnoisE program.
#' 
#' @param remove_singletons Logical. If TRUE this will remove the sequences that,
#' after dereplication when joining all samples (those will be denoised for when 
#' algorithm="DnoisE_SWARM" or algorithm="DnoisE"), have only one read.
#' 
#' @param remove_DMS Logical. If TRUE, it will delete all obidms objects that are
#' created during the process. This can save a lot of hard disk space. The FALSE 
#' option is useful for developing and debugging.
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
#' 
#' # Run ODIN
#' mjolnir4_ODIN(lib,cores,d=13,min_reads_MOTU=2,min_reads_ESV=2,alpha=5,entropy=c(0.47,0.23,1.02,313), algorithm="DnoisE_SWARM", remove_singletons = TRUE)


mjolnir4_ODIN <- function(lib,cores,d=13,min_reads_MOTU=2,min_reads_ESV=2,alpha=5,entropy=c(0.47,0.23,1.02,313),#entropy=F/c("auto_sample",313)/c("auto_dataset")
                          algorithm="DnoisE_SWARM",obipath="",python_packages="", swarmpath=NULL, dnoise_path=NULL, 
                          remove_singletons = TRUE,remove_DMS=T){
  #####
  # 0: define variables
  #####
  algorithm = tolower(algorithm)
  if (is.null(dnoise_path)) dnoise_path <- "~/DnoisE/src/"
  dnoise_path <- path.expand(dnoise_path)
  if (is.null(swarmpath)) swarmpath <- "~/swarm/bin/"
  swarmpath <- path.expand(swarmpath)
  swarm <- paste0(swarmpath,"swarm")
  if (is.null(obipath)) obipath <- "~/obi3-env/bin/"
  obipath <- path.expand(obipath)
  dnoise <- paste0("python3 ", dnoise_path, "DnoisE.py") # Change this to where the Dnoise executable is
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(python_packages,old_path, obipath, dnoise_path, sep = ":"))

  suppressPackageStartupMessages(library(parallel))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))

  if (!(algorithm=="dnoise_swarm" | algorithm=="dnoise" | algorithm=="swarm_dnoise" | algorithm=="swarm")) {
    message("ERROR: algorithm has to be one of the following:\nDnoisE_SWARM\nSWARM_DnoisE\nSWARM\nDnoisE")
    quit()
  }
  run_entropy <- !is.logical(entropy)

  #####
  # 1: D and DS -> denoise the fasta files
  #####

  # inputs
  sample_list <- gsub("_HELA_nonchimeras.fasta","",list.files(pattern="^[a-zA-Z0-9]{4}_[a-zA-Z0-9]{4}_sample_[a-zA-Z0-9]{3}_HELA_nonchimeras.fasta$"))
  # ouputs
  # file"_ODIN_Adcorr_denoised_ratio_d.fasta | file"_ODIN_denoised_ratio_d.fasta
  if (algorithm=="dnoise_swarm" | algorithm=="dnoise") {
    message("ODIN will denoise each sample file")
    if (file.exists("summary_HELA.RData")){
      load("summary_HELA.RData")
      before_1_ODIN <- after_HELA
    } else {
      before_1_ODIN <- mclapply(sample_list,function(file){
        output <- system(paste0("grep '>' ",file,"_HELA_nonchimeras.fasta | wc -l"),intern = T,wait = T)
        value <- as.numeric(output)
        return(data.frame(file=paste0(file,"_HELA_nonchimeras.fasta"),
                          num_seqs=value))
      },mc.cores = cores)
    }
    if(entropy[1]=="auto_dataset") {
      system(paste0("cat *_HELA_nonchimeras.fasta >file_for_entropy.fasta ; ",
                    dnoise," --fasta_input file_for_entropy.fasta -g"),intern = T, wait = T)
      entropies  <- read.csv("file_for_entropy.fasta_entropy_values.csv")
      entropy <- c(entropies[1,4:6],entropies[1,1])
    }
    for (file in sample_list) {
      seq_counts <- before_1_ODIN$num_seqs[before_1_ODIN$file==paste0(file,"_HELA_nonchimeras.fasta")]
      if (seq_counts==1) {
          if (run_entropy) {
            system(paste0("cp ",file,"_HELA_nonchimeras.fasta ",file,"_ODIN_Adcorr_denoised_ratio_d.fasta "),intern = T, wait = T)
          } else  {
            system(paste0("cp ",file,"_HELA_nonchimeras.fasta ",file,"_ODIN_denoised_ratio_d.fasta "),intern = T, wait = T)
          }
        } else {
          if (is.logical(entropy)) {
            entropy_file <- ""
          } else if(entropy[1]=="auto_sample") {
            entropy_file <- paste0(" -y -m ",entropy[2])
          } else {
            entropy_file <- paste0(" -y -e ",paste0(entropy[1:3],collapse = ",")," -m ",entropy[4])
          }

          if (run_entropy) {
            message(paste("ODIN will denoise",file))
            # message(paste0(dnoise," --fasta_input ",file,"_HELA_nonchimeras.fasta ",
            #                "--fasta_output ",file,"_ODIN ",
            #                "-a ",alpha," -c ",cores,entropy," -r ",min_reads_ESV))
            system(paste0(dnoise," --fasta_input ",file,"_HELA_nonchimeras.fasta ",
                          "--fasta_output ",file,"_ODIN ",
                          "-a ",alpha," -c ",cores,entropy_file," -r ",min_reads_ESV),intern = T, wait = T)
          } else  {
            message(paste("ODIN will denoise",file))
            system(paste0(dnoise," --fasta_input ",file,"_HELA_nonchimeras.fasta ",
                          "--fasta_output ",file,"_ODIN ",
                          "-a ",alpha," -c ",cores," -r ",min_reads_ESV),intern = T, wait = T)
          }
        }
    }
  }

  #####
  # 2: D,DS,S,S -> cat all fasta and perform obi uniq, annotate
  #####

  # input
  # sample_list + _ODIN_Adcorr_denoised_ratio_d.fasta | _ODIN_denoised_ratio_d.fasta | _HELA_nonchimeras.fasta
  # final outputs
  # DMS as lib"_ODIN/"

  X <- NULL
  for (file in sample_list) {
    if (algorithm=="dnoise_swarm" | algorithm=="dnoise") {
      if (run_entropy) {
        input_file <- paste0(file,"_ODIN_Adcorr_denoised_ratio_d.fasta")
      } else {
        input_file <- paste0(file,"_ODIN_denoised_ratio_d.fasta")
      }
    } else {
      input_file <- paste0(file,"_HELA_nonchimeras.fasta")
    }
    X <- c(X,paste0(
      "sed -i 's/;size/; COUNT/g' ",input_file," ; ",
      "obi import --fasta-input ",input_file, " ", file,"_ODIN/sample ; ",
      "obi annotate -S sample:\"", gsub("^[a-zA-Z0-9]{4}_", "", file), "\" ", file,"_ODIN/sample  ", file,"_ODIN/sample_name "))

  }
  clust <- makeCluster(cores)
  clusterExport(clust, list("X","old_path","obipath"),envir = environment())
  clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))})
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  stopCluster(clust)

  for (i in 1:length(sample_list)) {
    file <- sample_list[i]
    if (i==1) {
      system(paste0("obi cat -c ",
                    file,"_ODIN/sample_name",
                    " ", lib,"_ODIN/version",i),
             intern = T, wait = T)
    } else if (i==length(sample_list)) {
      system(paste0("obi cat -c ",
                    file,"_ODIN/sample_name",
                    " -c ", lib,"_ODIN/version",c(i-1),
                    " ", lib,"_ODIN/samples ; ",
                    "obi rm ", lib,"_ODIN/version",c(i-1)),
             intern = T, wait = T)
    } else {
      system(paste0("obi cat -c ",
                    file,"_ODIN/sample_name",
                    " -c ", lib,"_ODIN/version",c(i-1),
                    " ", lib,"_ODIN/version",i," ; ",
                    "obi rm ", lib,"_ODIN/version",c(i-1)),
             intern = T, wait = T)
    }
  }
  system(paste0("obi uniq --merge 'sample' ",
                lib,"_ODIN/samples ",
                lib,"_ODIN/samples_uniq"),
         intern = T, wait = T)
  if (remove_singletons) {
    system(paste0("obi grep -p \"sequence[\'COUNT\'] > 1\" ",
                  lib,"_ODIN/samples_uniq ",
                  lib,"_ODIN/seq_nosing"),
           intern = T, wait = T)
    system(paste0("obi annotate --seq-rank ",
                  lib,"_ODIN/seq_nosing ",
                  lib,"_ODIN/seq_rank"),
           intern = T, wait = T)
  } else {
    system(paste0("obi annotate --seq-rank ",
                  lib,"_ODIN/samples_uniq ",
                  lib,"_ODIN/seq_rank"),
           intern = T, wait = T)

  }
  system(paste0("obi annotate --set-identifier ",
                "\'\"\'",lib,"\'_%09d\" % sequence[\"seq_rank\"]\' ",
                lib,"_ODIN/seq_rank ",
                lib,"_ODIN/seq_id"),
         intern = T, wait = T)

  output <- system(paste0("obi ls ",lib,"_ODIN | grep 'Line count'"),intern = T,wait = T)
  values <- as.numeric(gsub(".*count: ","",output))
  version <- gsub(".*# ","",gsub(": Date.*","",output))
  after_2_ODIN <-data.frame(algorithm=algorithm,
                                         version=version,
                                         num_seqs=values)



  #####
  # 3: D,DS,S,S -> export csv file of sequences (if denoised, they are ESV)
  #####

  # input
  # sample_list + _ODIN_Adcorr_denoised_ratio_d.fasta | _ODIN_denoised_ratio_d.fasta | _HELA_nonchimeras.fasta
  # final outputs
  if (algorithm == "dnoise"){
    filetab <- paste0(lib, "_ODIN_ESV.csv")
  } else{
    filetab <- paste0(lib, "_ODIN_seqs.csv")
  }

  system(paste0("obi export --tab-output --sep ','  -o ",
                filetab, " ",
                lib, "_ODIN/seq_id"),
         intern = T, wait = T)


  #####
  # 4a: DS,SD,S -> do SWARM and create csv files of MOTUs and ESV or seqs. Also create fasta files from MOTU csv
  # 4ab: SD -> run DnoisE over the csv files
  # 4b: D -> divide fasta file into different parts
  #####

  # input
  # filetab (4a) | lib"_ODIN/seq_id" (4b)
  # final outputs
  outfile <- paste0(lib,"_ODIN")# "_part_",sprintf("%02d",part),".fasta"
  if (algorithm=="swarm_dnoise" | algorithm=="swarm" | algorithm=="dnoise_swarm"){
    outfile_MOTU <- paste0(outfile,"_counts.tsv")
    outfile_ESV <- paste0(outfile,"_ESV.tsv")
    outfile_preDnoisE <- paste0(outfile,"_seqs_to_dnoise.tsv") # + _Adcorr_denoised_ratio_d.csv | _Adcorr_denoised_ratio_d.csv
    # intermediate input/ouput
    fileswarm <- paste0(outfile,"_SWARM_output")
  }



  if (algorithm=="swarm_dnoise" | algorithm=="swarm" | algorithm=="dnoise_swarm"){ # 4a
    message("ODIN will cluster sequences into MOTUs with SWARM.")
    system(paste0("obi export --fasta-output --only-keys \"COUNT\" ",lib,"_ODIN/seq_id > ",outfile,".fasta ; ",
                  "sed -i 's/COUNT/size/g' ",outfile,".fasta ; ",
                  "sed -i 's/;//g' ",outfile,".fasta ; ",
                  "sed -E -i 's/(size=[0-9]*).*/\\1;/g' ",outfile,".fasta ; ",
                  "sed -i 's/ /;/g' ",outfile,".fasta "))
    system(paste0(swarm," -d ",d," -z -t ",cores," -o ",outfile,"_SWARM_output -s ",outfile,"_SWARM",d,"nc_stats -w ",outfile,"_SWARM_seeds.fasta ",outfile,".fasta"),intern=T,wait=T)

    get_swarm_size <- function(cadena="="){
      # it gets the position of the "=" and gets the items after it
      return(as.numeric(gsub(";","",substr(cadena,gregexpr("=",cadena)[[1]][[1]]+1,nchar(cadena)))))
    }

    message("ODIN will recount abundances for every MOTU after Swarm.")

    # Read cluster list database
    message("1. ODIN is reading SWARM results...")
    swarm_db <- readLines(fileswarm)
    total_swarms <- length(swarm_db)
    message("2. ODIN has read ", total_swarms," total MOTUs.")

    message(paste("3. ODIN will now perform a remove of MOTUs with sequences of less than",
                  min_reads_MOTU,"reads."))
    # Calculate reads in each cluster and reduce the dataset if min_reads_MOTU>0
    clusters <- strsplit(swarm_db,"; ")
    if (min_reads_MOTU>0) {
      # first reduction of a proportion of MOTUs
      if (min_reads_MOTU>9) i <- 9 else i <- min_reads_MOTU-1
      clusters <- clusters[!(grepl(paste0("size=[0-",i,"]"),clusters) & lengths(clusters)==1)]
      cluster_reads <- NULL
      # second reduction
      cluster_reads <- mclapply(clusters,function(x) sum(as.numeric(lapply(X=(x),FUN=get_swarm_size))), mc.cores = cores)
      # for (i in 1:length(clusters)) cluster_reads[i] <- sum(as.numeric(lapply(X=(clusters[[i]]),FUN=get_swarm_size)))
      clusters <- clusters[cluster_reads>=min_reads_MOTU]
      total_swarms_reduced <- length(clusters)
      after_4a_ODIN  <- data.frame(version = c("min_reads_MOTU","num_clusters","num_clusters_reduced","cluster_reads","cluster_reads_reduced"),
                                   value = c(min_reads_MOTU,total_swarms,total_swarms_reduced,sum(unlist(cluster_reads)),sum(unlist(cluster_reads[cluster_reads>=min_reads_MOTU]))))
    } else{
      total_swarms_reduced <- total_swarms
      after_4a_ODIN  <- data.frame(version = c("min_reads_MOTU","num_clusters","num_clusters_reduced","cluster_reads","cluster_reads_reduced"),
                                   value = c(min_reads_MOTU,total_swarms,total_swarms_reduced,NA,NA))
    }

    message("4. ODIN is keeping only information of the sequences that form each cluster.")
    clusters <- mclapply(clusters,function(x){sub(";.*","",x)}, mc.cores = cores)
    names(clusters) <- mclapply(clusters, function(x) x[[1]], mc.cores = cores)

    message("5. ODIN kept only ", total_swarms_reduced," MOTUs of size greater than or equal to ",min_reads_MOTU," reads.")
    motu_seqs_names <- stack(clusters) %>% rename(ID = values, MOTU = ind)
    # motu_seqs_names <- unlist(clusters, use.names=F)

    # # Generate a file with the list of ids of non-singleton clusters
    # motulist <- file(paste0(lib,"_non_singleton_motu_list.txt"),"wt")
    # writeLines(id,motulist)
    # message("ODIN has created the file ",paste0(lib,"_non_singleton_motu_list.txt")," with the list of identifiers of non-singleton MOTUs.")

    # Read counts database and keep only the needed clusters
    message("6. ODIN is reading the abundance database. This could take Him a while, since He has just one eye left, after all.")
    db <- read.table(filetab,sep=",",head=T)
    numseqs <- nrow(db)
    db <- merge(motu_seqs_names,db,by="ID")
    # db <- db[db$ID %in% motu_seqs_names,]
    numseqs_reduced <- nrow(db)
    samples <- sum(grepl("sample",names(db)))
    message("7. ODIN finished reading the Database, which includes ", numseqs," total unique sequences and ",samples," samples.\n",
            "ODIN kept only ", numseqs_reduced," sequences for calculations.")

    message("8. ODIN will now calculate the number of reads in every sample for each MOTU.")
    db.total <- split(db[,grepl("sample",names(db))], db$MOTU)
    db.total <- mclapply(db.total,function(x)as.data.frame(t(as.matrix(c(COUNT=sum(x),colSums(x),CLUST_WEIGHT=dim(x)[1])))), mc.cores = cores)
    db.total <- do.call(rbind,db.total)
    db.total <- cbind(data.frame(ID=rownames(db.total)), db.total)
    db.total <- merge(db.total,db[,grepl("ID|NUC_SEQ",names(db))],by = "ID")

    # order the columns
    col_order <- c("ID", "COUNT", "MOTU", names(db)[grepl("sample",names(db))], "NUC_SEQ" )
    db <- db[,col_order]
    s_opt <- min(grep("sample",names(db)))
    z_opt <- max(grep("sample",names(db)))

    # print datasets
    write.table(db.total,outfile_MOTU,sep="\t",quote=F,row.names=F)
    write.table(db,ifelse(algorithm=="swarm_dnoise",outfile_preDnoisE,outfile_ESV),sep="\t",quote=F,row.names=F)

    # divide dataset into different files for THOR
    db.total <- db.total[,c("ID","NUC_SEQ")]
    db.total <- paste(paste0(">",db.total$ID),db.total$NUC_SEQ,sep="\n")
    db.total <- split(db.total, factor(sort(1:length(db.total)%%cores)))
    for (part in 1:length(db.total)) {
      writeLines(paste0(db.total[[part]],collapse = "\n"),paste0(outfile,"_part_",sprintf("%02d",part),".fasta"))
    }

    message("File ", outfile, " written")

    if (algorithm=="swarm_dnoise") { # 4ab
      message("ODIN will generate now a list of ESVs for every non-singleton MOTU, using DnoisE.")
      if(entropy[1]=="auto_dataset" | entropy[1]=="auto_sample") {
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
        system(paste0(dnoise," --csv_input ",outfile_preDnoisE," ",
                      "--csv_output ",outfile_ESV," ",
                      "-a ",alpha," -c ",cores," -n 'COUNT' -p 1 -q 'NUC_SEQ' ",
                      "-s ", s_opt," -z ", z_opt," ",
                      entropy," -w 'MOTU' -r ",min_reads_ESV, " ; ",
                      "mv ",outfile_ESV,"_Adcorr_denoised_ratio_d.csv ",outfile_ESV),intern = T, wait = T)
        remaining_ESV <- as.numeric(system(paste0("wc -l ",outfile_ESV," |  cut -f1 -d ' ' "),intern = T, wait = T))-1
      } else  {
        system(paste0(dnoise," --csv_input ",outfile_preDnoisE," ",
                      "--csv_output ",outfile_ESV," ",
                      "-a ",alpha," -c ",cores," -n 'COUNT' -p 1 -q 'NUC_SEQ' ",
                      "-s ", s_opt," -z ", z_opt," ",
                      "-w 'MOTU' -r ",min_reads_ESV, " ; ",
                      "mv ",outfile_ESV,"_denoised_ratio_d.csv ",outfile_ESV),intern = T, wait = T)
        remaining_ESV <- as.numeric(system(paste0("wc -l ",outfile_ESV," |  cut -f1 -d ' ' "),intern = T, wait = T))-1
      }
      after_4a_ODIN <- rbind(after_4a_ODIN,data.frame(version = c("remaining_ESV"),
                                                      value = c(remaining_ESV)))
      message("")
    }

  } else { # 4b

    if (cores == 1) {
      system(paste0("obi export --fasta-output --only-keys \"COUNT\" ",lib,"_ODIN/seq_id > ",paste0(outfile,"_part_01.fasta")), intern = T, wait = T)
    } else {
      system(paste0("obi export --fasta-output --only-keys \"COUNT\" ",lib,"_ODIN/seq_id > ",outfile,".fasta "), intern = T, wait = T)
      fasta_file <- readLines(paste0(outfile,".fasta"))

      seqs <- grep(">",file_read)
      seq_rank <- sort(seqs%%cores)

      for (i in 1:cores) {
        if (i == 1) {
          start_line <- 1
          end_line <- seqs[grep(unique(seq_rank)[i+1],seq_rank)[1]]-1
        } else if (i == cores) {
          start_line <- seqs[grep(unique(seq_rank)[i],seq_rank)[1]]
          end_line <- length(fasta_file)
        } else {
          start_line <- seqs[grep(unique(seq_rank)[i],seq_rank)[1]]
          end_line <- seqs[grep(unique(seq_rank)[i+1],seq_rank)[1]]-1
        }
        writeLines(paste0(fasta_file[start_line:end_line],
                          collapse = "\n"),
                   paste0(outfile,"_part_",sprintf("%02d",i),".fasta"))
      }
    }
  }
  save(file = "summary_ODIN.RData",list = c(c("alpha","min_reads_MOTU","min_reads_ESV","algorithm","run_entropy","entropy","after_2_ODIN"),
                                     c("before_1_ODIN")[exists("before_1_ODIN")],c("after_4a_ODIN")[exists("after_4a_ODIN")]))
  if (remove_DMS) {
    system(paste0("rm -r *ODIN.obidms "),intern=T,wait=T)
  }
  message("ODIN is done.")
}

