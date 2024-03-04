#' RAGNAROC: Replace AGnomens with Names And Recover Original Codification
#' 
#' Final step of the MJOLNIR3 pipeline to apply the last filters to the abundance data
#' 
#' @details 
#' RAGNAROC consists on different contamination removals and filtering steps 
#' as follows:
#' 
#' Removal of Bacteria: this removed the Units tagged as "Prokaryota" or "root" in the <EXPX>_LOKI_Cutated.tsv
#' 
#' Removal of contaminations: this step removes the taxa specified in the "contaminaion_file"
#' 
#' NUMT removal: this step is design for the Leray-XT COI marker. It deletes all
#' sequences that do not have a 313 (plus/minus multiple of 3 equivalent to a codon)
#' bp length. Then removes sequences with stop codons and those metazoan 
#' sequences that do not translate for 5 specific conservative aminoacids.
#' 
#' @param experiment Character string. Acronym for the experiment. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. However they can be the same.
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
#' @param output_file Character string specifying the outputfile name
#' 
#' @param output_file_ESV Character string specifying the outputfile name 
#' for ESVs abundances if required.
#' 
#' @param min_reads Number of the minimum number of reads allowed for each MOTU/ESV
#' or ESV within MOTU.
#'  
#' @param remove_bacteria Logical. If TRUE it will apply the bacteria removal 
#' filtering (see Details).
#' 
#' @param remove_contamination Logical. If TRUE it will apply the contamination removal 
#' filtering (see Details).
#' 
#' @param contamination_file Character string specifying the name of the contamination file. (see Details)
#' 
#' @param ESV_within_MOTU Logical. If TRUE this will take into account the ESV 
#' that were clustered into MOTUs in ODIN if algorithm was set to "DnoisE_SWARM" 
#' or "SWARM_DnoisE" and apply all filters to both data.
#' 
#' @param remove_numts Logical whether to apply the NUMT filter (TRUE) or not (FALSE)
#' 
#' @param cores Numeric. Number of threads for parallel processing during NUMT 
#' removal.
#' 
#' @export 
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
#' 
#' # set the directory where the database is stored
#' tax_db <- "~/taxo_NCBI/DUFA_COI"
#' 
#' # Run THOR
#' mjolnir5_THOR(experiment = experiment, cores = cores, 
#'               tax_db = tax_db, run_ecotag = T)
#' 
#' # Run FRIGGA
#' mjolnir6_FRIGGA(experiment = experiment)
#' 
#' # Run LOKI
#' mjolnir7_LOKI(experiment = experiment, min_id = .84)
#' 
#' # Run RAGNAROC
#' mjolnir8_RAGNAROC(experiment = experiment, remove_bacteria = T, ESV_within_MOTU = T,
#'                   remove_numts = F, cores = 1)

mjolnir8_RAGNAROC <- function(experiment = NULL, metadata_table = "",
                              output_file = "", output_file_ESV = "",
                              min_reads = 0, remove_bacteria = T,
                              remove_contamination = F,
                              contamination_file = "contaminants.txt",
                              ESV_within_MOTU = T,
                              remove_numts = F, cores = 1, ...){

  
  if (exists("lib") && is.null(experiment)) {
    # Use lib as experiment
    experiment <- lib
    # Print deprecation warning
    warning("The 'lib' argument is deprecated. Please use 'experiment' instead.")
  }
  
  message("RAGNAROC is coming. Original sample names will be recovered.")
  suppressPackageStartupMessages(library(tidyr))
  if (output_file == "") {
    output_file <- paste0(experiment,"_RAGNAROC_final_dataset.tsv")
  }
  if (output_file_ESV == "") {
    output_file_ESV <- paste0(experiment,"_RAGNAROC_final_dataset_ESV.tsv")
  }

  if (file.exists("summary_FREYJA.RData")) {
    load("summary_FREYJA.RData")
    rownames(variables_FREYJA) <- variables_FREYJA$variable
    FREYJA <- T
  } else {
    FREYJA <- F
  }
  if (file.exists("summary_HELA.RData")) {
    load("summary_HELA.RData")
    HELA <- T
  } else {
    HELA <- F
  }
  if (file.exists("summary_ODIN.RData")) {
    load("summary_ODIN.RData")
    ODIN <- T
  } else {
    ODIN <- F
  }
  if (file.exists("summary_LOKI.RData")) {
    load("summary_LOKI.RData")
    LOKI <- T
  } else {
    LOKI <- F
  }

  # Load the data set
  # two data sets can be used:
  # MOTU/ESV + taxa info with LULU/without LULU
  # ESV within MOTU only available for ODIN_algorithm DnoisE_SWARM/SWARM_DnoisE

  if (LOKI || file.exists(paste0(experiment,"_LOKI_Curated.tsv"))) {
    db <- read.csv(paste0(experiment,"_LOKI_Curated.tsv"), sep="\t",head=T,stringsAsFactors = F)
  } else {
    db <- read.csv(paste0(experiment,"_FRIGGA.tsv"), sep="\t",head=T,stringsAsFactors = F)
  }
  if (ESV_within_MOTU) {
    ESV_data_initial <- read.csv(paste0(experiment,"_ODIN_ESV.tsv"), sep="\t",head=T,stringsAsFactors = F)
  }

  # Remove bacteria
  if (remove_bacteria) {
    message("RAGNAROC is removing bacterial MOTUs now.")
    bacteria_removed <- sum(c(db$superkingdom_name == "Prokaryota" | db$SCIENTIFIC_NAME == "root"),na.rm = T)
    db <- db[(!grepl('Prokaryota',db$superkingdom_name) & !grepl('Prokaryota',db$superkingdom_name)),]
  }

  # Remove contamination
  if (remove_contamination){
    message("RAGNAROC is removing contaminant MOTUs now.")
    contamination <- readLines(contamination_file)
    db <- db[!((db$SCIENTIFIC_NAME %in% contamination) |
                         (db$phylum_name %in% contamination) |
                         (db$class_name %in% contamination) |
                         (db$order_name %in% contamination) |
                         (db$family_name %in% contamination) |
                         (db$genus_name %in% contamination)) ,]
  }

  # Load the metadata_table
  if (metadata_table=="") metadata_table <- paste0(experiment,"_metadata.tsv")
  sample_db <- read.table(metadata_table,sep="\t",head=T,stringsAsFactors = F)

  message("data loaded")

  # if ESV_within_MOTU remove the MOTUs that LOKI has removed and bacteria filter
  if (ESV_within_MOTU) {
    ESV_data_initial <- ESV_data_initial[ESV_data_initial$MOTU%in%db$id,]
  }

  # Select sample abundance columns
  sample_cols <- grep("sample",names(db))
  # initial_no_sample_cols <- length(sample_cols)
  sample_names <- names(db[sample_cols])

  # Change agnomens by original names 
  # Those names that are not present in the metadata will be tagged as
  # EMPTY and then removed
  new_sample_names <- sample_db$original_samples[match(sample_names,sample_db$mjolnir_agnomens)]
  new_sample_names[is.na(new_sample_names)] <- gsub("^","EMPTY",as.character(c(1:sum(is.na(new_sample_names)))))
  colnames(db)[sample_cols] <- new_sample_names

  db <- db[,!grepl("EMPTY",colnames(db))]

  sample_names <- sample_names[!grepl("EMPTY",sample_names)]
  sample_cols <- match(sample_names,colnames(db))
  sample_cols <- sample_cols[!is.na(sample_cols)]
  
  # same for ESVs
  if (ESV_within_MOTU) {
    # Select sample abundance columns
    sample_cols_ESV <- grep("sample",names(ESV_data_initial))
    sample_names_ESV <- gsub("MERGED_sample.","",names(ESV_data_initial[sample_cols_ESV]))

    # Change agnomens by original names
    new_sample_names_ESV <- sample_db$original_samples[match(sample_names_ESV,sample_db$mjolnir_agnomens)]
    new_sample_names_ESV[is.na(new_sample_names_ESV)] <- gsub("^","EMPTY",as.character(c(1:sum(is.na(new_sample_names_ESV)))))
    names(ESV_data_initial)[sample_cols_ESV] <- new_sample_names_ESV
  }
  
  if (!ESV_within_MOTU) {
    rownames(db) <- db$id

    db$COUNT <- rowSums(db[,sample_cols])

    message("RAGNAROC is removing MOTUs with less than ",min_reads," total reads.")
    db <- db[db$COUNT >= min_reads,]
  } else {
    rownames(ESV_data_initial) <- ESV_data_initial$ID

    ESV_data_initial$COUNT <- rowSums(ESV_data_initial[,sample_cols_ESV])

    message("RAGNAROC is removing MOTUs with less than ",min_reads," total reads.")
    ESV_data_initial <- ESV_data_initial[ESV_data_initial$COUNT >= min_reads,]

    # remove numts
    if (remove_numts) {
      message("numts will be removed")
      # no_ESV_before_numts <- dim(ESV_data_initial)[1]
      lengths <- nchar(as.vector(ESV_data_initial$NUC_SEQ))
      ESV_data_initial <- ESV_data_initial[(lengths-313)%%3 == 0,]
      lengths <- nchar(as.vector(ESV_data_initial$NUC_SEQ))

      # no_numts_data <- c()
      # numts_seqs <- c()
  
      number_of_motus <- length(unique(ESV_data_initial$MOTU))
      motu_taxa <- data.frame("id" = db$id, "Metazoa" = c(db$kingdom_name == "Metazoa" & !is.na(db$kingdom_name)))
      numts_ESV <- parallel::mclapply(1:number_of_motus,function(i,ESV_data_initial,motu_taxa){
        motu <- unique(ESV_data_initial$MOTU)[i]
        datas <- ESV_data_initial[ESV_data_initial$MOTU==motu,]
        is_metazoa <- motu_taxa$Metazoa[motu_taxa$id==as.character(motu)]
        datas_length <- nchar(as.vector(datas$NUC_SEQ))
        newlist <- numts(datas, is_metazoa = is_metazoa, motu = motu, datas_length = datas_length)
        return(newlist)
      },ESV_data_initial=ESV_data_initial,motu_taxa=motu_taxa,mc.cores = cores)
      numts_ESV <- do.call("rbind",numts_ESV)
      ESV_data_initial <- ESV_data_initial[!ESV_data_initial$ID %in% numts_ESV$id,]
      message("numts removed")
    }
    db <- db[db$id %in% unique(ESV_data_initial$MOTU),]
    # compute new abundances of MOTUs from ESV
    motu_abund <- parallel::mclapply(db$id, FUN = function(x,ESV_data_initial,sample_cols_ESV){
      data_motu <- ESV_data_initial[as.character(ESV_data_initial$MOTU) == x,sample_cols_ESV]
      if (dim(data_motu)[1]>1) {
        return(colSums(data_motu))
      } else {
        return(data_motu)
      }
    },ESV_data_initial=ESV_data_initial,sample_cols_ESV=sample_cols_ESV,mc.cores = cores)
    db[,sample_cols] <- do.call("rbind",motu_abund)
  }

  # Write final table
  write.table(db,output_file,row.names = F,sep="\t",quote = F)
  if (ESV_within_MOTU){
    write.table(ESV_data_initial,output_file_ESV,row.names = F,sep="\t",quote = F)
  }
  message("After RAGNAROC, MJOLNIR is done. File: ",output_file, " written with ",nrow(db), " MOTUs and ",sum(db$total_reads)," total reads.")
  
  #####
  # RAGNAROC REPORT
  #####
  # RAGNAROC_report <- paste("Dear friend,\n",
  #            "you have succesfully arrived at the end of RAGNAROC. You've meet gods and took their help to twist the data to your will.\n",
  #            "After RAGNAROC the rest is up to you. Don't lose the faith in your experiment, the end is near but new paths will open below your feet.\n",
  #            "Please don't forget to cite and thank the two dwarfs Cindri and Brok, AKA Owen and Adria, for the forge of my self.\n",
  #            "MJOLNIR.\n",
  #            "P.S.: See below for a small summary of your journey.\n")
  # if (FREYJA) {
  #   if (as.logical(variables_FREYJA["demultiplexed",2])) {
  #     RAGNAROC_report <- paste(RAGNAROC_report, "You started FREYJA with your samples allready demultiplexed and with the following sequences for each file \n")
  #     do.call("rbind",before_FREYJA)
  #   } else{
  #     RAGNAROC_report <- paste(RAGNAROC_report, "You started FREYJA with the following sequences for each file \n")
  #     RAGNAROC_report <- paste0(RAGNAROC_report,paste(apply(do.call("rbind",lapply(before_FREYJA,function(x)do.call("rbind",x))),1,paste, collapse =" : "), collapse = "\n"),"\n")
  #   }
  #   RAGNAROC_report <-  paste(RAGNAROC_report, "You used",variables_FREYJA["cores",2]," cores to aling your sequences. You choosed those sequences with a quality score of more than",variables_FREYJA["score_obialign",2],".\n",
  #              "You assign each sequence to a sample name and removed the primer's sequences.\n",
  #              "Finally in FREYJA you just kept those sequences with A, G, T or C's and with a sequence length between",variables_FREYJA["Lmin",2],"and",variables_FREYJA["Lmax",2],"bp.\n",
  #              "The resulting files had the following stats:\n")
  #   RAGNAROC_report <- paste0(RAGNAROC_report,"sample\tfiltered sequences\tuniq sequences\n",
  #                             paste(apply(as.data.frame(pivot_wider(do.call("rbind",after_FREYJA),names_from = "version",values_from = "num_seqs")),1,paste0, collapse ="\t"), collapse = "\n"),
  #                             "\n")
  # } else {
  #   RAGNAROC_report <-  paste(RAGNAROC_report, "Sorry but I couldn't find a summary of your FREYJA process \n")
  # }
  # if (HELA) {
  #   RAGNAROC_report <-  paste(RAGNAROC_report, "HELA removed the chimeras with the uchime_denovo algorithm and kept for each sample the following number of non-chimeras:\n")
  #   RAGNAROC_report <-  paste0(RAGNAROC_report,"sample\tsequences\n",
  #                              paste(apply(after_HELA,1,paste0, collapse ="\t"), collapse = "\n"),
  #                              "\n")
  #   
  # } else {
  #   RAGNAROC_report <-  paste(RAGNAROC_report, "Sorry but I couldn't find a summary of your HELA process \n")
  # }
  # if (ODIN) {
  #   RAGNAROC_report <-  paste(RAGNAROC_report, "ODIN was used to obtain meaningful units. In your case you chose the ",algorithm," algorithm.\n")
  #   if (algorithm=="dnoise_swarm" || algorithm=="dnoise") {
  #     RAGNAROC_report <-  paste(RAGNAROC_report, "ODIN used DnoisE to obtain the ESV's of your samples running within them with the following options:\n")
  #     if (!is.logical(entropy)) {
  #       RAGNAROC_report <-  paste(RAGNAROC_report, "Entropy correction with sequences delimited to a multiple of 313bp, alpha",alpha,"and minimum number of reads of",min_reads_ESV,"\n")
  #     } else {
  #       RAGNAROC_report <-  paste(RAGNAROC_report, "Alpha",alpha,"and minimum number of reads of",min_reads_ESV,"\n")
  #     }
  #   }
  #   if (algorithm=="dnoise_swarm"  | algorithm=="swarm" | algorithm=="swarm_dnoise") {
  #     if (exists("after_2_ODIN")&exists("after_4a_ODIN")) {
  #       RAGNAROC_report <-  paste(RAGNAROC_report, "ODIN joined all the sequences, obtained the unique ones and applied swarm to obtain the MOTUs. Before SWARM you had",after_2_ODIN$values[after_2_ODIN$version=="seq_id"],"sequences and at the end you obtained the following stats: \n")
  #     }
  #   } else{
  #     if (exists("after_2_ODIN")) {
  #       RAGNAROC_report <-  paste(RAGNAROC_report, "The samples were then grouped and the unique sequences obtained being",after_2_ODIN$values[after_2_ODIN$version=="seq_id"],"sequences in total.\n")
  #     }
  #   }
  #   if (algorithm=="swarm_dnoise") {
  #     RAGNAROC_report <-  paste(RAGNAROC_report, "ODIN used DnoisE to obtain the ESV's of your samples running within them with the following options:\n")
  #     if (run_entropy) {
  #       RAGNAROC_report <-  paste(RAGNAROC_report, "Entropy correction (",entropy,"), alpha ",alpha," and minimum number of reads of ",min_reads_ESV,"\n")
  #     } else {
  #       RAGNAROC_report <-  paste(RAGNAROC_report, "Alpha ",alpha," and minimum number of reads of ",min_reads_ESV,"\n")
  #     }
  #   }
  # } else {
  #   RAGNAROC_report <-  paste(RAGNAROC_report, "Sorry but I couldn't find a summary of your ODIN process \n")
  # }
  # if (LOKI) {
  #   RAGNAROC_report <-  paste(RAGNAROC_report, "LOKI used LULU to search for potential pseudogenes and found ",num_discarded," OTUs that were discarded.\n")
  # } else {
  #   RAGNAROC_report <-  paste(RAGNAROC_report, "Sorry but I couldn't find a summary of your LOKI process\n")
  # }
  # RAGNAROC_report <-  paste(RAGNAROC_report, "During RAGNAROC some filters were applied.\n")
  # if (remove_bacteria) RAGNAROC_report <-  paste(RAGNAROC_report, bacteria_removed," bacteria were removed\n")
  # if (remove_contamination) RAGNAROC_report <-  paste(RAGNAROC_report, "contaminations were removed\n")
  # if (dim(neg_samples)[2]>0) RAGNAROC_report <-  paste(RAGNAROC_report, dim(data_neg_filt_deleted)[1],ifelse(ESV_within_MOTU," ESV"," MOTU")," were removed by neg/blank filter\n")
  # RAGNAROC_report <-  paste(RAGNAROC_report, "The relative abundance filter of ",min_relative," within samples had effect on", dim(relabund_changed)[2],"id's\n")
  # if (ESV_within_MOTU&remove_numts) RAGNAROC_report <-  paste(RAGNAROC_report, "The numts filter found",dim(numts_ESV)[1],"numts that were removed.\n")
  # sink("RAGNAROC_summary_report.txt")
  # cat(RAGNAROC_report)
  # sink()
}

#' numts
#' 
#' Internal function of RAGNAROC to detect numt ESV within MOTUs
#' 
#' @details 
#' numt function will:
#' 1: keep sequences of the same length as the seed
#' 2: remove missaligned sequences
#' 3: search for codon stops using all the mitochondrial genome codes with the 
#' Biostrings package and keep the one with less reads with stop codons. 
#' 4: if the MOTU is assigned to metazoans and the seed has 313 bp then check it 
#' the following aminoacid positions have the conserved ones in metazoans:
#' 20 = "H"; 23 ="G"; 32 = "N"; 81 ="D"; 95 ="G" 
#' 
#' @param data Data Frame of the sequences in the MOTU.
#' 
#' @param is_metazoa Logical. whether if the MOTU is assigned to metazoa or not
#' 
#' @param motu Character. Name of the MOTU evaluated
#' 
#' @param datas_length Numeric vector. Vector with the number of characters in 
#' each sequence


numts<-function(datas, is_metazoa=FALSE, motu, datas_length) {
  suppressPackageStartupMessages(library(Biostrings))
  suppressPackageStartupMessages(library(stringr))
  
  compare.DNA <- function(x,y){
    as.integer(x) == as.integer(y)
  }

  # compare only mitochondrial genetic code
  mitochondrial_GC <- c(2,3,4,5,7,11,12,14,15,16,17,18)
  # START
  motu_name = motu
  datas$NUC_SEQ<-as.character(datas$NUC_SEQ)

  # remove sequences with different length than the seed
  if (sum(datas$ID==motu)==0) { # if the seed has been deleted in previous steps take the first more abundant
    motu = datas$ID[which(datas$COUNT==max(datas$COUNT,na.rm = TRUE))[1]]
  }
  correct_length <- datas_length[datas$ID==motu]
  datas <- datas[datas_length==correct_length,]

  # remove misaligned sequences (more than 30 differences between a sequence and
  # the seed)
  misaligned_seqs <- c()
  motu_seq <- DNAString(datas$NUC_SEQ[datas$ID == motu])
  for (i in 1:dim(datas)[1]) {
    if(sum(!compare.DNA(motu_seq,DNAString(datas$NUC_SEQ[i])))>=30){
      misaligned_seqs <- c(misaligned_seqs,i)
    }
  }
  if (length(misaligned_seqs)>0) {
    datas <- datas[-misaligned_seqs,]
  }

  # look for the best genetic code, this is the one with less stop codons in all sequences
  # the number of codons stop is multiplied by the number of count of the sequence.
  stops<-matrix(NA,dim(datas)[1],20)
  aa_xung<-matrix(NA,dim(datas)[1],20)
  seq<-DNAStringSet(datas$NUC_SEQ)
  seq<-DNAStringSet(seq,start=2,end=nchar(datas$NUC_SEQ[1]))

  for (qq in mitochondrial_GC){
    code<-getGeneticCode(as.character(GENETIC_CODE_TABLE$id[qq]))
    trans<-translate(seq,genetic.code=code)

    # for (k in 1:dim(datas)[1]){
    # nstops <- stringr::str_count(as.character(trans[k]),fixed("*"))
    # stops[k,qq] <- nstops * datas$COUNT[k]
    # }
    nstops <- apply(data.frame(as.character(trans)), 1, function(x){stringr::str_count(x,fixed("*"))})
    stops[,qq] <- nstops * datas$COUNT

  }

  goodcodes<-which(colSums(stops)==min(colSums(stops),na.rm = T))

  # if more than one code have been chosen as good code choose the first as the best.
  # However, if the MOTU is a Metazoan and has 313bp length we check the 5 well preserved aa
  # and the code with less changes per read is the chosen. Also remove sequences
  # with changes in thees positions

  if (is_metazoa & (correct_length == 313)){

    for (qq in 1:length(goodcodes))
    {
      code<-getGeneticCode(as.character(GENETIC_CODE_TABLE$id[goodcodes[qq]]))
      trans<-translate(seq,genetic.code=code)
      for (k in 1:dim(datas)[1])
      {
        aa<-strsplit(as.character(trans[k]),split="")
        aa<-unlist(aa)
        bad_aa<-0
        if (aa[20]!="H") {bad_aa<-bad_aa+datas$COUNT[k]} # the number of errors counted are the same as the number of counts of the seq.
        if (aa[23]!="G") {bad_aa<-bad_aa+datas$COUNT[k]}
        if (aa[32]!="N") {bad_aa<-bad_aa+datas$COUNT[k]}
        if (aa[81]!="D") {bad_aa<-bad_aa+datas$COUNT[k]}
        if (aa[95]!="G") {bad_aa<-bad_aa+datas$COUNT[k]}
        aa_xung[k,goodcodes[qq]]<-bad_aa
      }
    }

    goodcodes <- goodcodes[which(colSums(aa_xung)[goodcodes]==min(colSums(aa_xung)[goodcodes],na.rm = T))]
    bestcode <- goodcodes[1]
    bestcodename <- GENETIC_CODE_TABLE$name[bestcode]
    goodcodesnames <- GENETIC_CODE_TABLE$name[goodcodes]
    flag <- stops[,bestcode]>0 | aa_xung[,bestcode]>0
  } else {
    bestcode<-which(colSums(stops)==min(colSums(stops),na.rm = T))[1]
    bestcodename<-GENETIC_CODE_TABLE$name[bestcode]
    goodcodesnames <- GENETIC_CODE_TABLE$name[goodcodes]
    flag <- stops[,bestcode]>0
  }

  # numts
  if (sum(flag)>0) {
    numts_seqs <- data.frame("motu" = motu_name, "id" = datas$ID[flag],
                             "genetic_code" = bestcodename,
                             "similar_codes" = paste(goodcodesnames, collapse = " | "))
  } else {
    numts_seqs <- c()
  }

  # # remove numts
  # datas <- datas[(flag==FALSE),]
  #
  # newlist <- list("no_numts_data" = datas, "numts_seqs" = numts_seqs)

  return(numts_seqs)
}
