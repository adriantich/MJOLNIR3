#' LOKI: LULU Overseeing with Kinship Identification
#' 
#' LOKI is a convenient wrapper of LULU for the MJOLNIR3 metabarcoding pipeline.
#' 
#' @details 
#' LOKI starts from the combined dataset of abundances and taxonomy from the 
#' previous step (FRIGGA): <EXPX>_FRIGGA.tsv.
#' A match list of representative MOTU sequences is created using VSEARCH and saved as a txt file.
#' Then Units that are potential errors based on co-occurrence patterns
#' are labelled and removed using LULU.
#' The output is called <EXPX>_LOKI_Curated.tsv
#' 
#' @param lib Character string. Acronym for the experiment. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. However they can be the same.
#' 
#' @param min_id Numeric. Equivalent to the --id option from vsearch --usearch_global
#' 
#' From vsearch manual: Reject the sequence match if the pairwise identity is 
#' lower than min_id (value ranging from 0.0 to 1.0 included). The search process 
#' sorts target sequences by decreasing number of k-mers they hav e in common 
#' with the query sequence, using that information as a proxy for sequence 
#' similarity. That efficient pre-filtering also prevents pairwise alignments 
#' with weakly matching targets, as there needs to be at least 6 shared kmers 
#' to start the pairwise alignment, and at least one out of every 16 k-mers 
#' from the query needs to match the target. Consequently, using values lower 
#' than --id 0.5 is not likely to capture more weakly matching targets. The 
#' pairwise identity is by default defined as the number of 
#' (matching columns) / (alignment length - terminal gaps)
#' 
#' @param vsearchpath Character string specifying the PATH to vsearch.
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
#' 
#' # set the directory where the database is stored
#' tax_dir <- "~/taxo_NCBI/"
#' tax_dms_name <- "DUFA_COI"
#' 
#' # Run THOR
#' mjolnir5_THOR(lib, cores, tax_dir=tax_dir, tax_dms_name=tax_dms_name, run_ecotag=T)
#' 
#' # Run FRIGGA
#' mjolnir6_FRIGGA(lib)
#' 
#' # Run LOKI
#' mjolnir7_LOKI(lib, min_id=.84)

mjolnir7_LOKI <- function(lib,min_id = .84, vsearchpath = NULL){

  # lib is the name of the library to be processed. Usually a four-character uppercase name.
  # Input file name must be in the format: LIBR.All_MOTUs.tsv. Then lib must be = "LIBR"
  # min_id is the minimum identity between two sequences to be kept in the match_list output. Default: 0.84

  if (is.null(vsearchpath)) vsearchpath <- "~/vsearch-2.22.1/bin/"
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(old_path, path.expand(vsearchpath), sep = ":"))

  message("LOKI will produce a pairwise match list for LULU.")
  system(paste0("cat ",lib,"_ODIN_part_??.fasta > ",lib,"_LOKI.fasta"))
  system(paste0("vsearch --usearch_global ",lib,"_LOKI.fasta --db ",lib,"_LOKI.fasta --self --id ",min_id," --iddef 1 --userout ",lib,"_LOKI_match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"),intern=T,wait=T)
  message("LOKI will now remove the pseudogenes with LULU.")

  suppressPackageStartupMessages(library("lulu"))
  suppressPackageStartupMessages(library("dplyr"))

  #Load the dataset
  db <- read.table(paste0(lib,"_FRIGGA.tsv"),sep="\t",head=T,stringsAsFactors = F)
    # Select sample abundance columns
  sample_cols <- grep("sample",names(db))
  otutable_name <-db[,sample_cols]
  rownames(otutable_name) <- db$id

  #Load the matchlist
  matchlist_name <- read.csv(paste0(lib,"_LOKI_match_list.txt"),sep="\t",head=F,stringsAsFactors = F)

  #Run LULU
  curated_result <- lulu(otutable_name, matchlist_name)

  #Get discarded table:
  discarded_db <- db[db$id %in% curated_result$discarded_otus,]
  otus_discarded <- curated_result$discarded_otus
  num_discarded <- nrow(discarded_db)
  write.table(discarded_db,paste0(lib,"_LOKI_Discarded.tsv"),row.names = F,sep="\t",quote = F)

  #Get curated table:
  curated_db <- db[db$id %in% curated_result$curated_otus,]
  curated_db <- curated_db[order(curated_db$id),]
  curated_db[,sample_cols] <- curated_result$curated_table
  curated_db$total_reads <- rowSums(curated_result$curated_table)
  write.table(curated_db,paste0(lib,"_LOKI_Curated.tsv"),row.names = F,sep="\t",quote = F)

  #Get fate of deleted taxa
  deleted_otu_fate <- (curated_result$otu_map[curated_result$otu_map$curated=="merged",])
  deleted_otu_fate$original <- ""
  deleted_otu_fate$id_removed <- rownames(deleted_otu_fate)
  for (i in 1:nrow(deleted_otu_fate)){
    deleted_otu_fate$original[i] <- db$SCIENTIFIC_NAME[db$id==rownames(deleted_otu_fate)[i]]
  }
  parents_info <- db[db$id%in%deleted_otu_fate$parent_id,c("id","SCIENTIFIC_NAME","superkingdom_name","kingdom_name","phylum_name","class_name","order_name","family_name")]
  names(parents_info)[1:2] <- c("parent_id","parent_taxo")
  deleted_otu_fate <- left_join(deleted_otu_fate,parents_info,by="parent_id")
  write.table(deleted_otu_fate,paste0(lib,"_LOKI_Deleted_fate.tsv"),row.names = F,sep="\t",quote = F)

  save(file = "summary_LOKI.RData",list = c("num_discarded","otus_discarded"))

  message("LOKI is done. He kept ",nrow(curated_db)," MOTUs in the curated database, which He stored in file: ",paste0(lib,"_LOKI_Curated.tsv"))
  message("LOKI discarded ",nrow(discarded_db)," MOTUs, which He stored in file: ",paste0(lib,"_LOKI_Discarded.tsv"))
  message("LOKI stored the fate of the discarded MOTUs in file: ",paste0(lib,"_LOKI_Discarded.tsv"))
}
