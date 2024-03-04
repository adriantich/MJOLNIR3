#' THOR: Taxonomy with Higher-than-Order Ranks
#'
#' This is a wrapper of ecotag
#'
#' @details 
#' After assignment with ecotag, higher taxa at ranks higher than order are added
#' from cust. The database used can be download or build using the NJORDR package
#' (see https://github.com/adriantich/NJORDR-MJOLNIR3)
#' 
#' For vsearch assignment, the database must be in fasta format and the taxonomy
#' for the output CSV will be only in one single column.
#' 
#' @param experiment Character string. Acronym for the experiment. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. However they can be the same.
#' 
#' @param cores Numeric. Number of threads for parallel processing.
#' 
#' @param tax_db Character string specifying de PATH to the reference database. 
#' In case of using ecotag, the database must be an .obidms object. Also, when 
#' using ecotag it is important to have the following files in the same directory:
#' order.complete.csv ; family_to_order.csv ; genus_to_family.csv
#' 
#' @param run_ecotag Logical. Whether to run (TRUE, default) the ecotag taxonomic
#' assignment or not (FALSE). The latter could take place when alternative taxonomic
#' assignament software is applied but adding higher taxonomic ranks is desired.
#' 
#' @param vsearch Logical. Whether to run (TRUE) the vsearch taxonomic
#' assignment or not (FALSE, default). If vsearch has been selected, even if 
#' run_ecotag=TRUE, vsearch will be used instead of ecotag.
#' 
#' @param remove_DMS Logical. If TRUE, it will delete all obidms objects that are
#' created during the process. This can save a lot of hard disk space. The FALSE 
#' option is useful for developing and debugging.
#' 
#' @param minimum_circle Numeric. For ecotag: Minimum identity considered for 
#' the assignment circle (sequence is assigned to the LCA of all sequences 
#' within a similarity circle of the best matches; the threshold for this circle 
#' is the highest value between <CIRCLE_THRESHOLD> and the best assignment score 
#' found for the query sequence). Give value as a normalized identity, e.g. 0.95 
#' for an identity of 95%. For vsearch it is equivalent to the minimum id required
#' for a match. Default: 0.7
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

mjolnir5_THOR <- function(experiment = NULL, cores = 1,
                          tax_db = NULL,
                          run_ecotag = T, vsearch = F, remove_DMS = T, 
                          minimum_circle = 0.7, ...){
  
  if (exists("lib") && is.null(experiment)) {
    # Use lib as experiment
    experiment <- lib
    # Print deprecation warning
    warning("The 'lib' argument is deprecated. Please use 'experiment' instead.")
  }
  if ((exists("tax_dms_name") | exists("tax_dms_name")) && is.null(tax_db)) {
    # Print deprecation warning
    warning("The 'tax_dms_name' and the 'tax_dir' arguments are deprecated. Please use only 'tax_db' instead.")
    tax_db <- paste0(normalizePath(tax_dir), "/", tax_dms_name)
  } 

  tax_db_dir <- normalizePath(dirname(tax_db))
  tax_db_name <- basename(tax_db)
  divide_fasta <- function(cores = cores, experiment = experiment) {
    if (cores > 1) {
      message("THOR will use ", cores, " cores for parallel processing.")
      fasta_file <- readLines(paste0(experiment,"_ODIN.fasta"))
      
      seqs <- grep(">",fasta_file)
      seq_rank <- sort(seqs %% cores)
      
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
                   paste0(experiment, "_THOR_pretaxa_part_",sprintf("%02d",i),".fasta"))
      }
    } else {
      message("THOR will use 1 core for processing.")
    }
  }
  


  if (!vsearch) {
    if (run_ecotag) {
      message("THOR will assign the taxonomy to the order level with ecotag.")
      message("THOR will first clear the battle field.")
      message("Any directory or file containing the word THOR will be removed.")
      system("rm -r *THOR*",intern = T, wait = T)
      
      divide_fasta(cores = cores, experiment = experiment)
      
      # check whether the taxonomic folder exists
      tax_dir <- list.dirs(path = tax_db_dir, recursive = FALSE, full.names = TRUE)
      tax_dir <- tax_dir[grep(tax_db_name,tax_dir)]
      if (length(tax_dir) == 0) {
        message(paste0("The taxonomic db '", tax_db_name,"' does not exist in '", tax_db_dir,"'. Please check the path."))
        stop()
      } else if (length(grep("\\.obidms",tax_dir)) == 0) {
        message(paste0("The taxonomic db '", tax_dir,"' must be an .obidms object. Please check the path."))
        stop()
      } else {
        tax_dms_name <- gsub("\\.obidms", "", basename(tax_dir[grep("\\.obidms",tax_dir)]))
        tax_dir <- tax_db_dir
      }

      
      # it is necessary to run ecotag within a new directory for each part.
      # this is because dms can not be calles from two processes at the same time
      # and can not change the name of the dms so make a copy in each directory and
      # run there the ecotag
      X <- NULL
      create_dirs <- NULL
      for (i in 1:cores) {
        if (cores == 1) {
          fasta_file_name <- paste0(experiment,"_ODIN.fasta")
        } else {
          fasta_file_name <- paste0(experiment, "_THOR_pretaxa_part_",sprintf("%02d",i),".fasta")
        }
        create_dirs <- c(create_dirs,paste0("mkdir ",experiment, "_THOR_",sprintf("%02d",i)," ; cp -r ",tax_dir,"/",tax_dms_name,".obidms ",experiment, "_THOR_",sprintf("%02d",i),"/. ; "))
        X <- c(X,paste0("cd ",experiment, "_THOR_",sprintf("%02d",i)," ; ",
                        "obi import --fasta-input ../",fasta_file_name, " ",tax_dms_name,"/seqs ; ",
                        "obi ecotag -c ",minimum_circle," --taxonomy ",tax_dms_name,"/taxonomy/my_tax -R ",tax_dms_name,"/ref_db ",tax_dms_name,"/seqs ",tax_dms_name,"/assigned_seqs ; ",
                        "obi annotate --taxonomy ",tax_dms_name,"/taxonomy/my_tax ",
                        " --with-taxon-at-rank superkingdom ",
                        " --with-taxon-at-rank kingdom ",
                        " --with-taxon-at-rank phylum ",
                        " --with-taxon-at-rank class ",
                        " --with-taxon-at-rank order ",
                        " --with-taxon-at-rank family ",
                        " --with-taxon-at-rank genus ",
                        " --with-taxon-at-rank species ",
                        tax_dms_name,"/assigned_seqs ",tax_dms_name,"/assigned_seqs_add ; ",
                        "obi export --fasta-output ",tax_dms_name,"/assigned_seqs_add >../",experiment,"_THOR_part_",sprintf("%02d",i),".fasta"))
      }
      
      suppressPackageStartupMessages(library(parallel))
      mclapply(create_dirs, function(x) system(x,intern=T,wait=T), mc.cores = cores)
      mclapply(X, function(x) system(x,intern=T,wait=T), mc.cores = cores)
    }
    
    message("THOR will add higher taxonomic ranks now.")
    filefasta <-paste0(experiment,"_THOR.fasta")
    system(paste0("cat ",experiment,"_THOR_part_??.fasta | grep '>' > ",filefasta),intern=T,wait=T)
    outfile <- paste0(experiment,"_THOR_annotated.tsv")
    # Here old owi_add_taxonomy starts
    suppressPackageStartupMessages(library("tidyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("parallel"))
    length_id <- 14 # This is the total length of the MOTU IDs in filefasta. It can be changed if needed.
    # Load external information on higher taxa:
    taxo_names <- read.table(paste0(tax_dir,"/","order.complete.csv"),sep=",",head=T)
    class_to_sk <- unique(taxo_names[,2:5])
    phylum_to_sk <- unique(class_to_sk[,2:4])
    kingdom_to_sk <-unique(phylum_to_sk[,2:3])
    # showlines <- 1
    
    family_to_order <- read.table(paste0(tax_dir,"/","family_to_order.csv"),sep=",",head=T)
    genus_to_family <- read.table(paste0(tax_dir,"/","genus_to_family.csv"),sep=",",head=T)
    
    # The following are taxonomic ranks that can be found in the "scientific_name" field and should be treated as exceptions, by function fix_exceptions.
    exceptions <- c("Acanthomorphata","Acari","Acariformes","Aciculata","Actinopterygii","Agaricomycetesincertaesedis","Agaricomycetidae",
                    "Agaricomycotina","Alveolata","Amniota","Amoebozoa","Anystina","Apusozoa","Armophorea","asterids","Asterozoa","Bacillariophycidae",
                    "Bacteria","Batoidea","Bdelloidea","betaproteobacteriumCB","Biddulphiophycidae","Bilateria","Boreoeutheria","Bryophyta",
                    "Buccinoidea","Caenogastropoda","campanulids","Centroheliozoa","Cephalaspidea","Cercozoa","Chelicerata","Choreotrichia",
                    "Chrysophyceaesp.","Ciliophora","CirripFedia","Clitellata","Clupeocephala","Collembola","Corallinophycidae","Coscinodiscophycidae",
                    "Crustacea","Cryptophyta","Ctenosquamata","Cymbellales","Dactylopodida","Decapodiformes","Deuterostomia","Digenea","Dikarya",
                    "Dorylaimia","Ecdysozoa","Echinacea","Echinozoa","Eleutherozoa","Embryophyta","Endopterygota","Enoplia","Entomobryomorpha",
                    "Euacanthomorphacea","Eucarida","eudicotyledons","Euechinoidea","Euheterodonta","Eumalacostraca","Eumetazoa","Euopisthobranchia",
                    "Eupercaria","Euphyllophyta","Eurotiomycetidae","Euteleosteomorpha","Euteleosteomorpha","Euteleostomi","Euthyneura",
                    "Fragilariophycidae","Fungiincertaesedis","Galeoidea","Gnathostomata","Gregarinasina","Haptophyceae","Heterobranchia",
                    "Heteroconchia","Heteroscleromorpha","Hexacorallia","Hexanauplia","Hexapoda","Hoplocarida","Hydracarina","Hydroidolina",
                    "hypocreomyceta","Hypocreomycetidae","Hypsogastropoda","Imparidentia","Intramacronucleata","Jakobida","Labyrinthulomycetes",
                    "lamiids","Lecanoromycetidae","leotiomyceta","Leotiomycetesincertaesedis","Leotiomycetidae","Littorinimorpha","Littorinoidea",
                    "Lophotrochozoa","Magnoliophyta","Mandibulata","Mesangiospermae","Moniliformopses","Myxogastria","Nemaliophycidae","Neocopepoda",
                    "Neogastropoda","Neoloricata","Neoptera","Nudipleura","Ochrophyta","Octocorallia","Oligochaeta","Oligotrichia","Ophiuridea",
                    "Opisthokonta","Palaeoheterodonta","Palaeonemertea","Paleoptera","Palpata","Panarthropoda","Pancrustacea","Panpulmonata",
                    "Paraneoptera","Patellogastropoda","Pentapetalae","Peracarida","Percomorphaceae","Peritrichia","Peritrichia","Petrosaviidae",
                    "Pezizomycotina","Phascolosomatidea","Picobiliphytesp.MS584-11","Pinidae","Piroplasmida","Podoplea","Poduroidea","Poduromorpha",
                    "Prostigmata","Protostomia","Prymnesiophyceae","Pteriomorphia","Pterygota","Pucciniomycotina","PXclade","Rhabditophora","Rhodophyta",
                    "Rhodymeniophycidae","rosids","saccharomyceta","Sacoglossa","Sar","Sarcoptiformes","Scleractinia","Scolecida","Scuticociliatia",
                    "Sedentaria","Seriata","Silicofilosea","Sipuncula","sordariomyceta","Sordariomycetidae","Spermatophyta","Stichotrichia","Stramenopiles",
                    "Streptophytina","Stylommatophora","Symphypleona","Synurophyceae","Taphrinomycotina","Teleostei","Thoracica","Trachylinae",
                    "Trichostomatia","Trombidiformes","Tubulinea","Tubulinida","Tunicata","unculturedactinobacterium","unculturedbacterium",
                    "Ustilaginomycotina","Vetigastropoda","Xanthophyceae"
    )
    
    fix_exceptions <- function(scientific_name){
      if (scientific_name %in% c("asterids","campanulids","lamiids","rosids","Pentapetalae","Mesangiospermae","eudicotyledons","Magnoliophyta")) {
        matrix.data["class_name",2] <- "Magnoliopsida"
        matrix.data["phylum_name",2] <- "Tracheophyta"
        matrix.data["kingdom_name",2] <- "Viridiplantae"
        matrix.data["superkingdom_name",2] <- "Archaeplastida"
      }
      if (scientific_name %in% c("Pinidae")) {
        matrix.data["class_name",2] <- "Pinopsida"
        matrix.data["phylum_name",2] <- "Tracheophyta"
        matrix.data["kingdom_name",2] <- "Viridiplantae"
        matrix.data["superkingdom_name",2] <- "Archaeplastida"
      }
      if (scientific_name %in% c("Moniliformopses","Euphyllophyta")) {
        matrix.data["class_name",2] <- NA
        matrix.data["phylum_name",2] <- "Tracheophyta"
        matrix.data["kingdom_name",2] <- "Viridiplantae"
        matrix.data["superkingdom_name",2] <- "Archaeplastida"
      }
      if (scientific_name %in% c("Bacteria","uncultured bacterium","uncultured actinobacterium")) {
        matrix.data["kingdom_name",2] <- "Bacteria"
        matrix.data["superkingdom_name",2] <- "Prokaryota"
      }
      if (scientific_name %in% c("beta proteobacterium CB")) {
        matrix.data["class_name",2] <- "Betaproteobacteria"
        matrix.data["phylum_name",2] <- "Proteobacteria"
        matrix.data["kingdom_name",2] <- "Bacteria"
        matrix.data["superkingdom_name",2] <- "Prokaryota"
      }
      if (scientific_name %in% c("Stichotrichia","Peritrichia","Oligotrichia","Choreotrichia","Ciliophora","Trichostomatia")) {
        matrix.data["class_name",2] <- NA
        matrix.data["phylum_name",2] <- "Ciliophora"
        matrix.data["kingdom_name",2] <- "Alveolata"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name %in% c("Petrosaviidae")) {
        matrix.data["class_name",2] <- "Liliopsida"
        matrix.data["phylum_name",2] <- "Tracheophyta"
        matrix.data["kingdom_name",2] <- "Viridiplantae"
        matrix.data["superkingdom_name",2] <- "Archaeplastida"
      }
      if (scientific_name %in% c("Centroheliozoa")) {
        matrix.data["class_name",2] <- "Centrohelea"
        matrix.data["phylum_name",2] <- "Heliozoa"
        matrix.data["kingdom_name",2] <- "Hacrobia"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name %in% c("Sar")) {
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name %in% c("Cryptophyta")) {
        matrix.data["class_name",2] <- NA
        matrix.data["phylum_name",2] <- "Cryptophyta"
        matrix.data["kingdom_name",2] <- "Hacrobia"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name %in% c("Prymnesiophyceae")) {
        matrix.data["class_name",2] <- "Prymnesiophyceae"
        matrix.data["phylum_name",2] <- "Haptophyta"
        matrix.data["kingdom_name",2] <- "Hacrobia"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name %in% c("Picobiliphyte sp. MS584-11")) {
        matrix.data["class_name",2] <- "Picomonadea"
        matrix.data["phylum_name",2] <- "Picozoa"
        matrix.data["kingdom_name",2] <- "Hacrobia"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name %in% c("Chrysophyceae sp.")) {
        matrix.data["class_name",2] <- "Chrysophyceae"
        matrix.data["phylum_name",2] <- "Ochrophyta"
        matrix.data["kingdom_name",2] <- "Stramenopiles"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name %in% c("Ochrophyta")) {
        matrix.data["phylum_name",2] <- "Ochrophyta"
        matrix.data["kingdom_name",2] <- "Stramenopiles"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
        matrix.data["rank",2] <- "phylum"
      }
      if (scientific_name %in% c("Cymbellales")) {
        matrix.data["class_name",2] <- "Bacillariophyceae"
        matrix.data["phylum_name",2] <- "Bacillariophyta"
        matrix.data["kingdom_name",2] <- "Stramenopiles"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name %in% c("Coscinodiscophycidae")) {
        matrix.data["class_name",2] <- "Coscinodiscophyceae"
        matrix.data["phylum_name",2] <- "Bacillariophyta"
        matrix.data["kingdom_name",2] <- "Stramenopiles"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name %in% c("Biddulphiophycidae")) {
        matrix.data["class_name",2] <- "Mediophyceae"
        matrix.data["phylum_name",2] <- "Bacillariophyta"
        matrix.data["kingdom_name",2] <- "Stramenopiles"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name %in% c("Enoplia","Dorylaimia")) {
        matrix.data["class_name",2] <- "Enoplea"
        matrix.data["phylum_name",2] <- "Nematoda"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Rhabditophora","Seriata")) {
        matrix.data["class_name",2] <- "Rhabditophora"
        matrix.data["phylum_name",2] <- "Platyhelminthes"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Gregarinasina")) {
        matrix.data["class_name",2] <- "Conoidasida"
        matrix.data["phylum_name",2] <- "Apicomplexa"
        matrix.data["kingdom_name",2] <- "Alveolata"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name %in% c("Tubulinea","Tubulinida")) {
        matrix.data["class_name",2] <- "Tubulinea"
        matrix.data["phylum_name",2] <- "Tubulinea"
        matrix.data["kingdom_name",2] <- "Lobosa"
        matrix.data["superkingdom_name",2] <- "Amoebozoa"
      }
      if (scientific_name %in% c("Dactylopodida")) {
        matrix.data["class_name",2] <- "Flabellinia"
        matrix.data["phylum_name",2] <- "Discosea"
        matrix.data["kingdom_name",2] <- "Lobosa"
        matrix.data["superkingdom_name",2] <- "Amoebozoa"
      }
      if (scientific_name %in% c( "Amoebozoa","Apusozoa")) {
        matrix.data["superkingdom_name",2] <- scientific_name
      }
      if (scientific_name %in% c("Pezizomycotina","Taphrinomycotina")) {
        matrix.data["phylum_name",2] <- "Ascomycota"
        matrix.data["kingdom_name",2] <- "Fungi"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("saccharomyceta")) {
        matrix.data["class_name",2] <- "Saccharomycetes"
        matrix.data["phylum_name",2] <- "Ascomycota"
        matrix.data["kingdom_name",2] <- "Fungi"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Lecanoromycetidae")) {
        matrix.data["class_name",2] <- "Lecanoromycetes"
        matrix.data["phylum_name",2] <- "Ascomycota"
        matrix.data["kingdom_name",2] <- "Fungi"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Agaricomycotina","Pucciniomycotina","Ustilaginomycotina")) {
        matrix.data["phylum_name",2] <- "Basidiomycota"
        matrix.data["kingdom_name",2] <- "Fungi"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Agaricomycetes incertae sedis","Agaricomycetidae")) {
        matrix.data["class_name",2] <- "Agaricomycetes"
        matrix.data["phylum_name",2] <- "Basidiomycota"
        matrix.data["kingdom_name",2] <- "Fungi"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Bryophyta")) {
        matrix.data["phylum_name",2] <- "Bryophyta"
        matrix.data["kingdom_name",2] <- "Viridiplantae"
        matrix.data["superkingdom_name",2] <- "Archaeplastida"
      }
      if (scientific_name %in% c("Armophorea")) {
        matrix.data["class_name",2] <- "Armophorea"
        matrix.data["phylum_name",2] <- "Ciliophora"
        matrix.data["kingdom_name",2] <- "Alveolata"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name == "Cercozoa") {
        matrix.data["class_name",2] <- NA
        matrix.data["phylum_name",2] <- "Cercozoa"
        matrix.data["kingdom_name",2] <- "Rhizaria"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name == "Silicofilosea") {
        matrix.data["class_name",2] <- "Imbricatea"
        matrix.data["phylum_name",2] <- "Cercozoa"
        matrix.data["kingdom_name",2] <- "Rhizaria"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name == "Jakobida") {
        matrix.data["order_name",2] <- "Jakobida"
        matrix.data["class_name",2] <- "Jakobea"
        matrix.data["phylum_name",2] <- "Loukozoa"
        matrix.data["kingdom_name",2] <- "Loukozoa"
        matrix.data["superkingdom_name",2] <- "Excavata"
      }
      if (scientific_name == "Piroplasmida") {
        matrix.data["order_name",2] <- "Piroplasmorida"
        matrix.data["class_name",2] <- "Aconoidasida"
        matrix.data["phylum_name",2] <- "Apicomplexa"
        matrix.data["kingdom_name",2] <- "Alveolata"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name == "Intramacronucleata") {
        matrix.data["class_name",2] <- "Intramacronucleata"
        matrix.data["phylum_name",2] <- "Ciliophora"
        matrix.data["kingdom_name",2] <- "Alveolata"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name == "Scuticociliatia") {
        matrix.data["class_name",2] <- "Oligohymenophorea"
        matrix.data["phylum_name",2] <- "Ciliophora"
        matrix.data["kingdom_name",2] <- "Alveolata"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
      }
      if (scientific_name == "Heteroscleromorpha") {
        matrix.data["class_name",2] <- "Demospongiae"
        matrix.data["phylum_name",2] <- "Porifera"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name == "Palaeonemertea") {
        matrix.data["order_name",2] <- "Palaeonemertea"
        matrix.data["class_name",2] <- "Anopla"
        matrix.data["phylum_name",2] <- "Nemertea"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
        matrix.data["rank",2] <- "order"
      }
      if (scientific_name %in%  c("Phascolosomatidea")) {
        matrix.data["order_name",2] <- "Phascolosomatiformes"
        matrix.data["class_name",2] <- "Sipuncula"
        matrix.data["phylum_name",2] <- "Annelida"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in%  c("Sipuncula")) {
        matrix.data["class_name",2] <- "Sipuncula"
        matrix.data["phylum_name",2] <- "Annelida"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in%  c("Bdelloidea")) {
        matrix.data["class_name",2] <- "Eurotatoria"
        matrix.data["phylum_name",2] <- "Rotifera"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name == "Digenea") {
        matrix.data["class_name",2] <- "Trematoda"
        matrix.data["phylum_name",2] <- "Platyhelminthes"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name == "Alveolata") {
        matrix.data["class_name",2] <- NA
        matrix.data["phylum_name",2] <- NA
        matrix.data["kingdom_name",2] <- "Alveolata"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
        matrix.data["rank",2] <- "kingdom"
      }
      if (scientific_name == "Haptophyceae") {
        matrix.data["class_name",2] <- "Haptophyceae"
        matrix.data["phylum_name",2] <- "Haptophyta"
        matrix.data["kingdom_name",2] <- "Hacrobia"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
        matrix.data["rank",2] <- "class"
      }
      if (scientific_name == "Synurophyceae") {
        matrix.data["class_name",2] <- "Chrysophyceae"
        matrix.data["phylum_name",2] <- "Ochrophyta"
        matrix.data["kingdom_name",2] <- "Stramenopiles"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
        matrix.data["rank",2] <- "class"
      }
      if (scientific_name == "Labyrinthulomycetes") {
        matrix.data["class_name",2] <- "Labyrinthulea"
        matrix.data["phylum_name",2] <- "Bigyra"
        matrix.data["kingdom_name",2] <- "Stramenopiles"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
        matrix.data["rank",2] <- "class"
      }
      if (scientific_name == "Xanthophyceae") {
        matrix.data["class_name",2] <- "Xanthophyceae"
        matrix.data["phylum_name",2] <- "Ochrophyta"
        matrix.data["kingdom_name",2] <- "Stramenopiles"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
        matrix.data["rank",2] <- "class"
      }
      if (scientific_name %in% c("Fungi incertae sedis","Dikarya")) {
        matrix.data["class_name",2] <- NA
        matrix.data["phylum_name",2] <- NA
        matrix.data["kingdom_name",2] <- "Fungi"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
        matrix.data["rank",2] <- "kingdom"
      }
      if (scientific_name == "Opisthokonta") {
        matrix.data["class_name",2] <- NA
        matrix.data["phylum_name",2] <- NA
        matrix.data["kingdom_name",2] <- NA
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
        matrix.data["rank",2] <- "superkingdom"
      }
      if (scientific_name %in% c("Streptophytina","Embryophyta")) {
        matrix.data["class_name",2] <- NA
        matrix.data["phylum_name",2] <- "Streptophyta"
        matrix.data["kingdom_name",2] <- "Viridiplantae"
        matrix.data["superkingdom_name",2] <- "Archaeplastida"
      }
      if (scientific_name %in% c("Spermatophyta")) {
        matrix.data["class_name",2] <- NA
        matrix.data["phylum_name",2] <- NA
        matrix.data["kingdom_name",2] <- "Viridiplantae"
        matrix.data["superkingdom_name",2] <- "Archaeplastida"
      }
      if (scientific_name %in% c("Bilateria","Ecdysozoa","Eumetazoa","Lophotrochozoa","Panarthropoda","Protostomia","Deuterostomia")) {
        matrix.data["class_name",2] <- NA
        matrix.data["phylum_name",2] <- NA
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Scolecida","Aciculata","Palpata","Sedentaria")) {
        matrix.data["class_name",2] <- "Polychaeta"
        matrix.data["phylum_name",2] <- "Annelida"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Clitellata")) {
        matrix.data["class_name",2] <- "Clitellata"
        matrix.data["phylum_name",2] <- "Annelida"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Oligochaeta")) {
        matrix.data["class_name",2] <- "Oligochaeta"
        matrix.data["phylum_name",2] <- "Annelida"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Eleutherozoa","Echinozoa","Asterozoa")) {
        matrix.data["class_name",2] <- NA
        matrix.data["phylum_name",2] <- "Echinodermata"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Euechinoidea","Echinacea","Gnathostomata")) {
        matrix.data["class_name",2] <- "Echinoidea"
        matrix.data["phylum_name",2] <- "Echinodermata"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Ophiuridea")) {
        matrix.data["class_name",2] <- "Ophiuroidea"
        matrix.data["phylum_name",2] <- "Echinodermata"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name == "Rhodophyta") {
        matrix.data["phylum_name",2] <- "Rhodophyta"
        matrix.data["kingdom_name",2] <- "Rhodophyta"
        matrix.data["superkingdom_name",2] <- "Archaeplastida"
        matrix.data["rank",2] <- "phylum"
      }
      if (scientific_name %in% c("Corallinophycidae","Rhodymeniophycidae","Nemaliophycidae")) {
        matrix.data["class_name",2] <- "Florideophyceae"
        matrix.data["phylum_name",2] <- "Rhodophyta"
        matrix.data["kingdom_name",2] <- "Rhodophyta"
        matrix.data["superkingdom_name",2] <- "Archaeplastida"
        matrix.data["rank",2] <- "phylum"
      }
      if (scientific_name == "PX clade") {
        matrix.data["phylum_name",2] <- "Ochrophyta"
        matrix.data["kingdom_name",2] <- "Stramenopiles"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
        matrix.data["rank",2] <- "phylum"
      }
      if (scientific_name == "Stramenopiles") {
        matrix.data["kingdom_name",2] <- "Stramenopiles"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
        matrix.data["rank",2] <- "kingdom"
      }
      if (scientific_name %in% c("Euteleostomi","Teleostei","Euteleosteomorpha","Actinopterygii","Clupeocephala","Eupercaria","Percomorphaceae",
                                 "Acanthomorphata","Euacanthomorphacea","Ctenosquamata","Euteleosteomorpha")) {
        matrix.data["class_name",2] <- "Actinopterygii"
        matrix.data["phylum_name",2] <- "Chordata"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Batoidea","Galeoidea")) {
        matrix.data["class_name",2] <- "Chondrichthyes"
        matrix.data["phylum_name",2] <- "Chordata"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Boreoeutheria")) {
        matrix.data["class_name",2] <- "Mammalia"
        matrix.data["phylum_name",2] <- "Chordata"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Amniota","Tunicata")) {
        matrix.data["phylum_name",2] <- "Chordata"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Chelicerata","Crustacea","Pancrustacea","Mandibulata")) {
        matrix.data["class_name",2] <- NA
        matrix.data["phylum_name",2] <- "Arthropoda"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Eumalacostraca","Peracarida","Eucarida","Hoplocarida")) {
        matrix.data["class_name",2] <- "Malacostraca"
        matrix.data["phylum_name",2] <- "Arthropoda"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Neocopepoda","Podoplea","Cirripedia","Thoracica","Hexanauplia")) {
        matrix.data["class_name",2] <- "Maxillopoda"
        matrix.data["phylum_name",2] <- "Arthropoda"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Collembola","Entomobryomorpha","Poduromorpha","Symphypleona","Poduroidea")) {
        matrix.data["order_name",2] <- "Collembola"
        matrix.data["class_name",2] <- "Hexapoda"
        matrix.data["phylum_name",2] <- "Arthropoda"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
        matrix.data["rank",2] <- "order"
      }
      if (scientific_name %in% c("Hexapoda")) {
        matrix.data["class_name",2] <- "Hexapoda"
        matrix.data["phylum_name",2] <- "Arthropoda"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
        matrix.data["rank",2] <- "class"
      }
      if (scientific_name %in% c("Acariformes","Acari")) {
        matrix.data["class_name",2] <- "Arachnida"
        matrix.data["phylum_name",2] <- "Arthropoda"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
        matrix.data["rank",2] <- "order"
      }
      if (scientific_name=="Sarcoptiformes") {
        matrix.data["order_name",2] <- "Sarcoptiformes"
        matrix.data["class_name",2] <- "Arachnida"
        matrix.data["phylum_name",2] <- "Arthropoda"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
        matrix.data["rank",2] <- "order"
      }
      if (scientific_name %in% c("Trombidiformes","Anystina","Hydracarina","Prostigmata")) {
        matrix.data["order_name",2] <- "Trombidiformes"
        matrix.data["class_name",2] <- "Arachnida"
        matrix.data["phylum_name",2] <- "Arthropoda"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
        matrix.data["rank",2] <- "order"
      }
      if (scientific_name %in% c("Endopterygota","Neoptera","Pterygota","Paleoptera","Paraneoptera")) {
        matrix.data["class_name",2] <- "Insecta"
        matrix.data["phylum_name",2] <- "Arthropoda"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Stylommatophora","Euthyneura","Nudipleura","Panpulmonata","Heterobranchia",
                                 "Euopisthobranchia","Hypsogastropoda","Caenogastropoda")) {
        matrix.data["class_name",2] <- "Gastropoda"
        matrix.data["phylum_name",2] <- "Mollusca"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in%  c("Littorinimorpha","Littorinoidea")) {
        matrix.data["order_name",2] <- "Littorinimorpha"
        matrix.data["class_name",2] <- "Gastropoda"
        matrix.data["phylum_name",2] <- "Mollusca"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in%  c("Sacoglossa")) {
        matrix.data["order_name",2] <- "Sacoglossa"
        matrix.data["class_name",2] <- "Gastropoda"
        matrix.data["phylum_name",2] <- "Mollusca"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
        matrix.data["rank",2] <- "order"
      }
      if (scientific_name %in%  c("Patellogastropoda")) {
        matrix.data["order_name",2] <- "Patellogastropoda"
        matrix.data["class_name",2] <- "Gastropoda"
        matrix.data["phylum_name",2] <- "Mollusca"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in%  c("Vetigastropoda")) {
        matrix.data["order_name",2] <- "Vetigastropoda"
        matrix.data["class_name",2] <- "Gastropoda"
        matrix.data["phylum_name",2] <- "Mollusca"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in%  c("Neogastropoda","Buccinoidea")) {
        matrix.data["order_name",2] <- "Neogastropoda"
        matrix.data["class_name",2] <- "Gastropoda"
        matrix.data["phylum_name",2] <- "Mollusca"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in%  c("Neoloricata")) {
        matrix.data["class_name",2] <- "Polyplacophora"
        matrix.data["phylum_name",2] <- "Mollusca"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in%  c("Cephalaspidea")) {
        matrix.data["order_name",2] <- "Cephalaspidea"
        matrix.data["class_name",2] <- "Gastropoda"
        matrix.data["phylum_name",2] <- "Mollusca"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in%  c("Pteriomorphia","Imparidentia","Heteroconchia","Euheterodonta","Palaeoheterodonta")) {
        matrix.data["class_name",2] <- "Bivalvia"
        matrix.data["phylum_name",2] <- "Mollusca"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Decapodiformes")) {
        matrix.data["class_name",2] <- "Cephalopoda"
        matrix.data["phylum_name",2] <- "Mollusca"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name == "Myxogastria") {
        matrix.data["class_name",2] <- "Myxomycetes"
        matrix.data["phylum_name",2] <- "Mycetozoa"
        matrix.data["kingdom_name",2] <- "Conosa"
        matrix.data["superkingdom_name",2] <- "Amoebozoa"
        matrix.data["rank",2] <- "class"
        matrix.data["scientific_name",2] <- "Myxomycetes"
      }
      if (scientific_name == "Fragilariophycidae") {
        matrix.data["class_name",2] <- "Fragilariophyceae"
        matrix.data["phylum_name",2] <- "Bacillariophyta"
        matrix.data["kingdom_name",2] <- "Stramenopiles"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
        matrix.data["rank",2] <- "class"
        matrix.data["scientific_name",2] <- "Bacillariophyceae"
      }
      if (scientific_name == "Bacillariophycidae") {
        matrix.data["class_name",2] <- "Bacillariophyceae"
        matrix.data["phylum_name",2] <- "Bacillariophyta"
        matrix.data["kingdom_name",2] <- "Stramenopiles"
        matrix.data["superkingdom_name",2] <- "Chromalveolata"
        matrix.data["rank",2] <- "class"
        matrix.data["scientific_name",2] <- "Bacillariophyceae"
      }
      if (scientific_name %in% c("Octocorallia","Hexacorallia","Scleractinia")) {
        matrix.data["class_name",2] <- "Anthozoa"
        matrix.data["phylum_name",2] <- "Cnidaria"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Hydroidolina")) {
        matrix.data["class_name",2] <- "Hydrozoa"
        matrix.data["phylum_name",2] <- "Cnidaria"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("Trachylinae")) {
        matrix.data["order_name",2] <- "Trachylina"
        matrix.data["class_name",2] <- "Hydrozoa"
        matrix.data["phylum_name",2] <- "Cnidaria"
        matrix.data["kingdom_name",2] <- "Metazoa"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("sordariomyceta","Sordariomycetidae")) {
        matrix.data["class_name",2] <- "Sordariomycetes"
        matrix.data["phylum_name",2] <- "Ascomycota"
        matrix.data["kingdom_name",2] <- "Fungi"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("leotiomyceta","Leotiomycetidae","Leotiomycetes incertae sedis")) {
        matrix.data["class_name",2] <- "Leotiomycetes"
        matrix.data["phylum_name",2] <- "Ascomycota"
        matrix.data["kingdom_name",2] <- "Fungi"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("eurotiomyceta","Eurotiomycetidae")) {
        matrix.data["class_name",2] <- "Eurotiomycetes"
        matrix.data["phylum_name",2] <- "Ascomycota"
        matrix.data["kingdom_name",2] <- "Fungi"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      if (scientific_name %in% c("hypocreomyceta","Hypocreomycetidae")) {
        matrix.data["class_name",2] <- "Sordariomycetes"
        matrix.data["phylum_name",2] <- "Ascomycota"
        matrix.data["kingdom_name",2] <- "Fungi"
        matrix.data["superkingdom_name",2] <- "Opisthokonta"
      }
      return(matrix.data)
    }
    
    #Read  fasta file
    message("Reading ecotagged fasta file")
    info <- readLines(filefasta)
    info <- gsub(">", "", info)
    #Initialize
    inicio <- T
    message("Read ",length(info)," records")
    
    complete_taxa <- function(fila,info,genus_to_family,family_to_order,class_to_sk,phylum_to_sk,kingdom_to_sk,exceptions) {
      
      infofasta <- info[fila]
      # if order_name is repeated remove one
      if (length(gregexpr("order_name",infofasta)[[1]])==2) infofasta <- sub(pattern="; order_name=",replacement="",x=infofasta)
      # if order is repeated remove one
      if (length(gregexpr("order=",infofasta)[[1]])==2) infofasta <- sub(pattern="order=",replacement="",x=infofasta)
      # get id of the sequence
      id <- strsplit(infofasta,split=c(" "))[[1]][1]
      # get the rest of the info and split it
      infofasta2 <- substr(infofasta,nchar(id)+2,nchar(infofasta))
      infofasta2 <- strsplit(sub(pattern=";'",replacement="'",x=infofasta2),";")
      infofasta2 <- gsub("^ ", "", unlist(infofasta2))
      infofasta2 <- infofasta2[grepl("=",infofasta2)]
      matrix.data <- data.frame(X1=infofasta2)
      matrix.data <- separate(matrix.data,X1,into = c("X1","X2"),sep = "=")
      
      columns <- c("id","rank","SCIENTIFIC_NAME","BEST_IDENTITY",
                   "superkingdom_name","kingdom_name","phylum_name","class_name",
                   "order_name","family_name","genus_name","species_name",
                   "BEST_MATCH_IDS","TAXID")
      
      matrix.data <- suppressMessages(left_join(data.frame(X1=columns),matrix.data))
      rownames(matrix.data) <- matrix.data$X1
      matrix.data["id",2] <- substr(id,1,length_id)
      
      db.out <- data.frame(matrix(nrow = 1,ncol = length(columns)))
      colnames(db.out) <- columns
      
      # if only genus is filled, if the current genus is found in the genus_to_family.csv fill the family and order
      # else add "correct manually"
      if (is.na(matrix.data["order_name",2]) & is.na(matrix.data["family_name",2]) & !is.na(matrix.data["genus_name",2])){
        if (matrix.data["genus_name",2] %in% genus_to_family$genus_name){
          matrix.data["family_name",2] <- as.character(genus_to_family$family_name[genus_to_family$genus_name==matrix.data["genus_name",2]])
          matrix.data["order_name",2] <- as.character(family_to_order$order_name[family_to_order$family_name==matrix.data["family_name",2]])
        } else {
          matrix.data["order_name",2] <- "Correct_manually"
          matrix.data["family_name",2] <- "Correct_manually"
        }
      }
      
      # if only family is filled, if the current family is found in the family_to_order.csv fill order
      # else add "correct manually"
      if (is.na(matrix.data["order_name",2]) & !is.na(matrix.data["family_name",2])){
        if (matrix.data["family_name",2] %in% family_to_order$family_name){
          matrix.data["order_name",2] <- as.character(family_to_order$order_name[family_to_order$family_name==matrix.data["family_name",2]])
        } else {
          matrix.data["order_name",2] <- "Correct_manually"
        }
      }
      
      # add higher than order ranks for those orders that the user has delimited this info
      # in the taxo_names
      if (matrix.data["order_name",2]%in%taxo_names$order_name) {
        matrix.data[match(c("superkingdom_name"),rownames(matrix.data)),2] <- taxo_names[taxo_names$order_name==matrix.data["order_name",2],c("superkingdom_name")][1]
        matrix.data[match(c("kingdom_name"),rownames(matrix.data)),2] <- taxo_names[taxo_names$order_name==matrix.data["order_name",2],c("kingdom_name")][1]
        matrix.data[match(c("phylum_name"),rownames(matrix.data)),2] <- taxo_names[taxo_names$order_name==matrix.data["order_name",2],c("phylum_name")][1]
        matrix.data[match(c("class_name"),rownames(matrix.data)),2] <- taxo_names[taxo_names$order_name==matrix.data["order_name",2],c("class_name")][1]
      }
      
      if (!is.na(matrix.data["species_name",2])) {
        matrix.data["rank",2] <- "species"
      } else if (!is.na(matrix.data["genus_name",2])) {
        matrix.data["rank",2] <- "genus"
      } else if (!is.na(matrix.data["family_name",2])) {
        matrix.data["rank",2] <- "family"
      } else if (!is.na(matrix.data["order_name",2])) {
        matrix.data["rank",2] <- "order"
      } else if (!is.na(matrix.data["class_name",2])) {
        matrix.data["rank",2] <- "class"
      } else if (!is.na(matrix.data["phylum_name",2])) {
        matrix.data["rank",2] <- "phylum"
      } else if (!is.na(matrix.data["kingdom_name",2])) {
        matrix.data["rank",2] <- "kingdom"
      } else if (!is.na(matrix.data["superkingdom_name",2])) {
        matrix.data["rank",2] <- "superkingdom"
      } else {
        matrix.data["rank",2] <- "root"
      }
      taxa <- c("superkingdom","kingdom","phylum","class",
                "order","family","genus","species")
      taxa_name <- c("superkingdom_name","kingdom_name","phylum_name","class_name",
                     "order_name","family_name","genus_name","species_name")
      # trace those taxa with na info in higher taxa than assigned. If root there's nothing higher
      if(matrix.data["rank",2]=="root") {
        na_taxa <- 0
      } else {
        na_taxa <- taxa_name[1:which(matrix.data["rank",2]==taxa)]
      }
      na_taxa <- na_taxa[is.na(matrix.data[na_taxa,2])]
      
      matrix.data[na_taxa,2] <- "Correct_manually"
      
      # inspect if higher rank than order can be filled if they haven't in previous steps
      
      if (length(na_taxa)>0) {
        if (("phylum_name" %in% na_taxa) & !("class_name" %in% na_taxa)) {
          if (matrix.data["class_name",2] %in% taxo_names$class_name) {
            matrix.data[match(c("superkingdom_name"),rownames(matrix.data)),2] <- taxo_names[taxo_names$class_name==matrix.data["class_name",2],c("superkingdom_name")][1]
            matrix.data[match(c("kingdom_name"),rownames(matrix.data)),2] <- taxo_names[taxo_names$class_name==matrix.data["class_name",2],c("kingdom_name")][1]
            matrix.data[match(c("phylum_name"),rownames(matrix.data)),2] <- taxo_names[taxo_names$class_name==matrix.data["class_name",2],c("phylum_name")][1]
          }
        } else if (("kingdom_name" %in% na_taxa) & !("phylum_name" %in% na_taxa)) {
          if (matrix.data["phylum_name",2] %in% taxo_names$phylum_name) {
            matrix.data[match(c("superkingdom_name"),rownames(matrix.data)),2] <- taxo_names[taxo_names$phylum_name==matrix.data["phylum_name",2],c("superkingdom_name")][1]
            matrix.data[match(c("kingdom_name"),rownames(matrix.data)),2] <- taxo_names[taxo_names$phylum_name==matrix.data["phylum_name",2],c("kingdom_name")][1]
          }
        } else if (("superkingdom_name" %in% na_taxa) & !("kingdom_name" %in% na_taxa)) {
          if (matrix.data["kingdom_name",2] %in% taxo_names$kingdom_name) {
            matrix.data[c("superkingdom_name"),2] <- taxo_names[taxo_names$kingdom_name==matrix.data["kingdom_name",2],c("superkingdom_name")][1]
          }
        }
      }
      
      if (matrix.data["scientific_name",2] %in% exceptions) matrix.data <- fix_exceptions(matrix.data["scientific_name",2],matrix.data = matrix.data)
      
      # make sure superkingdom is correct
      if (matrix.data["kingdom_name",2] %in% kingdom_to_sk$kingdom_name) {
        matrix.data["superkingdom_name",2] <- kingdom_to_sk$superkingdom_name[kingdom_to_sk$kingdom_name==matrix.data["kingdom_name",2]]
      }
      
      db.out[1,] <- matrix.data[match(columns,rownames(matrix.data)),2]
      
      # print(db.out) }
      return(db.out)
    }
    
    db.out <- do.call("rbind", mclapply(seq_along(info), complete_taxa,
                                       info=info,genus_to_family=genus_to_family,family_to_order=family_to_order,
                                       class_to_sk=class_to_sk,phylum_to_sk=phylum_to_sk,kingdom_to_sk=kingdom_to_sk,exceptions=exceptions,
                                       mc.cores = cores))
    
    # Substitute commas by "|" in species_list and remove special characters
    db.out$BEST_MATCH_IDS <- gsub("\'","",db.out$BEST_MATCH_IDS)
    db.out$BEST_MATCH_IDS <- gsub("#","",db.out$BEST_MATCH_IDS)
    
    db.out$species_name <- gsub("#","",db.out$species_name)
    db.out$genus_name <- gsub("#","",db.out$genus_name)
    db.out$SCIENTIFIC_NAME <- gsub("#","",db.out$SCIENTIFIC_NAME)
    
    db.out$species_name <- gsub("\'","",db.out$species_name)
    db.out$genus_name <- gsub("\'","",db.out$genus_name)
    db.out$SCIENTIFIC_NAME <- gsub("\'","",db.out$SCIENTIFIC_NAME)
    
    write.table(db.out,outfile,row.names=F,col.names=T,sep="\t",quote = FALSE)

    if (remove_DMS) {
      system(paste0("rm -r ",experiment, "_THOR_?? "),intern = T, wait = T)
    }
  
    message("THOR is done. He wrote the output file ",outfile," with ",length(info), " assigned sequences.")
    
  } else {
    if(run_ecotag) {
      message("run_ecotag is set to T. However, as vsearch has been selected, THOR will run vsearch instead of ecotag")
    }
    message("Running vsearch")
    divide_fasta(cores = cores, experiment = experiment)

    # check whether the taxonomic database exists
    tax_dir <- list.dirs(path = tax_db_dir, recursive = FALSE, full.names = TRUE)
    tax_dir <- tax_dir[grep(tax_db_name,tax_dir)]
    if (length(tax_dir) == 0) {
      message(paste0("The taxonomic db '", tax_db_name,"' does not exist in '", tax_db_dir,"'. Please check the path."))
      stop()
    } else if (length(grep("\\.fasta",tax_dir)) == 0) {
      message(paste0("The taxonomic db '", tax_dir,"' must be a fasta file. Please check the path."))
      stop()
    } else {
      tax_db <- tax_dir[grep("\\.fasta",tax_dir)]
    }

    outfile <- paste0(experiment,"_THOR_annotated.tsv")

    # Run vsearch
    if (cores>1) {
      for(i in 1:cores) {
        system(paste0("vsearch --usearch_global ",
                      experiment,"_ODIN_part_",sprintf("%02d",i),".fasta",
                      " --db ",tax_db," --id ", minimum_circle,
                      " --threads ",cores,
                      " --blast6out ",experiment,"_THOR_vsearch_part_",sprintf("%02d",i),".tsv"),
               intern = T, wait = T)
      }
    } else { 
      system(paste0("vsearch --usearch_global ",
                      experiment,"_ODIN_part_",sprintf("%02d",i),".fasta",
                      " --db ",tax_db," --id ", minimum_circle,
                      " --threads ",cores,
                      " --uc_allhits --maxaccept 3",
                      " --blast6out ",experiment,"_THOR_vsearch_part_01.tsv"),
               intern = T, wait = T)
    }
    # create a single file with all the results. This file is a tablewith the columns:
    # id, SCIENTIFIC_NAME, TAXONOMY, BEST_IDENTITY, BEST_MATCH_IDS, TAXID.
    # to create this file we need to read all the files created by vsearch and merge them.
    # then the file is sorted by id and the best hit is selected for each id. In case of
    # ties, we will create a consensus TAXONOMY which is the lowest common taxonomic rank.
    # from the output file of vsearch, the 9th column is the id, the 4th is the IDENTITY of each match,
    # and in the 10th column we have the taxonomy of each match. The taxonomy is a string with two parts separated
    # by a space. The first part is the MATCH_ID and the second part is the taxonomy. 
    
    # read all the files created by vsearch
    vsearch_assignment <- c()
    for (i in 1:cores) {
      vsearch_assignment <- rbind(vsearch_assignment,
                                  read.table(paste0(experiment,"_THOR_vsearch_part_",sprintf("%02d",i),".tsv"),
                                             sep="\t",header=F,fill = T))
    }
    
    # the function will be applied to each id using the apply function.
    # the output of the apply function will be a list of data frames.
    # the list will be converted to a data frame and written to a file.
    
    vsearch_assignment <- lapply(unique(vsearch_assignment$V9),FUN = create_taxonomy_assignment,vsearch_assignment = vsearch_assignment)
    vsearch_assignment <- do.call(rbind, vsearch_assignment)

    write.table(vsearch_assignment,outfile,row.names=F,col.names=T,sep="\t",quote = FALSE)
    message("THOR is done. He wrote the output file ",outfile," with ",sum(vsearch_assignment$TAXONOMY!="unidentified"), " assigned sequences.")
  }
}

create_taxonomy_assignment <- function(id,vsearch_assignment = NULL) {
  # for each id, select the best hit or create the consensus
  # create a function that will be applied to each id.
  # the function will select the best hit or create the consensus
  # if there is a tie in the best hit, the consensus will be created
  # by selecting the lowest common taxonomic rank.
  # the function will return a data frame with the columns:
  # id, SCIENTIFIC_NAME, TAXONOMY, BEST_IDENTITY, BEST_MATCH_IDS.
  # select the rows of the vsearch_assignment file that correspond to the id
  vsearch_assignment_id <- vsearch_assignment[vsearch_assignment$V9==id,]

  if (nrow(vsearch_assignment_id)==1) {
    if(vsearch_assignment_id$V1 == "N"){
      vsearch_assignment_id <- data.frame(id = id,
                                        SCIENTIFIC_NAME = "root",
                                        TAXONOMY = "unidentified",
                                        BEST_IDENTITY = 0,
                                        BEST_MATCH_IDS = "")
    } else {
      TAXONOMY <- strsplit(vsearch_assignment_id$V10," ")[[1]][2]
      SCIENTIFIC_NAME <- strsplit(TAXONOMY,';')[[1]]
      SCIENTIFIC_NAME <- SCIENTIFIC_NAME[length(SCIENTIFIC_NAME)]
      BEST_IDENTITY <- as.numeric(vsearch_assignment_id$V4)/100
      BEST_MATCH_IDS <- strsplit(vsearch_assignment_id$V10," ")[[1]][1]

      vsearch_assignment_id <- data.frame(id = id,
                                          SCIENTIFIC_NAME = SCIENTIFIC_NAME,
                                          TAXONOMY = TAXONOMY,
                                          BEST_IDENTITY = BEST_IDENTITY,
                                          BEST_MATCH_IDS = BEST_MATCH_IDS)
    }
  } else {
    # keep the highest MATCH_ID
    BEST_IDENTITY <- max(vsearch_assignment_id$V4)
    MATCH_TAXONOMY <- do.call(rbind,
                        strsplit(vsearch_assignment_id$V10[vsearch_assignment_id$V4==BEST_IDENTITY]," "))
    BEST_MATCH_IDS <- paste0(MATCH_TAXONOMY[,1],collapse="|")
    TAXONOMY <- Reduce(function(x, y) consensus_taxonomy(x, y), strsplit(MATCH_TAXONOMY[,2],";"))
    SCIENTIFIC_NAME <- TAXONOMY[length(TAXONOMY)]
    TAXONOMY <- paste0(TAXONOMY, collapse=";")
    
    vsearch_assignment_id <- data.frame(id = id,
                                        SCIENTIFIC_NAME = SCIENTIFIC_NAME,
                                        TAXONOMY = TAXONOMY,
                                        BEST_IDENTITY = as.numeric(BEST_IDENTITY)/100,
                                        BEST_MATCH_IDS = BEST_MATCH_IDS)
  }
  return(vsearch_assignment_id)
}

consensus_taxonomy <- function(vec1, vec2, n) {
  common_words <- character(0)
  for (i in seq_len(min(length(vec1), length(vec2)))) {
    if (vec1[i] == vec2[i]) {
      common_words <- c(common_words, vec1[i])
    } else {
      break
    }
  }
  common_words
}

