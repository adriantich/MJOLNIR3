#' THOR: Taxonomy with Higher-than-Order Ranks
#'
#' This is a wrapper of ecotag
#'
#' @details 
#' After assignment with ecotag, higher taxa at ranks higher than order are added from custom CSV files.
#' 
#' The database used can be download or build using the NJORDR package (see https://github.com/adriantich/NJORDR-MJOLNIR3)
#' 
#' @param lib Character string. Acronym for the experiment. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. However they can be the same.
#' 
#' @param cores Numeric. Number of threads for parallel processing.
#' 
#' @param tax_dir String specifying de PATH to the directory were the taxonomic 
#' information is stored
#' 
#' @param tax_dms_name Character string specifying the name of the obidms object
#' without the ".obidms" extension. 
#' 
#' @param obipath Character string specifying the PATH to the obi binary.
#' 
#' @param run_ecotag Logical. Whether to run (TRUE, default) the ecotag taxonomic
#' assignment or not (FALSE). The latter could take place when alternative taxonomic
#' assignament software is applied but adding higher taxonomic ranks is desired.
#' 
#' @param remove_DMS Logical. If TRUE, it will delete all obidms objects that are
#' created during the process. This can save a lot of hard disk space. The FALSE 
#' option is useful for developing and debugging.
#' 
#' @param minimum_circle Numeric. from ecotag: Minimum identity considered for 
#' the assignment circle (sequence is assigned to the LCA of all sequences 
#' within a similarity circle of the best matches; the threshold for this circle 
#' is the highest value between <CIRCLE_THRESHOLD> and the best assignment score 
#' found for the query sequence). Give value as a normalized identity, e.g. 0.95 
#' for an identity of 95%.
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

mjolnir5_THOR <- function(lib,cores,tax_dir,tax_dms_name=NULL,obipath="",run_ecotag=T,remove_DMS=T, minimum_circle=0.7){

  old_path <- Sys.getenv("PATH")
  if (is.null(obipath)) obipath <- "~/obi3-env/bin/"
  obipath <- path.expand(obipath)
  Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))

  if (run_ecotag) {
    message("THOR will assign the taxonomy to the order level with ecotag.")

    if (is.null(tax_dms_name)) {
      message("obidms with taxonomic database is required")
      exit()
      # system(paste0("obi import ",tax_dir,"/",ref_db," ",lib, "_THOR_taxo/ref_seqs ; ",
      #               "obi import --taxdump ",tax_dir,"/",tax_db," ",lib, "_THOR_taxo/taxonomy/my_tax ; ",
      #               "obi grep --require-rank=species --require-rank=genus --require-rank=family --taxonomy ",lib, "_THOR_taxo/taxonomy/my_tax ",lib, "_THOR_taxo/ref_seqs ",lib, "_THOR_taxo/ref_seqs_clean ; ",
      #               "obi uniq --taxonomy ",lib, "_THOR_taxo/taxonomy/my_tax ",lib, "_THOR_taxo/ref_seqs_clean ",lib, "_THOR_taxo/ref_seqs_uniq ; ",
      #               "obi build_ref_db -t 0.95 --taxonomy ",lib, "_THOR_taxo/taxonomy/my_tax ",lib, "_THOR_taxo/ref_seqs_uniq ",lib, "_THOR_taxo/ref_db "),
      #        intern = T, wait = T)
      # tax_dms_name <- paste0(lib, "_THOR_taxo")
    }

    # it is necessary to run ecotag within a new directory for each part.
    # this is because dms can not be calles from two processes at the same time
    # and can not change the name of the dms so make a copy in each directory and
    # run there the ecotag
    X <- NULL
    for (i in 1:cores) {
      system(paste0("mkdir ",lib, "_THOR_",sprintf("%02d",i)," ; cp -r ",tax_dir,"/",tax_dms_name,".obidms ",lib, "_THOR_",sprintf("%02d",i),"/. ; "),intern = T, wait = T)
      X <- c(X,paste0("cd ",lib, "_THOR_",sprintf("%02d",i)," ; ",
                      "obi import --fasta-input ../",lib,"_ODIN_part_",sprintf("%02d",i),".fasta ",tax_dms_name,"/seqs ; ",
                      # "obi import --fasta-input ",lib,"_ODIN_part_",sprintf("%02d",i),".fasta ",lib, "_THOR_",sprintf("%02d",i),"/seqs ; ",
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
                      "obi export --fasta-output ",tax_dms_name,"/assigned_seqs_add >../",lib,"_THOR_part_",sprintf("%02d",i),".fasta"))
      # "obi ecotag --taxonomy ",lib, "_THOR_",sprintf("%02d",i),"/taxonomy/my_tax -R ",lib, "_THOR_",sprintf("%02d",i),"/ref_db ",lib, "_THOR_",sprintf("%02d",i),"/seqs ",lib, "_THOR_",sprintf("%02d",i),"/assigned_seqs"))
    }

    suppressPackageStartupMessages(library(parallel))
    no_cores <- cores
    clust <- makeCluster(no_cores)
    clusterExport(clust, list("X","old_path","obipath"),envir = environment())
    clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))})
    parLapply(clust,X, function(x) system(x,intern=T,wait=T))
    stopCluster(clust)
  }

  message("THOR will add higher taxonomic ranks now.")
  filefasta <-paste0(lib,"_THOR.fasta")
  system(paste0("cat ",lib,"_THOR_part_??.fasta | grep '>' > ",filefasta),intern=T,wait=T)
  outfile <- paste0(lib,"_THOR_annotated.tsv")
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
  showlines <- 1

  family_to_order <- read.table(paste0(tax_dir,"/","family_to_order.csv"),sep=",",head=T)
  genus_to_family <- read.table(paste0(tax_dir,"/","genus_to_family.csv"),sep=",",head=T)

  # get_rank <- function(cadena){
  #   if (gregexpr("rank=",cadena)[[1]]>0)
  #   {cadena <- substr(cadena,gregexpr("rank=",cadena)[[1]]+5,nchar(cadena))
  #   cadena <- substr(cadena,1,regexpr(";",cadena)[[1]]-1)
  #   }  else cadena <- NA
  #   return(cadena)
  # }

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
    # for (fila in 1:length(info)) {

    infofasta <- info[fila]
    # message(infofasta)
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
      # matrix.data[match(c("superkingdom_name","kingdom_name","phylum_name","class_name"),rownames(matrix.data)),2] <- taxo_names[taxo_names$order_name==matrix.data["order_name",2],c("superkingdom_name","kingdom_name","phylum_name","class_name")]
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
          # matrix.data[c("superkingdom_name","kingdom_name","phylum_name"),2] <- taxo_names[taxo_names$class_name==matrix.data["class_name",2],c("superkingdom_name","kingdom_name","phylum_name")]
        }
      } else if (("kingdom_name" %in% na_taxa) & !("phylum_name" %in% na_taxa)) {
        if (matrix.data["phylum_name",2] %in% taxo_names$phylum_name) {
          matrix.data[match(c("superkingdom_name"),rownames(matrix.data)),2] <- taxo_names[taxo_names$phylum_name==matrix.data["phylum_name",2],c("superkingdom_name")][1]
          matrix.data[match(c("kingdom_name"),rownames(matrix.data)),2] <- taxo_names[taxo_names$phylum_name==matrix.data["phylum_name",2],c("kingdom_name")][1]
          # matrix.data[c("superkingdom_name","kingdom_name"),2] <- taxo_names[taxo_names$phylum_name==matrix.data["phylum_name",2],c("superkingdom_name","kingdom_name")][1,]
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

  db.out <- do.call("rbind",mclapply(1:length(info), complete_taxa,
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
  # write.table(reordered_db.out,outfile,row.names=F,col.names=T,sep="\t",quote = FALSE)

  if (remove_DMS) {
    system(paste0("rm -r ",lib, "_THOR_?? "),intern = T, wait = T)
  }
  message("THOR is done. He wrote the output file ",outfile," with ",length(info), " assigned sequences.")
}

