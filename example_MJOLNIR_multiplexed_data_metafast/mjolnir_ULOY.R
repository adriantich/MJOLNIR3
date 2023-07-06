# Load MJOLNIR silently
suppressPackageStartupMessages(library(mjolnir))

# Define input fastq files (only names of R1 files are needed)
R1_filenames <-c("ULO1_R1.fastq.gz","ULO2_R1.fastq.gz","ULO3_R1.fastq.gz","ULO4_R1.fastq.gz")

# Input identifiers for the individual libraries to be used. It should be a 4-character name, matching the information in the ngsfilter files
lib_prefixes <- c("ULO1","ULO2","ULO3","ULO4")

# Input name for the final combined library (should be a 4-character name)
lib <- "ULOY"

mjolnir_PATH = "~/SOFT/MJOLNIR3/"

obipath="" 
vsearchpath = paste0(mjolnir_PATH,"vsearch/bin/")
swarmpath = paste0(mjolnir_PATH,"swarm/bin/")
dnoise_path = paste0(mjolnir_PATH,"DnoisE/src/")
tax_dir = paste0(mjolnir_PATH,"taxo_NCBI")
tax_dms_name = "DUFA_COI"

####################
# MJOLNIR pipeline #
####################

# Enter number of cores to be used in parallel. 
# Warning!:
# For RAN and FREYJA, until further optimisation, it is recommended to use between 5-7 cores
cores <- 7

# RAN will distribute R1 & R2 fastq files into equal-sized pieces for parallel computing
mjolnir1_RAN(R1_filenames = R1_filenames,cores = cores,lib_prefixes = lib_prefixes ,R1_motif="_R1",R2_motif="_R2")

# FREYJA will do the paired-end alignment, demultiplexing & length filtering. It will give individual filtered sample files as an output.
mjolnir2_FREYJA(lib_prefix = lib_prefixes,lib = lib,cores = cores,Lmin=299,Lmax=320,obipath = obipath,remove_DMS=F)

# If necessary, enter number of cores to be used in parallel. 
# You can enter the total number of cores available in the system, for full computing power.
# cores <- 16

# HELA will remove chimaeric sequences in a sample-by-sample basis, will change identifiers of remaining unique sequences & will generate a table of their abundances in each sample & a fasta file with unique sequences and their total abundance for ODIN
mjolnir3_HELA(lib = lib,cores = cores,vsearchpath = vsearchpath)

# ODIN will do the clustering & will generate a table with the abundances of each MOTU in each sample
mjolnir4_ODIN(lib = lib,cores = cores,d=13,min_reads_MOTU=2,min_reads_ESV=2,alpha=5, entropy=c(0.47,0.23,1.02,313),algorithm="DnoisE_SWARM",obipath=obipath, swarmpath=swarmpath, dnoise_path=dnoise_path, remove_singletons = TRUE,remove_DMS=F)

# THOR will asign the taxonomy to the representative sequence of each MOTU
mjolnir5_THOR(lib = lib,cores = cores, tax_dir=tax_dir, tax_dms_name=tax_dms_name, obipath=obipath, run_ecotag=T,remove_DMS=F)

# FRIGGA will integrate the information of MOTU abundances and taxonomy assignment from ODIN & THOR in a single table
mjolnir6_FRIGGA(lib = lib)

# LOKI kill remove the pseudogenes and will keep track of the taxonomic information of the removed MOTUs
mjolnir7_LOKI(lib = lib, min_id=.84,vsearchpath = vsearchpath)

# RAGNAROC will change the names of the samples to recover the original names and will remove unnecessary columns
mjolnir8_RAGNAROC(lib = lib,min_reads=2,min_relative=1/50000,
		remove_bacteria=T,remove_contamination=F,
                ESV_within_MOTU=T,blank_col="BLANK",blank_tag=T,remove_numts=T,cores=cores,blank_relative=0.1)
