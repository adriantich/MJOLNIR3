library(mjolnir)

# Define input fastq files (only names of R1 files are needed)
R1_filenames <- c("ULO1_R1.fastq.gz", "ULO2_R1.fastq.gz", "ULO3_R1.fastq.gz",
                  "ULO4_R1.fastq.gz")

# Input identifiers for the individual libraries to be used.
# It should be a 4-character name, matching the information in the
# ngsfilter files.
lib_prefixes <- c("ULO1", "ULO2", "ULO3", "ULO4")

# experiment identifier
experiment <- 'ULOY'
# Enter number of cores to be used in parallel.
cores <- 7

mjolnir1_RAN(R1_filenames, lib_prefix = lib_prefixes, experiment = experiment,
             cores = cores, R1_motif = "_R1", R2_motif = "_R2")

# Run FREYJA
mjolnir2_FREYJA(experiment = experiment, cores = cores, Lmin=299, Lmax=320)

# Run HELA
mjolnir3_HELA(experiment, cores)

# Run ODIN
mjolnir4_ODIN(experiment, cores, d = 13,
              min_reads_MOTU = 2, min_reads_ESV = 2,
              min_relative = 1 / 50000, blank_relative = 0.1,
              metadata_table = "", blank_col = "BLANK", blank_tag = TRUE,
              alpha = 4, entropy = c(0.47, 0.23, 1.02, 313),
              algorithm = "DnoisE_SWARM")

# set the directory where the database is stored
tax_dir <- "~/taxo_NCBI/"
tax_dms_name <- "DUFA_COI"

# Run THOR
mjolnir5_THOR(experiment, cores,
              tax_dir = tax_dir, tax_dms_name = tax_dms_name,run_ecotag = T)

# Run FRIGGA
mjolnir6_FRIGGA(experiment)

# Run LOKI
mjolnir7_LOKI(experiment, min_id=.84)