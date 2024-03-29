% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mjolnir8_RAGNAROC.R
\name{mjolnir8_RAGNAROC}
\alias{mjolnir8_RAGNAROC}
\title{RAGNAROC: Replace AGnomens with Names And Recover Original Codification}
\usage{
mjolnir8_RAGNAROC(
  experiment = NULL,
  metadata_table = "",
  output_file = "",
  output_file_ESV = "",
  min_reads = 0,
  remove_bacteria = T,
  remove_contamination = F,
  contamination_file = "contaminants.txt",
  ESV_within_MOTU = T,
  remove_numts = F,
  cores = 1,
  ...
)
}
\arguments{
\item{experiment}{Character string. Acronym for the experiment. This
acronym must be of 4 characters in capital letters. Do not mix up library and
experiment acronyms. However they can be the same.}

\item{metadata_table}{tsv table. if not specified, the file must be named 
<EXPX>_metadata.tsv. This table must have: a column named "mjolnir_agnomens"
with the names given to the samples during the pipeline in FREYJA; a column 
named "original_samples" with the samples names that will be given to the 
samples at the end of the pipeline; and a column with the name specified in
in the "blank_col" parameter ("BLANK" by default) where blanks, negatives and
controls are tagged with a flag specified in the "blank_tag" parameter (T by
default).}

\item{output_file}{Character string specifying the outputfile name}

\item{output_file_ESV}{Character string specifying the outputfile name 
for ESVs abundances if required.}

\item{min_reads}{Number of the minimum number of reads allowed for each MOTU/ESV
or ESV within MOTU.}

\item{remove_bacteria}{Logical. If TRUE it will apply the bacteria removal 
filtering (see Details).}

\item{remove_contamination}{Logical. If TRUE it will apply the contamination removal 
filtering (see Details).}

\item{contamination_file}{Character string specifying the name of the contamination file. (see Details)}

\item{ESV_within_MOTU}{Logical. If TRUE this will take into account the ESV 
that were clustered into MOTUs in ODIN if algorithm was set to "DnoisE_SWARM" 
or "SWARM_DnoisE" and apply all filters to both data.}

\item{remove_numts}{Logical whether to apply the NUMT filter (TRUE) or not (FALSE)}

\item{cores}{Numeric. Number of threads for parallel processing during NUMT 
removal.}
}
\description{
Final step of the MJOLNIR3 pipeline to apply the last filters to the abundance data
}
\details{
RAGNAROC consists on different contamination removals and filtering steps 
as follows:

Removal of Bacteria: this removed the Units tagged as "Prokaryota" or "root" in the <EXPX>_LOKI_Cutated.tsv

Removal of contaminations: this step removes the taxa specified in the "contaminaion_file"

NUMT removal: this step is design for the Leray-XT COI marker. It deletes all
sequences that do not have a 313 (plus/minus multiple of 3 equivalent to a codon)
bp length. Then removes sequences with stop codons and those metazoan 
sequences that do not translate for 5 specific conservative aminoacids.
}
\examples{
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
mjolnir3_HELA(experiment = experiment, cores = cores)

# Run ODIN
mjolnir4_ODIN(experiment = experiment, cores = cores, d = 13, 
              min_reads_MOTU = 2, min_reads_ESV = 2,
              min_relative = 1 / 50000, blank_relative = 0.1, 
              metadata_table = "", blank_col = "BLANK", blank_tag = TRUE, 
              alpha = 4, entropy = c(0.47, 0.23, 1.02, 313), 
              algorithm = "DnoisE_SWARM")

# set the directory where the database is stored
tax_db <- "~/taxo_NCBI/DUFA_COI"

# Run THOR
mjolnir5_THOR(experiment = experiment, cores = cores, 
              tax_db = tax_db, run_ecotag = T)

# Run FRIGGA
mjolnir6_FRIGGA(experiment = experiment)

# Run LOKI
mjolnir7_LOKI(experiment = experiment, min_id = .84)

# Run RAGNAROC
mjolnir8_RAGNAROC(experiment = experiment, remove_bacteria = T, ESV_within_MOTU = T,
                  remove_numts = F, cores = 1)
}
