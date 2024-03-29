% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mjolnir2_FREYJA.R
\name{mjolnir2_FREYJA}
\alias{mjolnir2_FREYJA}
\title{FREYJA: Filtering of Reads, Enrollment, Yoke-reads Joining and Alignment}
\usage{
mjolnir2_FREYJA(
  experiment = NULL,
  cores = 1,
  Lmin = 299,
  Lmax = 320,
  score_obialign = 40,
  R1_motif = "_R1",
  R2_motif = "_R2",
  remove_DMS = T,
  run_on_tmp = F,
  ...
)
}
\arguments{
\item{experiment}{Character string. Acronym for the experiment. This
acronym must be of 4 characters in capital letters. Do not mix up library and
experiment acronyms. However they can be the same.}

\item{cores}{Numeric. Number of threads for parallel processing.}

\item{Lmin}{Numeric. Minimum bp length for a sequence to be accepted.}

\item{Lmax}{Numeric. Maximum bp length for a sequence to be accepted.}

\item{R1_motif}{Character string that distinguish the forward line file from
the reverse.}

\item{R2_motif}{Character string that distinguish the reverse line file from
the forward.}

\item{remove_DMS}{Logical. If TRUE, it will delete all obidms objects that are
created during the process. This can save a lot of hard disk space. The FALSE
option is useful for developing and debugging.}

\item{run_on_tmp}{Logical. If TRUE, the obidms objects will be created in
the /tmp location. This increases the speed as the communication within the
processor and the object that is being edited all the time is faster. However,
this method will consume much of the /tmp memory and it is recommended to have
three to four times the memory available in the /tmp directory than the original
forward files and remove_DMS=T}

\item{score_obilign}{Numeric. Minimum quality threshold to retain a sequence
in the qualiy filtering after pairalignment.}
}
\description{
FREYJA will use OBITools3 commands to merge paired-end reads, trim primer
sequences, filter by length, split sequences per sample and dereplicate
within each sample.
}
\details{
Input file fastq files are expected to be without primers sequence and all
forward sequences in the R1 file and all reverse sequences in the R2 file.
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
}
