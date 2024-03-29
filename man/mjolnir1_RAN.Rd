% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mjolnir1_RAN.R
\name{mjolnir1_RAN}
\alias{mjolnir1_RAN}
\title{RAN: Reads Allotment in N samples}
\usage{
mjolnir1_RAN(
  R1_filenames = "",
  lib_prefix = "",
  experiment = NULL,
  primer_F = "GGWACWRGWTGRACWNTNTAYCCYCC",
  primer_R = "TANACYTCNGGRTGNCCRAARAAYCA",
  cores = 1,
  R1_motif = "_R1",
  R2_motif = "_R2",
  ...
)
}
\arguments{
\item{R1_filenames}{Character vector with the names of the forward fastq or
fastq.gz files.}

\item{lib_prefix}{Character vector. Acronym for each sequencing library. This
acronym must be of 4 characters in capital letters. Do not mix up library and
experiment acronyms. The latter will be required in following steps. However
they can be the same.}

\item{experiment}{Character string. Acronym for the experiment. This
acronym must be of 4 characters in capital letters. Do not mix up library and
experiment acronyms. However they can be the same.}

\item{primer_F}{Character string of the Forward primer. Necessary when
samples are already demultiplexed.}

\item{primer_R}{Character string of the Reverse primer. Necessary when
samples are already demultiplexed.}

\item{cores}{Numeric. Number of parts into which the input files will be split
for the parallel processing of the FREYJA function.}

\item{R1_motif}{Character string that distinguishes the forward line file from
the reverse.}

\item{R2_motif}{Character string that distinguishes the reverse line file from
the forward.}
}
\description{
RAN function will prepare the FASTQ raw data for every sample without sample
tags and primers. The R1 output files will contain all forward sequences and
R2 the reverse sequences.
}
\details{
If samples are already demultiplexed primers need to be set or LERAY_XT primers
for COI will be used by default. Also RAN will read the
names of each individual R1 fastq files from a column in the LIBX_metadata.tsv
file, called fastq_name_R1. In the metadata table, each sample in the
original_samples column must have a their corresponding fastq_name_R1 and
mjolnir_agnomen (LIBX_sample_XXX).
If the fastq_name_R1. column is not in the metadata, RAN will use the lib_prefix.
Files required:
- ngsfilter file, needed only for multiplexed libraries
     For each library, a ngsfilter file is needed and must be named
     ngsfilter_<library identifier>.tsv. This must contain five tab-separated
     columns and no header. The first column with the library identifier, the
     second with the mjolnir_agnomens, the third with the sample tags, the
     fourth with the forward primers and the fifth with the reverse primers.
- metadata file
     a metadata file containing at least two columns required,
     'original_samples' and 'mjolnir_agnomens' (fastq_name_R1 for demultiplexed
      libraries), and named as <experiment identifier>_metadata.tsv.
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
}
