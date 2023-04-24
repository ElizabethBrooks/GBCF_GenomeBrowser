
## Design gRNAs for CRISPRko with the SpCas9 nuclease
## End-to-end gRNA design workflow
## Hsapiens Tutorial

## Installation

# Installing the core crisprVerse packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version="devel")
BiocManager::install("crisprVerse")

# Installing data packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version="devel")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor")

# start by loading the crisprVerse packages
library(crisprBase)
library(crisprDesign)
library(crisprDesignData)

# load the BSgenome package containing DNA sequences for the hg38 genome
library(BSgenome.Hsapiens.UCSC.hg38)

## Nuclease specification

# load the SpCas9 nuclease object from the crisprBase package
data(SpCas9, package="crisprBase")

# view the SpCas9 nuclease object
SpCas9


