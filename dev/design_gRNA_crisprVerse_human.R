#!/usr/bin/env Rscript

## Design gRNAs for CRISPRko with the SpCas9 nuclease
# End-to-end gRNA design workflow
# Tutorial to design gRNAs that knock out the human KRAS gene
# https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9

# set working directory
setwd("/Users/bamflappy/GBCF/genomeBrowser/humanTest")
#setwd("/home/ebrooks5/GenomeBrowser")

## Installation

# install the BiocManager and devtools
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#if (!requireNamespace("devtools", quietly = TRUE))
#  install.packages("devtools")

# sudo R
#installed.packages()[, c("Package", "LibPath")]
#remove.packages(c("boot", "class", "cluster", "codetools", "foreign", "lattice", "MASS", "Matrix", "mgcv", "nlme", "nnet", "rpart", "spatial", "survival"), lib="/usr/lib/R/library")

# Installing the core crisprVerse packages
#BiocManager::install(version="3.17", force = TRUE)
#BiocManager::install("crisprVerse")
#BiocManager::install("crisprBase")
#BiocManager::install("biomaRt")
#BiocManager::install("crisprDesign")
#devtools::install_github("crisprVerse/crisprDesignData")
#devtools::install.packages("crisprVerse/crisprViz")
BiocManager::install("crisprBwa")

# Installing data packages
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor")

# start by loading the crisprVerse packages
library(crisprBase)
library(crisprDesign)
library(crisprDesignData)
library(Rbowtie)
library(crisprBwa)

# load the BSgenome package containing DNA sequences for the hg38 genome
library(BSgenome.Hsapiens.UCSC.hg38)

# reference genome file path
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
fastaFile <- "/Data/hg38.fa.gz"

# prefix of the hg38 bowtie index
bowtie_index <- "hg38"

# path to VCF files for common SNPs (dbSNPs) downloaded from NCBI
# https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz
#vcf <- "/Users/fortinj2/crisprIndices/snps/dbsnp151.grch38/00-common_all_snps_only.vcf.gz"

# build reference genome index using bowtie
#bowtie_build(fastaFile,
#             outdir="./",
#             force=TRUE,
#             prefix="hg38")

# build reference genome index using bwa
#bwa_build_index(fastaFile,
#                index_prefix="hg38")


## Nuclease specification

# load the SpCas9 nuclease object from the crisprBase package
data(SpCas9, package="crisprBase")

# view the SpCas9 nuclease object
SpCas9

# inspect the protospacer construct 
prototypeSequence(SpCas9)


## Specification of the target DNA sequence (KRAS CDS)

# obtain from crisprDesignData a GRangesList object that defines the 
# genomic coordinates (in hg38 coordinates) of coding genes in the human genome
data(txdb_human, package="crisprDesignData")

# obtain a GRanges object containing the CDS coordinates of KRAS
gr <- queryTxObject(txObject=txdb_human,
                    featureType="cds",
                    queryColumn="gene_symbol",
                    queryValue="KRAS")

# only consider exons that constitute the primary transcript of KRAS
# transcript ID ENST00000311936
gr <- gr[gr$tx_id == "ENST00000311936"]

# optionally, adjust the arguments in our call to queryTxObject to retrieve those transcript-specific coordinates
#gr <- queryTxObject(txObject=txObject,
#                    featureType="cds",
#                    queryColumn="tx_id",
#                    queryValue="ENST00000311936")


## Finding spacer sequences targeting KRAS

# obtain all possible spacer sequences that target protospacers located in our target DNA sequence
# Full genomic sequences for Homo sapiens as provided by UCSC (genome hg38, based on assembly
# GRCh38.p14 since 2023/01/31)
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40/
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
guideSet <- findSpacers(gr,
                        bsgenome=bsgenome,
                        crisprNuclease=SpCas9)

# view the genomic coordinates (PAM sites) for all found spacer sequences
guideSet

# extract information about the spacer sequences
spacers(guideSet)
protospacers(guideSet)
pams(guideSet)
head(pamSites(guideSet))
head(cutSites(guideSet))


## Characterizing gRNA spacer sequences

# evaluate the spacer sequences
guideSet <- addSequenceFeatures(guideSet)

# view the results
head(guideSet)


## Off-target search with bowtie

# identify potentially more problematic off-targets, such as those 
# located in the CDS of another gene
guideSet <- addSpacerAlignments(guideSet,
                                aligner="bowtie",
                                aligner_index=bowtie_index,
                                bsgenome=BSgenome.Hsapiens.UCSC.hg38,
                                n_mismatches=2,
                                txObject=txdb_human)

# view the results
guideSet

# inspect individual on- and off-targets and their context
alignments(guideSet)

## Removing repeat elements

# remove promiscuous protospacer sequences that occur in repeats 
# or low-complexity DNA sequences
data("gr.repeats.hg38", package="crisprDesignData")
guideSet <- removeRepeats(guideSet,
                          gr.repeats=gr.repeats.hg38)





## Off-target scoring (MIT and CFD specificity scores)

# predict the likelihood of the nuclease to cut at the off-target 
# locations based on mismatch tolerance
guideSet <- addOffTargetScores(guideSet)

# view the results
guideSet

## On-target scoring (gRNA efficiency)

# add the DeepHF and DeepSpCas9 scores
# requires package python=3.6
# create directory /Users/bamflappy/Library/Caches/org.R-project.R/R/ExperimentHub
guideSet <- addOnTargetScores(guideSet,
                              methods=c("deephf", "deepspcas9"))

# view the results
guideSet

## Restriction enzymes

# flag gRNAs containing restriction sites for a user-defined set of enzymes
guideSet <- addRestrictionEnzymes(guideSet)

# retrieve added annotations for the following commonly used enzymes
head(enzymeAnnotation(guideSet))

## Gene annotation

# adds transcript- and gene-level context to gRNAs
guideSet <- addGeneAnnotation(guideSet,
                              txObject=txdb_human)

# retrieve gene annotations
geneAnnotation(guideSet)

## TSS annotation

# find which protospacer sequences are located within promoter 
# regions of known genes
data(tssObjectExample, package="crisprDesign")
guideSet <- addTssAnnotation(guideSet,
                             tssObject=tssObjectExample)

# view the results
tssAnnotation(guideSet)

## SNP annotation

# annotate gRNAs with respect to a reference database of SNPs
guideSet <- addSNPAnnotation(guideSet, vcf=vcf)

# view the results
snps(guideSet)

## Filtering and ranking gRNAs

# filter out any unwanted gRNAs
# for example, only keep gRNAs that have percent GC between 20% and 80% 
# and that do not contain a polyT stretch
guideSet <- guideSet[guideSet$percentGC>=20]
guideSet <- guideSet[guideSet$percentGC<=80]
guideSet <- guideSet[!guideSet$polyT]

# rank gRNAs based on a set of criteria
# for example, sort gRNAs by the DeepHF on-target score
# Creating an ordering index based on the DeepHF score:
# Using the negative values to make sure higher scores are ranked first:
o <- order(-guideSet$score_deephf) 

# Ordering the GuideSet
guideSet <- guideSet[o]

# view the results
head(guideSet)

# sort gRNAs using several annotation columns
# for example, sort gRNAs using the DeepHF score, but also by prioritizing 
# first gRNAs that have no 1-mismatch off-targets in coding regions
o <- order(guideSet$n1_c, -guideSet$score_deephf) 

# Ordering the GuideSet
guideSet <- guideSet[o]

# view the results
head(guideSet)

# implement recommended rankings for the SpCas9, enAsCas12a and CasRx nucleases
# take into account the position of the gRNA within the target CDS of the 
# transcript ID in the ranking procedure
tx_id <- "ENST00000311936"
guideSet <- rankSpacers(guideSet,
                        tx_id=tx_id)

# view the results
head(guideSet)
