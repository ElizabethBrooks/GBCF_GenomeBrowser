#!/usr/bin/env Rscript

## Design gRNAs for CRISPRko with the SpCas9 nuclease
# End-to-end gRNA design workflow
# Tutorial to design gRNAs that knock out the human KRAS gene
# https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9

## Installation

# install the BiocManager and devtools
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#if (!requireNamespace("devtools", quietly = TRUE))
#  install.packages("devtools")

# Installing the core crisprVerse packages
#BiocManager::install(version="3.17", force = TRUE)
#BiocManager::install("crisprVerse")
#BiocManager::install("crisprBase")
#BiocManager::install("biomaRt")
#BiocManager::install("VariantAnnotation")
#BiocManager::install("crisprDesign")
#devtools::install_github("crisprVerse/crisprDesignData")
#devtools::install_github("crisprVerse/crisprViz")
#BiocManager::install("GenomicFeatures")

# set working directory
setwd("/Users/bamflappy/GBCF/genomeBrowser/daphnia")
#setwd("/home/ebrooks5/GenomeBrowser")

# start by loading the crisprVerse packages
library(crisprBase)
library(crisprDesign)
#library(crisprDesignData)
library(Rbowtie)
library(GenomicFeatures)
library(BSgenome)
library(BSgenome.Dpulex.KAP4)

# reference genome or chromosome fasta file path
#fastaFile <- "GCF_021134715.1_ASM2113471v1_genomic.fna"
fastaFile <- "GCF_021134715.1_ASM2113471v1_genomic_NC_060017.1.fna"

# prefix of the bowtie index
bowtie_index <- "NC_060017"

## build reference genome index using bowtie
## it is preferable to build in advance outside R
##bowtie_build(fastaFile,
##             outdir="./",
##             force=TRUE,
##             prefix=bowtie_index)

# retrieve input gff file path
inputGFF = "genomic.gff"

# store gff3 file as TxDb object
# NCBI gff files have transcript records without IDs that are set to NA
# it is possible to remove these records or subset the gff (e.g., by gene ID)
txdb_daphnia <- getTxDb(file=inputGFF, organism="Daphnia melanica")
#txdb_daphnia <- makeTxDbFromGFF(file=inputGFF, format="gff3", organism="Daphnia melanica")

# create a GRangesList object
# Organism 'Daphnia melanica' not recognized in biomaRt.
#grList <- TxDb2GRangesList(txdb_daphnia, seqlevelsStyle="NCBI")


## TO-DO ##
# path to VCF files for common SNPs (dbSNPs) downloaded from NCBI


## Nuclease specification

# load the nuclease object from the crisprBase package
data(SpCas9, package="crisprBase")
#data(CasRx, package="crisprBase")

# view the nuclease object
SpCas9
#CasRx

# inspect the protospacer construct 
prototypeSequence(SpCas9)
#prototypeSequence(CasRx)

## Specification of the target DNA sequence for a gene of interest (LOC124188748)
## product=deoxyribodipyrimidine photo-lyase-like
## transcript_id=XM_046581578.1

# obtain from crisprDesignData a GRangesList object that defines the 
# genomic coordinates of coding genes
# https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.html
txdb_tx <- transcriptsBy(txdb_daphnia, by="gene")
txdb_exons <- exonsBy(txdb_daphnia, by="gene")

# obtain a GRanges object containing the CDS coordinates of the gene of interest
gr_tx <- txdb_tx[["LOC124188748"]]
gr_exons <- txdb_exons[["LOC124188748"]]

## TO-DO ##
## consider removing NAs
## https://support.bioconductor.org/p/9135988/
##keep <- !is.na(mcols(txdb_transcripts)$tx_name)
##gr <- txdb_transcripts[keep,]

# only consider exons that constitute the primary transcript of the gene of interest
# transcript ID XM_046581578.1
primary_tx <- head(gr_tx, 1)$tx_name
gr <- gr_exons[grep(primary_tx, gr_exons$exon_name), ]


## Finding spacer sequences targeting KRAS

# obtain all possible spacer sequences that target protospacers located in our target DNA sequence
# Full genomic sequences for Homo sapiens as provided by UCSC (genome hg38, based on assembly
# GRCh38.p14 since 2023/01/31)
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40/
#bsgenome_daphnia <- Dpulex
guideSet <- findSpacers(gr,
                        bsgenome=Dpulex,
                        crisprNuclease=SpCas9)
                        #crisprNuclease=CasRx)

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
                                bsgenome=Dpulex,
                                ## annotate genomic alignments with a gene model annotation
                                ## https://bioconductor.org/packages/devel/bioc/vignettes/crisprDesign/inst/doc/intro.html
                                ## Organism 'Daphnia melanica' not recognized in biomaRt.
                                ##txObject=txdb_daphnia,
                                n_mismatches=2)

# view the results
guideSet

# inspect individual on- and off-targets and their context
alignments(guideSet)


## TO-DO ##
## Removing repeat elements

## remove promiscuous protospacer sequences that occur in repeats 
## or low-complexity DNA sequences
## a GRanges object containing repeat elements regions
##gr.repeats.daphnia <-
##guideSet <- removeRepeats(guideSet,
##                          gr.repeats=gr.repeats.daphnia)


## TO-DO ##
## Off-target scoring (MIT and CFD specificity scores)

## predict the likelihood of the nuclease to cut at the off-target 
## locations based on mismatch tolerance
## errors if no off targets to score
##guideSet <- addOffTargetScores(guideSet)

## view the results
##guideSet

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
