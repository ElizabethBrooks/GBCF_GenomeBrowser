#!/usr/bin/env Rscript

# R script to generate a TxDb object from a gff3 file
# Usage: Rscript generateTxDb_fromGFF3.r

# install Rsubread using BiocManager, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GenomicFeatures")

# load the Rsubread library for featureCounts
library(GenomicFeatures)
library(crisprDesign)
library(crisprDesignData)

# retrieve input gff file path
inputGFF = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/genomic.gff"

# store gff3 file as TxDb object
txdb <- makeTxDbFromGFF(file=inputGFF, organism="Daphnia melanica")

# create a GRangesList object
grList <- TxDb2GRangesList(txdb)

# build a tssObject
tssObject <- getTssObjectFromTxObject(grList)
