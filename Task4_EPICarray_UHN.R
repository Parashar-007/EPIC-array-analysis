# Title: DNA Methylation  analysis on illumina's Methylation EPIC array data
# File: Task4_Methylation_UHN
# Project: UHN_Assignment

# Working Dir: setwd("C:/Users/surab/Desktop/UHN_Tasks/Methylation_data-20240221T021029Z-002")

#install and load the packages
BiocManager::install("minfi", force=TRUE)

#load libraries
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load EPIC array data
## Replace 'path_to_idat_files' with the directory containing your IDAT files

library(R.utils)
baseDir <- "Methylation_data/GSE153347_RAW/"

idat_files <- list.files(pattern = 'idat.gz', path = baseDir)
for (i in 1:length(idat_files)) {
  gunzip(filename = file.path(baseDir, idat_files[i]), 
         destname = gsub("[.]gz$", "", file.path(baseDir, idat_files[i])))
}

RGset <- readEPIC(idatPath = "baseDir")

# Pre-process data
# Background correction
RGset <- preprocessRaw(RGset)

# Subset to control probe data
control_probes <- getControlProbes(RGset)

# Noob background correction
RGset_noob <- preprocessNoob(RGset, control = control_probes)

# Quantile normalization
RGset_noob_quantile <- preprocessQuantile(RGset_noob)

# Quality control and filtering
Mset <- preprocessFilter(RGset_noob_quantile)

# DNA Methylation Analysis

# Design matrix
design <- model.matrix(~ factor(cases))

# Differential methylation analysis
fit <- lmFit(Mset, design)
fit <- eBayes(fit)

# Extract differentially methylated probes
results <- decideTests(fit)

# Visualization
volcanoplot(fit, coef = 2, highlight = 500, names = NULL)

