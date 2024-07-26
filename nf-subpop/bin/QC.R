#!/usr/bin/Rscript
## ---------------------------
##
## Script name: SVA analysis non batch-corrected data
## Purpose of script: Testing out doing SVA on the non batch-corrected data. Looking at ways to test/visualise relationship between SVs and batch/etc.
##
## Author: Kevin Ryan
##
## Date Created: 2023-02-01
##
## Email: k.ryan45@nuigalway.ie
##
## ---------------------------
##
## Notes:
##   Putting this QC script into an Rscript so that it can go into a Nextflow pipeline.
##   Need to run from the command line, without hardcoding paths
##   This script carries out initial QC, including outlier removal and batch correction


## ---------------------------

#######################################
## load packages                    ###
#######################################
#library(here)
library(biomaRt)
library(tximeta)
library(DESeq2)
library(ggplot2)
library(PCAtools)
library(sva)
library(optparse)

#######################################
## Read in command-line arguments    ##
#######################################

option_list = list(
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="metadata file", metavar="character"),
  make_option(c("-t", "--tx2gene"), type="character", default=NULL, 
              help="tx2gene file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="specify outdir name with a forward slash at end", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

metadata <- opt$metadata
out <- opt$out
tx2gene_file <- opt$tx2gene
#if (length(opt) < 2){
 # print_help(opt_parser)
#  stop("You must specify your metadata file (-m) and your outfile (-o)", call.=FALSE)
#}


##########################
#### Define functions #### 
##########################
filter_out_low_expressed <- function(dds){
  library(DESeq2)
  print(paste("no of genes before filtering...", nrow(dds)))
  # returns a vector of whether the total count of each gene is >= 10 (True or false)
  keep <- rowSums(counts(dds)) >= 10
  # only keep rows (genes) for which keep is TRUE
  dds <- dds[keep,]
  # at least X samples with a count of 10 or more, where X is 5% of samples
  X <- round(0.05*ncol(dds))
  keep <- rowSums(counts(dds) >= 10) >= X
  dds <- dds[keep,]
  print(paste("no of genes after filtering...", nrow(dds)))
  return(dds)
}

svaBatchCor <- function(dat, mmi, mm0,n.sv=NULL, svs.include = NULL){
  library(sva)
  dat <- as.matrix(dat)
  Y <- t(dat)
  if(is.null(n.sv)) n.sv <- num.sv(dat,mmi,method="leek")
  o <- sva(dat,mmi,mm0,n.sv=n.sv)
  if (is.null(svs.include)){
    W <- o$sv
  } else{
    W <- o$sv[,svs.include]
  }
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  o$corrected <- t(Y - W %*% alpha)
  return(o)
}


######################
#### Read in data ####
######################

# metadata file created with create_metadata.R
metadata <- read.table(metadata, row.names = 1, sep = "\t")
# get files from metadata directory
files <- file.path(metadata$directory, rownames(metadata), "quant.sf")
# make coldata given information for tximeta
coldata <- data.frame(files, names=rownames(metadata), Study = metadata$Study, 
                      Subpopulation = metadata$Subpopulation, 
                      Tumor_JuxtaTumor = metadata$Tumor_JuxtaTumor,
                      stringsAsFactors=FALSE)
# tx2gene generated with generate_tx2gene_table.R or supplied
tx2gene <- read.table(tx2gene_file, header = T)
# read in quant.sf file for each sample
se <- tximeta(coldata, skipMeta=TRUE, txOut=FALSE, tx2gene=tx2gene)
dds <- DESeqDataSet(se, design = ~1)
# carry out variance-stabilising transformation
vsd <- vst(dds, blind = TRUE)

#########################################
#### Exploratory data analysis - PCA ####
#########################################


