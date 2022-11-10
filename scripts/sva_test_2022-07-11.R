#!/usr/bin/Rscript
## ---------------------------
##
## Script name: SVA analysis non batch-corrected data
## Purpose of script: Testing out doing SVA on the non batch-corrected data. Looking at ways to test/visualise relationship between SVs and batch/etc.
##
## Author: Kevin Ryan
##
## Date Created: 2022-11-07
##
## Email: k.ryan45@nuigalway.ie
##
## ---------------------------
##
## Notes:
##   
##


## ---------------------------

#######################################
## load packages                    ###
#######################################
library(DESeq2)
library(sva)
library(tximeta)
library(biomaRt)
library(here)
library(stringr)
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

# metadata file created with 02-QC.Rmd
metadata <- read.table(here("intermediate_files/metadata/metadata_no_inhouse.txt"), row.names = 1, sep = "\t")
files <- file.path(metadata$directory, rownames(metadata), "quant.sf")
coldata <- data.frame(files, names=rownames(metadata), Study = metadata$Study, 
                      Subpopulation = metadata$Subpopulation, 
                      Tumor_JuxtaTumor = metadata$Tumor_JuxtaTumor,
                      stringsAsFactors=FALSE)
# tx2gene but using the hgnc symbol instead of ensembl gene id version
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://www.ensembl.org")
#host="uswest.ensembl.org")
tx2gene.ensg.version <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version"), mart = mart, useCache = FALSE)
se.ensg.non.batch.corrected <- tximeta(coldata, skipMeta=TRUE, txOut=FALSE, tx2gene=tx2gene.ensg.version)
dds.non.batch.corrected <- DESeqDataSet(se.ensg.non.batch.corrected, design = ~1)
dds.non.batch.corrected

# read in: total read count from samtools.stats.alignment, % reads aligned?. qualimap - 5'-3', %exonic
files.samtools.stats <- paste(metadata$directory, "/samtools_stats/", rownames(metadata), ".markdup.sorted.bam.flagstat", sep = "")
files.qualimap <- paste(metadata$directory, "/qualimap/", rownames(metadata), "/rnaseq_qc_results.txt", sep = "")
samtools.stats <- sapply(files.samtools.stats, FUN = read.table, nrows = 1, sep = "+")
bias.5.3 <- sapply(files.qualimap, FUN = read.table, skip = 39, nrows = 1, sep = "=")
LinesAll <- sapply(files.qualimap, readLines)
subLAll <- sapply(LinesAll, FUN = grep, pattern = "5'-3' bias =") 
inSub <- read.table(text=Lines[subL,,], header = FALSE)
#sapply(FUN = read.table, X = Lines, text = Lines[subL,,], header = FALSE )
samtools.stats.table <- data.frame(sample = rownames(metadata), nreads.samtools.stats = unlist(samtools.stats[seq(1, length(samtools.stats), 3)]))

#### write outputs ####
date <- Sys.Date()