#!/usr/bin/Rscript
## ---------------------------
##
## Script name: TPM singscore
##
## Purpose of script: Converts counts to TPM as recommended for using singscore, including filtering of genes with low counts
##
## Author: Kevin Ryan
##
## Date Created: 2022-11-04
##
## Email: k.ryan45@nuigalway.ie
##
## ---------------------------
##
## Notes:
##   
##The methods described in the Sigscore papers [@Foroutan2018],[@Bhuva2019] were used to apply the sigscore method. Steps to prepare data are:
#1. Read in count data
#2. Filter out genes with low counts
#3. Convert to TPM - here we use salmon effective gene lengths, which vary between samples
# This is carried out on both batch corrected data and non-batch corrected. Outlier samples are not included
# Here I am hard-coding in directories and file names, relative to "/home/kevin/Documents/PhD/subtypes/caf-subtype-analysis"
filter_out_low_expressed <- function(dds){
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

tpm.convert <- function(counts.mat, gene.length){
  x <- counts.mat / gene.length
  #Then with this matrix x, you do the following:
  tpm.mat <- t( t(x) * 1e6 / colSums(x) )
}

read_in_quant_genes <- function(file.name){
  x <- read.table(file = file.name, header = T)[,c(1,3)]
  #colnames(x)[2] <- as.character(sample.name)
  return(x)
}

getBM.call.distinct.ensg.version <- function(vec){
  library(biomaRt)
  library(dplyr)
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://www.ensembl.org")
  info <- getBM(attributes=c("hgnc_symbol",
                             "ensembl_gene_id_version",
                             "chromosome_name",
                             "start_position",
                             "end_position",
                             "strand",
                             "entrezgene_description",
                             "entrezgene_id"
  ),
  filters = "ensembl_gene_id_version",
  values = vec,
  mart = mart,
  useCache=FALSE)
  info.out <- distinct(info, ensembl_gene_id_version, .keep_all = TRUE)
}

## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(DESeq2)
library(tximeta)
library(here)
library(biomaRt)
library(SummarizedExperiment)
library(edgeR)
library(dplyr)
library(stringr)

# Define functions
filter_out_low_expressed <- function(dds){
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

tpm.convert <- function(counts.mat, gene.length){
x <- counts.mat / gene.length
#Then with this matrix x, you do the following:
tpm.mat <- t( t(x) * 1e6 / colSums(x) )
}

read_in_quant_genes <- function(file.name){
  x <- read.table(file = file.name, header = T)[,c(1,3)]
  return(x)
}

getBM.call.distinct.ensg.version <- function(vec){
  library(biomaRt)
  library(dplyr)
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://www.ensembl.org")
  info <- getBM(attributes=c("hgnc_symbol",
                               "ensembl_gene_id_version",
                               "chromosome_name",
                               "start_position",
                               "end_position",
                               "strand",
                               "entrezgene_description",
                               "entrezgene_id"
                               ),
                  filters = "ensembl_gene_id_version",
                  values = vec,
                  mart = mart,
                  useCache=FALSE)
  info.out <- distinct(info, ensembl_gene_id_version, .keep_all = TRUE)
}

##############################
## Non batch corrected data###
##############################

# Prepare to read in data

# metadata file created with create_metadata.R
metadata <- read.table(here("intermediate_files/metadata/metadata_all_samples.txt"), row.names = 1, sep = "\t")
dds.ensg.remove.outliers.batch.corrected <- readRDS("intermediate_files/dds_batch_corrected_group_tumor_ensembl_gene_id_version_2022-10-07.Rds")
files <- file.path(metadata$directory, rownames(metadata), "quant.sf")
files.samples <- str_split_fixed(files, pattern = "/", n = 11)[,10]
files.outliers.removed <- files[files.samples %in% colnames(dds.ensg.remove.outliers.batch.corrected)]
metadata.outliers.removed <- metadata[files.samples %in% colnames(dds.ensg.remove.outliers.batch.corrected),]
coldata <- data.frame(files.outliers.removed, names=rownames(metadata.outliers.removed), Study = metadata.outliers.removed$Study, 
                      Subpopulation = metadata.outliers.removed$Subpopulation, 
                      Tumor_JuxtaTumor = metadata.outliers.removed$Tumor_JuxtaTumor,
                      stringsAsFactors=FALSE)
colnames(coldata)[1] <- "files"

# tx2gene but using the hgnc symbol instead of ensembl gene id version
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://www.ensembl.org")
#host="uswest.ensembl.org")
tx2gene.ensg.version <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version"), mart = mart, useCache = FALSE)


# Read in data and filter out lowly expressed genes
# salmon was used in alignment mode so there is no salmon index, therefore there is no checksum to import the metadata 
# txOut = FALSE means to summarise to gene level (i.e. don't give out transcripts, give out gene level)
se.ensg <- tximeta(coldata, skipMeta=TRUE, txOut=FALSE, tx2gene=tx2gene.ensg.version)
annotation.data <- getBM.call.distinct.ensg.version(rownames(se.ensg))
stopifnot(all.equal(annotation.data$ensembl_gene_id_version, rownames(se.ensg)) == TRUE)
rowData(se.ensg) <- annotation.data
dge <- DGEList(counts = assay(se.ensg), genes = rowData(se.ensg))
# for each gene, get proportion of samples in which cpm is > 1
prop_expressed = rowMeans(cpm(dge) > 1)                                
# only keep genes in which this proportion s >0.5
keep = prop_expressed > 0.5  
dge <- dge[keep, , keep.lib.sizes = FALSE]
se.ensg <- se.ensg[keep,]
#rownames(se.ensg) <- gsub('\\.[0-9]*', '', rownames(se.ensg))

# getting effective lengths of genes
files.gene <- file.path(metadata.outliers.removed$directory, rownames(metadata.outliers.removed), "quant.genes.sf")
files.gene.samples <- str_split_fixed(files.gene, pattern = "/", n = 11)[,10]
files.gene.inhouse <- files.gene[files.gene.samples %in% colnames(se.ensg)]
genes.effective.lengths <- lapply(files.gene, read_in_quant_genes)

# combine all of these into a list of dataframes with the name of the sample
genes.effective.lengths.list <- list()
for (i in 1:length(genes.effective.lengths)){
  df <- genes.effective.lengths[[i]]
  colnames(df)[2] <- files.gene.samples[i]
  genes.effective.lengths.list[[i]] <- df
}
# get one big dataframe with the effective length of each gene in each sample
genes.effective.lengths.df <- genes.effective.lengths.list %>% purrr::reduce(inner_join, by = "Name")
tx2gene.ensg <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id_version"), mart = mart, useCache = FALSE)
colnames(tx2gene.ensg)[2] <- "Name"
# add hgnc symbol to gene length df
genes.effective.lengths.df.hgnc <- left_join(genes.effective.lengths.df, tx2gene.ensg,by = "Name")
df.se.ensg <- assay(se.ensg)
# filter out effective length genes that were already filtered out in terms of their CPM above
genes.effective.lengths.df.filtered <- genes.effective.lengths.df[which(genes.effective.lengths.df$Name %in% rownames(df.se.ensg)),]
# match the expression data matrix to the effective lengths
df.se.ensg.matched <- df.se.ensg[genes.effective.lengths.df.filtered$Name,]
rownames(genes.effective.lengths.df.filtered) <- genes.effective.lengths.df.filtered$Name
genes.effective.lengths.df.filtered$Name <- NULL
df.se.ensg.matched <- data.frame(df.se.ensg.matched)
genes.effective.lengths.df.filtered[1:5,1:5]
df.se.ensg.matched[1:5,1:5]
dim(genes.effective.lengths.df.filtered)
dim(df.se.ensg.matched)
stopifnot(dim(genes.effective.lengths.df.filtered) == dim(df.se.ensg.matched))
stopifnot(all.equal(rownames(genes.effective.lengths.df.filtered), rownames(df.se.ensg.matched)))
# convert to tpm
tpm.out <- tpm.convert(counts.mat =df.se.ensg.matched, gene.length = genes.effective.lengths.df.filtered )
colnames(tpm.out) <- gsub("X", "", colnames(tpm.out), )
tpm.out.remove.outliers.not.batch.corrected.genes.outfile <- paste("intermediate_files/singscore/caf_tpm_remove_outliers_engs_version_not_batch_corrected_remove_low_counts_for_singscore", Sys.Date(), ".txt", sep = "")
write.table(x = tpm.out, file = here(tpm.out.remove.outliers.not.batch.corrected.genes.outfile), sep = "\t", quote = F, row.names = T)

##############################
###  Batch-corrected data  ###
##############################


df.dds.ensg.remove.outliers.batch.corrected <- assay(dds.ensg.remove.outliers.batch.corrected)
dim(df.dds.ensg.remove.outliers.batch.corrected)

rownames(genes.effective.lengths.df) <- NULL
genes.effective.lengths.df.filtered <- genes.effective.lengths.df[which(genes.effective.lengths.df$Name %in% rownames(df.dds.ensg.remove.outliers.batch.corrected)),]
dim(genes.effective.lengths.df.filtered)
df.dds.ensg.remove.outliers.batch.corrected.matched <- df.dds.ensg.remove.outliers.batch.corrected[genes.effective.lengths.df.filtered$Name,]
rownames(genes.effective.lengths.df.filtered) <- genes.effective.lengths.df.filtered$Name
genes.effective.lengths.df.filtered$Name <- NULL
df.dds.ensg.remove.outliers.batch.corrected.matched <- data.frame(df.dds.ensg.remove.outliers.batch.corrected.matched)
genes.effective.lengths.df.filtered[1:5,1:5]
df.dds.ensg.remove.outliers.batch.corrected.matched[1:5,1:5]
tpm.out <- tpm.convert(counts.mat =df.dds.ensg.remove.outliers.batch.corrected.matched, gene.length = genes.effective.lengths.df.filtered )
tpm.out[1:5,1:5]
tpm.out.outfile <- paste("intermediate_files/singscore/caf_tpm_all_samples_outliers_removed_engs_version_batch_corrected_", Sys.Date(), ".txt", sep = "")
write.table(x = tpm.out, file = here(tpm.out.outfile), sep = "\t", quote = F, row.names = T)
q

