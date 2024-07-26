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
library(dplyr)
library(broom)
library(ggplot2)
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

metadata <- read.table(here("intermediate_files/metadata/metadata_all_samples.txt"), row.names = 1, sep = "\t")
dds.ensg.remove.outliers.batch.corrected <- readRDS("intermediate_files/dds_batch_corrected_group_tumor_ensembl_gene_id_version_2022-10-07.Rds")
dds.ensg.remove.outliers.batch.corrected.no.inhouse <- dds.ensg.remove.outliers.batch.corrected[,dds.ensg.remove.outliers.batch.corrected$Study != "InHouse"]
files <- file.path(metadata$directory, rownames(metadata), "quant.sf")
files.samples <- str_split_fixed(files, pattern = "/", n = 11)[,10]
files.outliers.removed <- files[files.samples %in% colnames(dds.ensg.remove.outliers.batch.corrected.no.inhouse)]

metadata.outliers.removed <- metadata[files.samples %in% colnames(dds.ensg.remove.outliers.batch.corrected.no.inhouse),]
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
dds.non.batch.corrected <- DESeqDataSet(se.ensg, design = ~1)


# read in: total read count from samtools.stats.alignment, % reads aligned?. qualimap - 5'-3', %exonic

files.samtools.stats <- paste(metadata.outliers.removed$directory, "/samtools_stats/", rownames(metadata.outliers.removed), ".markdup.sorted.bam.flagstat", sep = "")
files.qualimap <- paste(metadata$directory, "/qualimap/", rownames(metadata), "/rnaseq_qc_results.txt", sep = "")
samtools.stats <- sapply(files.samtools.stats, FUN = read.table, nrows = 1, sep = "+")
bias.5.3 <- sapply(files.qualimap, FUN = read.table, skip = 39, nrows = 1, sep = "=")
LinesAll <- sapply(files.qualimap, readLines)
subLAll <- sapply(LinesAll, FUN = grep, pattern = "5'-3' bias =") 
#inSub <- read.table(text=Lines[subL,,], header = FALSE)
#sapply(FUN = read.table, X = Lines, text = Lines[subL,,], header = FALSE )
samtools.stats.table <- data.frame(sample = rownames(metadata.outliers.removed), 
                                   nreads.samtools.stats = unlist(samtools.stats[seq(1, length(samtools.stats), 3)]), 
                                   Study = metadata.outliers.removed$Study, 
                                   tumour_juxta = metadata.outliers.removed$Tumor_JuxtaTumor,
                                   Subpopulation = metadata.outliers.removed$Subpopulation,
                                   Strandedness = metadata.outliers.removed$Strandedness)


#### carry out SVA ####
# r extract normalised counts and remove genes with low counts ENSG
dds.non.batch.corrected <- DESeq(dds.non.batch.corrected)
counts.ensg.remove.outliers.no.inhouse  <- counts(dds.non.batch.corrected, normalized = TRUE)
idx  <- rowMeans(counts.ensg.remove.outliers.no.inhouse) > 1
counts.ensg.remove.outliers.no.inhouse  <- counts.ensg.remove.outliers.no.inhouse[idx, ]
# SVA
mod <- model.matrix(~ Subpopulation + Tumor_JuxtaTumor, colData(dds.non.batch.corrected))
mod0 <- model.matrix(~ 1, colData(dds.non.batch.corrected))
svseq_no_inhouse <- svaseq(counts.ensg.remove.outliers.no.inhouse, mod, mod0, n.sv = 10) 
svseq_no_inhouse_svs <- svseq_no_inhouse$sv
cols <- paste("SV", seq(1,10), sep = "")
colnames(svseq_no_inhouse_svs) <- cols
sample_metadata_samtools_svs <- cbind.data.frame(samtools.stats.table, svseq_no_inhouse_svs)
sample_metadata_samtools_svs$Study <- as.factor(sample_metadata_samtools_svs$Study)
sample_metadata_samtools_svs$tumour_juxta <- as.factor(sample_metadata_samtools_svs$tumour_juxta)
sample_metadata_samtools_svs$Subpopulation <- as.factor(sample_metadata_samtools_svs$Subpopulation)
sample_metadata_samtools_svs$Strandedness <- as.factor(sample_metadata_samtools_svs$Strandedness)


glm.sv1 <- glm(sample_metadata_samtools_svs[,"SV1"] ~ sample_metadata_samtools_svs[,"nreads.samtools.stats"] +
                 sample_metadata_samtools_svs[,"Study"] +
                 sample_metadata_samtools_svs[,"tumour_juxta"] +
                 sample_metadata_samtools_svs[,"Subpopulation"]) 

glm.sv1 <- glm(SV1 ~ nreads.samtools.stats +
                 Study,
               data = sample_metadata_samtools_svs)
tidy(glm.sv1)

glm.tidy.output <- tibble()
for (i in 1:10){
  sv <- cols[i]
  glm.sv <- glm(get(sv) ~ nreads.samtools.stats +
                Study +
                tumour_juxta +
                Subpopulation, 
                data = sample_metadata_samtools_svs) 
                  
                
  glm.tidy <- tidy(glm.sv)
  glm.tidy$SV <- sv
  glm.tidy.output <- rbind.data.frame(glm.tidy.output, glm.tidy)
}
glm.tidy.output <- glm.tidy.output[which(!is.na(glm.tidy.output$p.value)),]

glm.tidy.output.sig <- glm.tidy.output[glm.tidy.output$p.value < 0.05 & !(glm.tidy.output$term == "(Intercept)"),]
glm.tidy.output.sig$term <- str_replace_all(glm.tidy.output.sig$term, "Study", "Study ")
View(glm.tidy.output.sig)
i <- 1
plots.out <- list()
for (i in 1:nrow(glm.tidy.output.sig)){
  if (grepl("Study", glm.tidy.output.sig$term[i]) == TRUE){
    plot.graph <- ggplot(sample_metadata_samtools_svs, aes(x=Study, y=glm.tidy.output.sig$SV[i])) +
      geom_jitter(position=position_jitter(0.2)) + 
      ggtitle(paste(glm.tidy.output.sig$term[i], "associated with", glm.tidy.output.sig$SV[i], sep = " ")) + 
      ylab(glm.tidy.output.sig$SV[i])
    plots.out[[i]] <- plot.graph
  } else{
    plot.graph <- ggplot(sample_metadata_samtools_svs, aes(x=tumour_juxta, y=SV1)) +
      geom_jitter(position=position_jitter(0.2)) + 
      ggtitle(paste(glm.tidy.output.sig$term[i], "associated with", glm.tidy.output.sig$SV[i], sep = " ")) +
      ylab(glm.tidy.output.sig$SV[i])
    plots.out[[i]] <- plot.graph
    
  }
}
plots.out[[1]]
plots.out[[2]]
plots.out[[3]]
plots.out[[4]]
plots.out[[5]]
plots.out[[6]]

ggplot(sample_metadata_samtools_svs, aes(x=Study, y=SV1)) +
  geom_jitter(position=position_jitter(0.2)) + 
  ggtitle("Study EGAD00001005744 associated with SV1")

ggplot(sample_metadata_samtools_svs, aes(x=Study, y=SV2)) +
  geom_jitter(position=position_jitter(0.2)) + 
  ggtitle("Study EGAD00001005744 and EGAD00001006144 associated with SV2")

ggplot(sample_metadata_samtools_svs, aes(x=Study, y=SV3)) +
  geom_jitter(position=position_jitter(0.2)) + 
  ggtitle("Study EGAD00001005744 and EGAD00001006144 associated with SV2")

ggplot(sample_metadata_samtools_svs, aes(x=Study, y=SV3)) +
  geom_jitter(position=position_jitter(0.2)) + 
  ggtitle("Study EGAD00001005744 and EGAD00001006144 associated with SV2")
#### write outputs ####
date <- Sys.Date()