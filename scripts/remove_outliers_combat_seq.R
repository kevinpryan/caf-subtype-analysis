#!/usr/bin/Rscript

library(optparse)
library(biomaRt)
library(tximeta)
library(DESeq2)
library(PCAtools)
library(sva)
option_list = list(
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="metadata file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="specify outdir name with a forward slash at end", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 2){
  print_help(opt_parser)
  stop("You must specify your gtf file (-g) and your outfile (-o)", call.=FALSE)
}

metadata_file <- opt$metadata
outdir <- opt$out
# metadata file created with create_metadata.R
metadata <- read.table(metadata_file, sep = "\t")
files <- file.path(metadata$directory, rownames(metadata), "quant.sf")
coldata <- data.frame(files, names=rownames(metadata), Study = metadata$Study, 
                      Subpopulation = metadata$Subpopulation, 
                      Tumor_JuxtaTumor = metadata$Tumor_JuxtaTumor,
                      stringsAsFactors=FALSE)

# tx2gene but using the hgnc symbol instead of ensembl gene id version
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", #host="https://www.ensembl.org")
                host="uswest.ensembl.org")
tx2gene <- getBM(attributes = c("ensembl_transcript_id_version", "hgnc_symbol"), mart = mart, useCache = FALSE)
# read in data
se <- tximeta(coldata, skipMeta=TRUE, txOut=FALSE, tx2gene=tx2gene)

dds <- DESeqDataSet(se, design = ~1)
# returns a vector of whether the total count of each gene is >= 10 (True or false)
keep <- rowSums(counts(dds)) >= 10
# only keep rows (genes) for which keep is TRUE
dds <- dds[keep,]
# at least X samples with a count of 10 or more, where X is 5% of samples
X <- round(0.05*ncol(dds))
#X <- 7
keep <- rowSums(counts(dds) >= 10) >= X
dds <- dds[keep,]
vsd <- vst(dds, blind = TRUE)
rv <- rowVars(assay(vsd))
pc <- prcomp(t(assay(vsd)[head(order(-rv),nrow(assay(vsd))),]))
# -110 found by visual inspection - see README
idx <- pc$x[,1] < -110
print(paste("Number of samples being removed,", idx, sep = " "))
patient_samples <- rownames(metadata)
patient_samples <- patient_samples[!idx]
dds_remove_outliers <- dds[,!idx]
vsd_remove_outliers <- vst(dds_remove_outliers, blind = TRUE)
vsd_remove_outliers_mat <- assay(vsd_remove_outliers)
# remove batch effect
counts_matrix_remove_outliers <- assay(dds_remove_outliers)
batch_remove_outliers <- colData(dds_remove_outliers)$Study
print("batch correcting...")
ptm <- proc.time()
adjusted_remove_outliers <- ComBat_seq(counts = counts_matrix_remove_outliers, batch = batch_remove_outliers, group = colData(dds_remove_outliers)$Tumor_JuxtaTumor)
time_taken <- proc.time() - ptm
print("time taken...")
print(time_taken)
dds_batch_corrected <- DESeqDataSetFromMatrix(adjusted_all_samples, colData = metadata, design = ~1)
# not sure if blind should be TRUE or FALSE here
#vsd_batch_corrected <- vst(dds_batch_corrected, blind = TRUE)
date <- Sys.Date()
outfile <- paste(outdir, "dds_batch_corrected_", date, ".Rds", sep = "")
saveRDS(dds_batch_corrected, file = outfile)