#!/usr/bin/Rscript

library(optparse)
library(biomaRt)
library(tximeta)
library(DESeq2)
library(PCAtools)
library(sva)
#library(here)
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
  stop("You must specify your metadata file (-m) and your outfile (-o)", call.=FALSE)
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
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://www.ensembl.org")
                #host="uswest.ensembl.org")
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
print(paste("Number of samples being removed,", length(idx[which(idx == TRUE)]), sep = " "))
patient_samples <- rownames(metadata)
patient_samples <- patient_samples[!idx]
# create new metadata file 
metadata_reduced <- metadata[patient_samples,]
date <- Sys.Date()
metadata_outfile <- paste(outdir, "metadata_outliers_removed_", date, ".txt", sep = "")
print(paste("writing reduced metadata file to ", metadata_outfile, sep = ""))
write.table(metadata_reduced, file = metadata_outfile, sep = "\t", quote = F, row.names = T)
dds_remove_outliers <- dds[,!idx]
vsd_remove_outliers <- vst(dds_remove_outliers, blind = TRUE)
vsd_remove_outliers_mat <- assay(vsd_remove_outliers)
# remove batch effect
counts_matrix_remove_outliers <- assay(dds_remove_outliers)
batch_remove_outliers <- colData(dds_remove_outliers)$Study
print(paste("Number of genes before filtering: ", nrow(counts_matrix_remove_outliers), sep = ""))
# remove genes with 0 in > 1/3 of samples as per GitHub user benostendorf https://github.com/zhangyuqing/ComBat-seq/issues/20 
counts_filt <- counts_matrix_remove_outliers[apply(counts_matrix_remove_outliers, 1, function(x) sum(x == 0)) < ncol(counts_matrix_remove_outliers) / 3, ]
print(paste("Number of genes after filtering: ", nrow(counts_filt), sep = ""))
combat_seq_group <- colData(dds_remove_outliers)$Tumor_JuxtaTumor
print("batch correcting...")
ptm <- proc.time()
adjusted_remove_outliers <- ComBat_seq(counts = counts_filt, batch = batch_remove_outliers, group = combat_seq_group, full_mod = TRUE)
time_taken <- proc.time() - ptm
print("time taken...")
print(time_taken)
dds_batch_corrected <- DESeqDataSetFromMatrix(adjusted_remove_outliers, colData = colData(dds_remove_outliers), design = ~1)
# not sure if blind should be TRUE or FALSE here
#vsd_batch_corrected <- vst(dds_batch_corrected, blind = TRUE)
outfile <- paste(outdir, "dds_batch_corrected_group_tumor_", date, ".Rds", sep = "")
print(paste("writing batch corrected data to ", outfile, sep = ""))
saveRDS(dds_batch_corrected, file = outfile)