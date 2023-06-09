#!/usr/bin/Rscript
# Script used to generate tables for CIBERSORTx
library(optparse)
library(SummarizedExperiment)
library(tibble)
library(tximport)
library(here)

option_list = list(
  make_option(c("-m", "--metadata"), type="character", default="NULL", 
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
metadata <- read.table(metadata_file, row.names = 1)
idx_inhouse <- which(metadata$Study == "InHouse")
metadata_no_inhouse <- metadata[-c(idx_inhouse),]
files <- file.path(metadata$directory, rownames(metadata), "quant.sf")
coldata <- data.frame(files, names=rownames(metadata), Study = metadata$Study, 
                      Subpopulation = metadata$Subpopulation, 
                      Tumor_JuxtaTumor = metadata$Tumor_JuxtaTumor,
                      stringsAsFactors=FALSE)

# Get TPM expression for CIBERSORTx to make Signature matrix

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", #host="https://www.ensembl.org")
                host="uswest.ensembl.org")
tx2gene <- getBM(attributes = c("ensembl_transcript_id_version", "hgnc_symbol"), mart = mart, useCache = FALSE)

names = basename(dirname(coldata$files))
txi = tximport(coldata$files, type = "salmon", txOut = TRUE)
se <- SummarizedExperiment(assays = list(txi[["counts"]], abundance = txi[["abundance"]], length = txi[["length"]]), colData = coldata)
gi.s = summarizeToGene(txi, tx2gene = tx2gene,countsFromAbundance="scaledTPM")
gse.s = SummarizedExperiment(assays = list(counts = gi.s[["counts"]], abundance = gi.s[["abundance"]], length = gi.s[["length"]]),
                             colData = coldata)
abundance_tpm_data <- as.data.frame(assays(gse.s)[["abundance"]])
colnames(abundance_tpm_data) <- coldata$names
abundance_tpm_data <- abundance_tpm_data[-c(1),]
abundance_tpm_data_no_inhouse <- abundance_tpm_data[,colnames(abundance_tpm_data) %in% rownames(metadata_no_inhouse)]
abundance_tpm_data_inhouse <- abundance_tpm_data[,!(colnames(abundance_tpm_data) %in% rownames(metadata_no_inhouse))]
abundance_tpm_data_no_inhouse <- tibble::rownames_to_column(abundance_tpm_data_no_inhouse, "genes")
abundance_tpm_data_inhouse <- tibble::rownames_to_column(abundance_tpm_data_inhouse, "genes")
subpopulations_labelled <- metadata_no_inhouse$Subpopulation
colnames(abundance_tpm_data_no_inhouse)[-1] <- subpopulations_labelled
abundance_tpm_data_no_inhouse <- tibble(abundance_tpm_data_no_inhouse, .name_repair = "minimal")
s1_cols <- which(colnames(abundance_tpm_data_no_inhouse) == "S1")
s3_cols <- which(colnames(abundance_tpm_data_no_inhouse) == "S3")
s4_cols <- which(colnames(abundance_tpm_data_no_inhouse) == "S4")
abundance_tpm_data_no_inhouse_order <- tibble(genes = abundance_tpm_data_no_inhouse$genes, 
                                              abundance_tpm_data_no_inhouse[,s1_cols], 
                                              abundance_tpm_data_no_inhouse[,s3_cols],
                                              abundance_tpm_data_no_inhouse[,s4_cols],
                                              .name_repair = "minimal"
)
abundance_tpm_data_no_inhouse_order[1:5,1:5]
abundance_tpm_data_inhouse[1:5,1:5]
date <- Sys.Date()
outfile_for_sig_matrix <- paste(outdir, "caf_subtypes_tpm_for_sig_matrix_hgnc_", date, ".txt", sep = "")
print(paste("writing signature matrix to ", outfile_for_sig_matrix), sep = "")
write.table(abundance_tpm_data_no_inhouse_order, 
            file = outfile_for_sig_matrix, 
            sep = "\t", quote = F,
            row.names = F)

outfile_for_mixture_file <- paste(outdir, "caf_subtypes_tpm_mixture_hgnc_", date, ".txt", sep = "")
print(paste("writing mixture data to ", outfile_for_mixture_file), sep = "")

write.table(abundance_tpm_data_inhouse, 
            file = outfile_for_mixture_file, 
            sep = "\t", quote = F,
            row.names = F)

# create phenotype classes file
s1_cols_length <- length(s1_cols)
s3_cols_length <- length(s3_cols)
s4_cols_length <- length(s4_cols)
S1 <- c(rep(1,s1_cols_length), rep(2, (s3_cols_length+s4_cols_length)))
S3 <- c(rep(2,s1_cols_length), rep(1, s3_cols_length), rep(2, s4_cols_length))
S4 <- c(rep(2, (s1_cols_length+s3_cols_length)), rep(1, s4_cols_length))
if (length(S1) != length(colnames(abundance_tpm_data_no_inhouse_order))-1){
  stop("length s1 row incorrect")
} else if (length(S3) != length(colnames(abundance_tpm_data_no_inhouse_order))-1) {
  stop("length s3 row incorrect")
} else if (length(S4) != length(colnames(abundance_tpm_data_no_inhouse_order))-1) {
  stop("length s4 row incorrect")
}
out <- rbind.data.frame(S1, S3, S4)
colnames(out) <- NULL
row.names(out) <- c("S1", "S3", "S4")
out[,1:5]
outfile_for_phenotype_classes <- paste(outdir, "phenotype_classes_cibersort_", date, ".txt", sep = "")
print(paste("writing phenotype classes to ", outfile_for_phenotype_classes), sep = "")
write.table(out, file = outfile_for_phenotype_classes, sep = "\t", quote = F, col.names = F)
