#!/usr/bin/Rscript
# Script used to generate tables for CIBERSORTx - moving out of Rmd as it is distracting

library(SummarizedExperiment)
library(tibble)
library(tximport)
metadata <- read.table("/home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/metadata_caf_subtypes.txt", row.names = 1)
idx_6144 <- which(metadata$Study == "EGAD00001006144")
metadata <- metadata[-c(idx_6144),]
idx_inhouse <- which(metadata$Study == "InHouse")
metadata_no_inhouse <- metadata[-c(idx_inhouse),]
files <- file.path(metadata$directory, rownames(metadata), "quant.sf")
coldata <- data.frame(files, names=rownames(metadata), Study = metadata$Study, 
                      Subtype = metadata$Subtype, 
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
abundance_tpm_data_no_inhouse <- abundance_tpm_data[,colnames(abundance_tpm_data) %in% rownames(metadata_no_inhouse)]
abundance_tpm_data_inhouse <- abundance_tpm_data[,!(colnames(abundance_tpm_data) %in% rownames(metadata_no_inhouse))]
abundance_tpm_data_no_inhouse <- tibble::rownames_to_column(abundance_tpm_data_no_inhouse, "genes")
abundance_tpm_data_inhouse <- tibble::rownames_to_column(abundance_tpm_data_inhouse, "genes")
subtypes_labelled <- metadata_no_inhouse$Subtype
colnames(abundance_tpm_data_no_inhouse)[-1] <- subtypes_labelled
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
#write.table(abundance_tpm_data_no_inhouse_order, 
 #          file = "/home/kevin/Documents/PhD/cibersort/caf_subtypes_tpm_for_sig_matrix_no6144_hgnc.txt", 
  #        sep = "\t", quote = F,
   #      row.names = F)
#write.table(abundance_tpm_data_inhouse, 
    #       file = "/home/kevin/Documents/PhD/cibersort/caf_tpm_mixture_hgnc.txt", 
     #     sep = "\t", quote = F,
      #   row.names = F)

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
#write.table(out, file = "/home/kevin/Documents/PhD/cibersort/phenoclasses_caf_no6144.txt", sep = "\t", quote = F, col.names = F)
