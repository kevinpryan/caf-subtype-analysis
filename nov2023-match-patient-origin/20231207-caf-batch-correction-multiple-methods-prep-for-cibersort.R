library(tximeta)
library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(optparse)
option_list = list(
  make_option(c("-m", "--method"), type="character", default=NULL, 
              help="nextflow if running nf pipeline", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="counts_out.txt", 
              help="name of output counts matrix", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
params = parse_args(opt_parser);
if (params$method == ""){
  print_help(opt_parser)
  stop("must specify method at least", call.=FALSE)
} 

if (params$method == "nextflow") {
  #if (exists("params") == TRUE){
  metadata_path <- params$metadata
  print(paste("metadata file: ", metadata_path))
  tx2gene_file <- params$tx2gene
  print(paste("tx2gene file: ", tx2gene_file))
  #out <- params$out
  #print(paste("out file: ", out))
  #inhouse_metadata <- params$inhouse_metadata
  metadata_full <- read.csv(metadata_path)
  counts_matrix_in <- read.csv(params$counts_matrix, row.names = 1, header = T)
  out <- "batch_effect_removed.txt"
  #metadata_original <- read.table(params$metadata_without_patient)
  #} else if (method == "rstudio") {
} else {
  #metadata <- "/home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/nf-subpop/outdir/metadata/metadata_full.txt"
  metadata_path <- "~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/metadata_full_with_patient_20231127.csv"
  tx2gene_file <- "~/Documents/PhD/subtypes/caf-subtype-analysis/nf-subpop/outdir/tx2gene/tx2gene.txt"
  out <- "batch_effect_removed.txt"
  #inhouse_metadata <- "~/Documents/PhD/CAF_data/InHouse/reformat_samples.csv"
  metadata_full <- read.csv(metadata_path)
  counts_matrix_in <- read.csv("~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/20231127_counts_matrix_subpopulation.txt", row.names = 1, header = T)
  counts_matrix_inhouse <- read.csv("~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/20231207_counts_matrix_inhouse.txt")
  colnames(counts_matrix_inhouse)[1] <- "Gene"
  colnames(counts_matrix_inhouse) <- gsub(colnames(counts_matrix_inhouse), replacement = "", pattern = "X")
  sub <- "rstudio"
}

################################################
#  Approach:
#  
#  1. Use `limma::removebatcheffects` on variance-stabilised transformed data to remove batch effects
#  2. Get 2^ of this to unlog this
#  3. Do this for all 8 scenarios - combination of remove study batch effect, remove patient batch effect, remove outliers
#  4. Use this as input to CIBERSORT.
#  5. Use the 2^ method for the inhouse sample too - this will be our mixture
###################################################

source("~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/remove_batch_effects_prep_cibersort_functions.R")
binary_matrix <- expand.grid(lapply(1:2, function(x) c(0, 1)))
#colnames(binary_matrix) <- c("outliers", "study", "patient")
colnames(binary_matrix) <- c("study", "patient")

# want to read in subpop data
metadata <- read.csv(metadata_path)
metadata$directory <- gsub("kevin", sub, metadata$directory)
metadata_subpop <- metadata %>% dplyr::filter(Study != "InHouse")
metadata_inhouse <- metadata %>% dplyr::filter(Study == "InHouse")
files <- file.path(metadata_subpop$directory, 
                   metadata_subpop$Samplenames_old, 
                   "quant.sf") 
files_inhouse <- file.path(metadata_inhouse$directory, 
                           metadata_inhouse$Samplenames_old, 
                           "quant.sf") 
# make coldata given information for tximeta
coldata_subpop <- data.frame(files, 
                             names=metadata_subpop$Sample, 
                             Study = metadata_subpop$Study, 
                             Subpopulation = metadata_subpop$Subpopulation, 
                             Tumor_JuxtaTumor = metadata_subpop$Tumor_JuxtaTumor,
                             Patient = metadata_subpop$Patient,
                             stringsAsFactors=FALSE)

coldata_inhouse <- data.frame(files = files_inhouse, 
                             names=metadata_inhouse$Sample, 
                             Study = metadata_inhouse$Study, 
                             Subpopulation = metadata_inhouse$Subpopulation, 
                             Tumor_JuxtaTumor = metadata_inhouse$Tumor_JuxtaTumor,
                             Patient = metadata_inhouse$Patient,
                             stringsAsFactors=FALSE)

tx2gene <- read.table(tx2gene_file, header = T)

dds <- DESeqDataSetFromMatrix(countData=round(counts_matrix_in), colData=coldata_subpop, design=~ 1)
dds <- filter_out_low_expressed(dds)
vsd <- vst(dds, blind = TRUE)
vsd.mat <- assay(vsd)
rownames(metadata) <- metadata$Sample

# prep inhouse mixture data for cibersort
genes_counts_matrix_inhouse <- counts_matrix_inhouse$Gene
counts_matrix_inhouse <- round(counts_matrix_inhouse[,2:ncol(counts_matrix_inhouse)])
counts_matrix_inhouse <- cbind.data.frame(Genes = genes_counts_matrix_inhouse, counts_matrix_inhouse)
Tumor_JuxtaTumor <- factor(coldata_subpop$Tumor_JuxtaTumor, levels=c("juxtatumor","tumor"))
design <- model.matrix(~Tumor_JuxtaTumor)
#dds_inhouse <- DESeqDataSetFromMatrix(countData=counts_matrix_inhouse, colData=coldata_inhouse, design=~ 1, tidy = T)
dds_inhouse <- DESeqDataSetFromMatrix(countData=counts_matrix_inhouse, colData=coldata_inhouse, design=~ 1, tidy = T)
vsd_inhouse_non_log_for_cibersort <- inhouse_vsd_to_cibersort_format(dds_inhouse)
write.table(vsd_inhouse_non_log_for_cibersort, file = paste("~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/limma_batch_correction/", Sys.Date(), "_inhouse_data_for_cibersort.txt", sep = ""), quote = F, col.names = F, row.names = F, sep = "\t")

# remove batch effects with limma
# set up model matrix
Tumor_JuxtaTumor <- factor(coldata_subpop$Tumor_JuxtaTumor, levels=c("juxtatumor","tumor"))
design <- model.matrix(~Tumor_JuxtaTumor)

result_list <- lapply(1:nrow(binary_matrix), function(i) {
  row_vector <- as.vector(binary_matrix[i, ])
  limma_batch_effect(binary_vector = row_vector, 
                     vsd = vsd.mat, 
                     #model.matrix = matrix(1,nrow(design),1), 
                     model.matrix = design,
                     coldata = coldata_subpop)
})

for (i in 1:length(result_list)){
  write.table(result_list[[i]][[2]], file = result_list[[i]][[1]], sep = "\t", quote = F, row.names = F, col.names = F)
}

s1_cols <- which(as.vector(result_list[[1]][[2]][1,]) == "S1")
s3_cols <- which(as.vector(result_list[[1]][[2]][1,]) == "S3")
s4_cols <- which(as.vector(result_list[[1]][[2]][1,]) == "S4")

s1_cols_length <- length(s1_cols)
s3_cols_length <- length(s3_cols)
s4_cols_length <- length(s4_cols)

S1 <- c(rep(1,s1_cols_length), rep(2, (s3_cols_length+s4_cols_length)))
S3 <- c(rep(2,s1_cols_length), rep(1, s3_cols_length), rep(2, s4_cols_length))
S4 <- c(rep(2, (s1_cols_length+s3_cols_length)), rep(1, s4_cols_length))

if (length(S1) != length(colnames(vsd.mat))){
  stop("length s1 row incorrect")
} else if (length(S3) != length(colnames(vsd.mat))) {
  stop("length s3 row incorrect")
} else if (length(S4) != length(colnames(vsd.mat))) {
  stop("length s4 row incorrect")
}
outfile.name <- "~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/limma_batch_correction/phenoclasses_caf.txt"
out <- rbind.data.frame(S1, S3, S4)
colnames(out) <- NULL
row.names(out) <- c("S1", "S3", "S4")
write.table(out, file = outfile.name, sep = "\t", quote = F, col.names = F)

write.table(vsd_inhouse_mat_log, file = "~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/20231207_antilog_inhouse_vsd.txt", sep = "\t", quote = F, row.names = F, col.names = F)


