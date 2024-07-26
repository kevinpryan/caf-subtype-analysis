#!/usr/bin/Rscript
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

library(dplyr)
library(DESeq2)
library(tximeta)

if (params$method == "nextflow") {
  # this is the full metadata file with the patient of origin
  metadata_path <- params$metadata
  print(paste("metadata file: ", metadata_path))
  tx2gene_file <- params$tx2gene
  print(paste("tx2gene file: ", tx2gene_file))
  out <- params$out
  #print(paste("out file: ", out))
  
  # this is the reformat_samples.csv file
  inhouse_metadata <- params$inhouse_metadata
  metadata_original <- read.csv(metadata_path)
  sub <- "kevin"
} else if (params$method == "rstudio") {
  metadata_path <- "~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/metadata_full_with_patient_20231127.csv"
  tx2gene_file <- "~/Documents/PhD/subtypes/caf-subtype-analysis/nf-subpop/outdir/tx2gene/tx2gene.txt"
  inhouse_metadata <- "~/Documents/PhD/CAF_data/InHouse/reformat_samples.csv"
  metadata_full <- read.csv(metadata_path)
  out <- params$out
  sub <- "rstudio"
} else {
  # not sure about this
  metadata_path <- "~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/metadata_full_with_patient_20231127.csv"
  tx2gene_file <- "~/Documents/PhD/subtypes/caf-subtype-analysis/nf-subpop/outdir/tx2gene/tx2gene.txt"
  inhouse_metadata <- "~/Documents/PhD/CAF_data/InHouse/reformat_samples.csv"
  metadata_full <- read.csv(metadata_path)
  out <- params$out
  sub <- "kevin"
}

# want to read in subpop data
metadata <- read.csv(metadata_path)
metadata$directory <- gsub("kevin", sub, metadata$directory)
metadata_subpop <- metadata %>% dplyr::filter(Study != "InHouse")
metadata_inhouse <- metadata %>% dplyr::filter(Study == "InHouse")

files <- file.path(metadata_subpop$directory, 
                   metadata_subpop$Samplenames_old, 
                   "quant.sf") 

file_inhouse <- file.path(metadata_inhouse$directory, 
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
coldata_inhouse <- data.frame(files = file_inhouse, 
                              names=metadata_inhouse$Sample, 
                              Study = metadata_inhouse$Study, 
                              Subpopulation = metadata_inhouse$Subpopulation, 
                              Tumor_JuxtaTumor = metadata_inhouse$Tumor_JuxtaTumor,
                              Patient = metadata_inhouse$Patient,
                              stringsAsFactors=FALSE)
tx2gene <- read.table(tx2gene_file, header = T)
# read in quant.sf file for each sample
se <- tximeta(coldata_subpop, skipMeta=TRUE, txOut=FALSE, tx2gene=tx2gene)
se_inhouse <- tximeta(coldata_inhouse, skipMeta=TRUE, txOut=FALSE, tx2gene=tx2gene)
dds <- DESeqDataSet(se, design = ~1)
dds_inhouse <- DESeqDataSet(se_inhouse, design = ~1)
counts_matrix <- counts(dds) 
counts_matrix_inhouse <- counts(dds_inhouse)
out <- "~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/20231127_counts_matrix_subpopulation.txt"
out_inhouse <- "~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/20231207_counts_matrix_inhouse.txt"
write.csv(counts_matrix, file = out, quote = FALSE)
write.csv(counts_matrix_inhouse, file = out_inhouse, quote = FALSE)
