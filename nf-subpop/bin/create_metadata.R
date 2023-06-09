#!/usr/bin/Rscript

library(optparse)
library(stringr)
option_list = list(
  make_option(c("-f", "--out_full"), type="character", default="metadata_full.txt", 
              help="specify output file for all samples", metavar="character"),
  make_option(c("-r", "--out_reduced"), type="character", default="metadata_no_inhouse.txt",
              help="specify outfile name for reduced dataset without inhouse samples", metavar="character"),
  make_option(c("-b", "--basedir"), type="character", default="/home/kevin/Documents/PhD/CAF_data/",
              help="basedir with all datasets i.e. /home/kevin/Documents/PhD/CAF_data/", metavar="character")
); 

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

#if (length(opt) < 2){
#  print_help(opt_parser)
#  stop("You must specify your output file for all samples (-f) and the reduced dataset without inhouse samples (-r)", call.=FALSE)
#}

opt$out_full -> out_full
opt$out_reduced -> out_reduced
opt$basedir -> basedir
file_4810 <- paste(basedir, "/EGAD00001004810/delimited_maps/Run_Sample_meta_info.map", sep = "")
EGAD_4810 <- read.table(file_4810, 
                        sep = ";")
EGAD_4810_cancers <- str_split_fixed(EGAD_4810$V1, pattern = "=", n = 2)[,2]
EGAD_4810_keep <- which(EGAD_4810_cancers == "BC")
EGAD_4810_remove <- which(EGAD_4810_cancers != "BC")

print("removing samples: ")
print(EGAD_4810[EGAD_4810_remove,][,4])
EGAD_4810_removed <- EGAD_4810[EGAD_4810_remove,]

EGAD_4810_filtered <- EGAD_4810[EGAD_4810_keep,]
EGAD_4810_meta <- data.frame(
  Sample = str_split_fixed(EGAD_4810_filtered$V4, pattern = "=", n = 2)[,2],
  Study = "EGAD00001004810",
  Subpopulation = "S3",
  Tumor_JuxtaTumor = tolower(str_split_fixed(EGAD_4810_filtered$V3, pattern = " ", n = 2)[,2]),
  Strandedness = "unstranded",
  directory = paste(basedir, "nfcore_results/EGAD00001004810_nfcore_results/star_salmon", sep = ""),
  row.names = 1
)

file_3808 <- paste(basedir, "/EGAD00001003808/meta_CAF-S1_S4_BC_47samples.txt", sep = "")
EGAD_3808 <- read.table(file_3808, 
                        header = T, sep = "\t")
EGAD_3808_meta <- data.frame(
  Sample = EGAD_3808$Sample.Name,
  Study = "EGAD00001003808",
  Subpopulation = EGAD_3808$subset,
  Tumor_JuxtaTumor = EGAD_3808$Type,
  Strandedness = "unstranded",
  directory =   paste(basedir, "nfcore_results/EGAD00001003808_nfcore_results/star_salmon", sep = ""),
  row.names = 1
)

file_6144 <- paste(basedir, "/EGAD00001006144/meta_7samples.txt", sep = "")
EGAD_6144 <- read.table(file_6144, 
                        header = T,sep = "\t")
EGAD_6144_meta <- data.frame(
  Sample = paste("CAF_Culture_", EGAD_6144$Sample.Name, sep = ""),
  Study = "EGAD00001006144",
  Subpopulation = "S1",
  Tumor_JuxtaTumor = "tumor",
  Strandedness = "reverse",
  directory =   paste(basedir, "nfcore_results/EGAD00001006144_nfcore_results/star_salmon", sep = ""),
  row.names = 1
)

file_5744 <- paste(basedir, "/EGAD00001005744/metaData_Pelon_et_al.txt", sep = "")
EGAD_5744 <- read.table(file_5744, header =T, check.names = F)
EGAD_5744$Sample.Name <- gsub(pattern = "\\.", replacement = "-", x = EGAD_5744$Sample.Name )
EGAD_5744_filtered <- EGAD_5744[!(EGAD_5744$subset == "EPCAM+") & (EGAD_5744$Type == "T"),]
EGAD_5744_to_remove <- EGAD_5744[(EGAD_5744$subset == "EPCAM+") | (EGAD_5744$Type == "LN"),]

EGAD_5744_meta <- data.frame(
  Sample = EGAD_5744_filtered$Sample.Name,
  Study = "EGAD00001005744",
  Subpopulation = EGAD_5744_filtered$subset,
  Tumor_JuxtaTumor = "tumor",
  Strandedness = "unstranded",
  directory = paste(basedir, "nfcore_results/EGAD00001005744_nfcore_results/star_salmon", sep = ""),
  row.names = 1
)

file_barkley <- paste(basedir, "/InHouse/reformat_samples.csv", sep= "")
barkley_samples <- read.csv(file_barkley, header = T, row.names = "samples", check.names = F)
barkley_samples_meta <- data.frame(
  Sample = row.names(barkley_samples),
  Study = "InHouse",
  Subpopulation = "Unknown",
  Tumor_JuxtaTumor = ifelse(barkley_samples$Condition == "Tumour", "tumor", "juxtatumor"),
  Strandedness = "unstranded",
  directory =   paste(basedir, "nfcore_results/inhouse_data_nfcore_results_version_3_8_1/star_salmon", sep = ""),
  row.names = 1
)
metadata <- rbind.data.frame(EGAD_4810_meta, EGAD_3808_meta, 
                             EGAD_6144_meta, EGAD_5744_meta, barkley_samples_meta)
metadata$Tumor_JuxtaTumor <- gsub(x = metadata$Tumor_JuxtaTumor, pattern = "-", replacement = "")

samples_not_included_in_analysis <- c(EGAD_5744_to_remove$Sample.Name, str_split_fixed(EGAD_4810_removed$V4, n = 2, pattern = "=")[,2]) 
#write.table(samples_not_included_in_analysis, file = "/home/kevin/Documents/PhD/CAF_data/samples_not_included_in_analysis_2022-10-14.txt", quote = F, row.names = F, col.names = F)
paste(samples_not_included_in_analysis, collapse = "|")
print("directories...")
print(metadata$directory)
write.table(metadata, file = out_full, quote = F, sep = "\t", row.names = T)

metadata_no_inhouse <- rbind.data.frame(EGAD_4810_meta, EGAD_3808_meta, 
                                        EGAD_6144_meta, EGAD_5744_meta)
metadata_no_inhouse$Tumor_JuxtaTumor <- gsub(x = metadata_no_inhouse$Tumor_JuxtaTumor, pattern = "-", replacement = "")

metadata_no_inhouse_no6144 <- rbind.data.frame(EGAD_4810_meta, EGAD_3808_meta, EGAD_5744_meta)
name_without_6144 <- paste(str_split_fixed(out_reduced, pattern = "\\.", n = 2)[,1], "_without_6144.", str_split_fixed(out_reduced, pattern = "\\.", n = 2)[,2], sep = "")
write.table(metadata_no_inhouse, file = out_reduced, quote = F, sep = "\t", row.names = T)
write.table(metadata_no_inhouse_no6144, file = name_without_6144, quote = F, sep = "\t", row.names = T)
