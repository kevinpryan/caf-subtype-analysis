#!/usr/bin/Rscript
# Script to generate tx2gene file for given GTF file. Useful if tximeta does not have your version of the gtf file
library(optparse)
suppressPackageStartupMessages(library(GenomicFeatures))

option_list = list(
  make_option(c("-g", "--gtf"), type="character", default=NULL, 
              help="gencode gtf file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="specify outfile name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 2){
  print_help(opt_parser)
  stop("You must specify your gtf file (-g) and your outfile (-o)", call.=FALSE)
}


opt$gtf -> gtf_file
opt$out -> outfile
# read in gtf/gff file as TxDb object
txdb <- makeTxDbFromGFF(file=gtf_file)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
write.table(tx2gene, file = outfile, quote = F, row.names = F, sep = "\t")