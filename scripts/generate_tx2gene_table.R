#!/usr/bin/Rscript
# Script to generate tx2gene file for given GTF file. Useful if tximeta does not have your version of the gtf file
# if you want the hgnc_symbol then use the --hgnc parameter
library(optparse)
suppressPackageStartupMessages(library(GenomicFeatures))

option_list = list(
  make_option(c("-g", "--gtf"), type="character", default=NULL, 
              help="gencode gtf file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="specify outfile name", metavar="character"),
  make_option(c("-x", "--hgnc"), type ="logical", default = FALSE,
              help = "supply TRUE if you want the gene names to be in hgnc format")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 2){
  print_help(opt_parser)
  stop("You must specify your gtf file (-g) and your outfile (-o)", call.=FALSE)
}


opt$gtf -> gtf_file
opt$out -> outfile
opt$hgnc -> hgnc
# read in gtf/gff file as TxDb object
txdb <- makeTxDbFromGFF(file=gtf_file)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
if (hgnc == TRUE){
   suppressPackageStartupMessages(library(tidyr))
   suppressPackageStartupMessages(library(rtracklayer))
   suppressPackageStartupMessages(library(GenomicRanges))
   suppressPackageStartupMessages(library(dplyr))
   gff_in <- import.gff(gtf_file)
   gff_transcript_ensg_geneid <- data.frame(TXNAME = gff_in$transcript_id, GENEID = gff_in$gene_id, GENENAME = gff_in$gene_name)
   gff_transcript_ensg_geneid_notdup_transcript <- gff_transcript_ensg_geneid[!duplicated(gff_transcript_ensg_geneid$TXNAME),]
   tx2gene <- gff_transcript_ensg_geneid_notdup_transcript %>% drop_na() %>% dplyr::select(-GENEID)
}
write.table(tx2gene, file = outfile, quote = F, row.names = F, sep = "\t")
