#!/usr/bin/Rscript

library(optparse)
library(DESeq2)
library(PCAtools)

option_list = list(
  make_option(c("-d", "--dds"), type="character", default=NULL, 
              help="dds file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="specify outdir name with a forward slash at end", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 2){
  print_help(opt_parser)
  stop("You must specify your metadata file (-m) and your outfile (-o)", call.=FALSE)
}
dds_path <- opt$dds
outdir <- opt$out
# get last element of list
#x[[-1]]
# readRDS(file = "intermediate_files/dds_batch_corrected_group_tumor_2022-09-15.Rds")

dds <- readRDS(dds_path)
vsd <- vst(dds, blind = FALSE)
# Biplot
biplot(p, showLoadings = T, lab = NULL)

# PCA plot colour by Subpopulation
plotPCA(vsd, intgroup = "Subpopulation", ntop = nrow(vsd))
# PCA plot colour by Study
plotPCA(vsd, intgroup = "Study", ntop = nrow(vsd))

peigencor <- eigencorplot(p,
                          components = getComponents(p, 1:10),
                          metavars = colnames(colData(dds)),
                          col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
                          cexCorval = 0.7,
                          colCorval = 'black',
                          fontCorval = 2,
                          posLab = 'bottomleft',
                          rotLabX = 45,
                          posColKey = 'top',
                          cexLabColKey = 1.5,
                          scale = TRUE,
                          corFUN = 'pearson',
                          corUSE = 'pairwise.complete.obs',
                          corMultipleTestCorrection = 'none',
                          main = 'PCs clinical correlations',
                          colFrame = 'white',
                          plotRsquared = TRUE)