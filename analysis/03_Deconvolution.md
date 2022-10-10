Deconvolution file preparation
================

- <a href="#cibersort" id="toc-cibersort">CIBERSORT</a>
  - <a href="#generate-inputs-to-cibersort"
    id="toc-generate-inputs-to-cibersort">Generate inputs to CIBERSORT</a>

# CIBERSORT

CIBERSORT requires TPM as input, therefore we need to convert our batch
corrected counts to TPM. To get the gene lengths for this, we will use
the effective gene lengths from the Salmon output. Therefore, we must
use the ENGS gene id version.

## Generate inputs to CIBERSORT

The function `tpm.convert` is used to convert the counts to TPM.

``` r
tpm.convert <- function(counts.mat, gene.length){
x <- counts.mat / gene.length
#Then with this matrix x, you do the following:
tpm.mat <- t( t(x) * 1e6 / colSums(x) )
}

read_in_quant_genes <- function(file.name){
  x <- read.table(file = file.name, header = T)[,c(1,3)]
  #colnames(x)[2] <- as.character(sample.name)
  return(x)
}
```

``` r
#dds.ensg.remove.outliers.batch.corrected <- readRDS("../intermediate_files/dds_batch_corrected_group_tumor_ensembl_gene_id_version_2022-10-03.Rds")
dds.ensg.remove.outliers.batch.corrected <- readRDS("../intermediate_files/dds_batch_corrected_group_tumor_ensembl_gene_id_version_2022-10-07.Rds")
```

``` r
metadata <- read.table(here("intermediate_files/metadata/metadata_all_samples.txt"), row.names = 1, sep = "\t")
files.gene <- file.path(metadata$directory, rownames(metadata), "quant.genes.sf")
files.gene.samples <- str_split_fixed(files.gene, pattern = "/", n = 11)[,10]
files.gene.outliers.removed <- files.gene[files.gene.samples %in% colnames(dds.ensg.remove.outliers.batch.corrected)]
sample.names.outliers.removed <- files.gene.samples[files.gene.samples %in% colnames(dds.ensg.remove.outliers.batch.corrected)]
```

``` r
genes.effective.lengths <- lapply(files.gene.outliers.removed, read_in_quant_genes)
length(genes.effective.lengths)
```

    ## [1] 102

``` r
genes.effective.lengths.list <- list()
for (i in 1:length(genes.effective.lengths)){
  df <- genes.effective.lengths[[i]]
  colnames(df)[2] <- sample.names.outliers.removed[i]
  genes.effective.lengths.list[[i]] <- df
}
```

``` r
genes.effective.lengths.df <- genes.effective.lengths.list %>% purrr::reduce(inner_join, by = "Name")
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", #host="https://www.ensembl.org")
                host="useast.ensembl.org")
```

    ## Warning: Ensembl will soon enforce the use of https.
    ## Ensure the 'host' argument includes "https://"

``` r
tx2gene.ensg <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id_version"), mart = mart, useCache = FALSE)
colnames(tx2gene.ensg)[2] <- "Name"
genes.effective.lengths.df.hgnc <- left_join(genes.effective.lengths.df, tx2gene.ensg,by = "Name")
```

``` r
df.dds.ensg.remove.outliers.batch.corrected <- assay(dds.ensg.remove.outliers.batch.corrected)[2:nrow(assay(dds.ensg.remove.outliers.batch.corrected)),]
dim(df.dds.ensg.remove.outliers.batch.corrected)
```

    ## [1] 14251   102

``` r
rownames(genes.effective.lengths.df) <- NULL
genes.effective.lengths.df.filtered <- genes.effective.lengths.df[which(genes.effective.lengths.df$Name %in% rownames(df.dds.ensg.remove.outliers.batch.corrected)),]
dim(genes.effective.lengths.df.filtered)
```

    ## [1] 8824  103

``` r
df.dds.ensg.remove.outliers.batch.corrected.matched <- df.dds.ensg.remove.outliers.batch.corrected[genes.effective.lengths.df.filtered$Name,]
```

``` r
rownames(genes.effective.lengths.df.filtered) <- genes.effective.lengths.df.filtered$Name
genes.effective.lengths.df.filtered$Name <- NULL
df.dds.ensg.remove.outliers.batch.corrected.matched <- data.frame(df.dds.ensg.remove.outliers.batch.corrected.matched)
genes.effective.lengths.df.filtered[1:5,1:5]
```

    ##                      B86T3   B86T13   B86T16   B86T22 B103T103
    ## ENSG00000210194.1   16.772   16.823   16.720   16.922   15.630
    ## ENSG00000198695.2  276.634  251.026  246.299  232.295  227.839
    ## ENSG00000198886.2 1123.290 1092.570 1084.950 1064.020 1057.720
    ## ENSG00000212907.2  104.876   94.202   91.811   90.686   84.268
    ## ENSG00000198938.2  530.090  500.089  492.650  472.705  466.896

``` r
df.dds.ensg.remove.outliers.batch.corrected.matched[1:5,1:5]
```

    ##                    B86T3 B86T13 B86T16 B86T22 B103T103
    ## ENSG00000210194.1     31      7      8      3        3
    ## ENSG00000198695.2  15435   9371  10259   3881     4549
    ## ENSG00000198886.2 229714 117804 157657  98906    77018
    ## ENSG00000212907.2  15838   4750   5595   4565     3971
    ## ENSG00000198938.2 145999  59078  80861  65075    60425

``` r
tpm.out <- tpm.convert(counts.mat =df.dds.ensg.remove.outliers.batch.corrected.matched, gene.length = genes.effective.lengths.df.filtered )
tpm.out[1:5,1:5]
```

    ##                        B86T3      B86T13      B86T16      B86T22    B103T103
    ## ENSG00000210194.1   258.4799    63.58747    59.63943    26.86029    25.18767
    ## ENSG00000198695.2  7802.8100  5704.84946  5191.84982  2531.30774  2620.07417
    ## ENSG00000198886.2 28598.6497 16477.36348 18112.70336 14083.61073  9555.36400
    ## ENSG00000212907.2 21119.0404  7705.67065  7596.00276  7626.78847  6183.90652
    ## ENSG00000198938.2 38516.8024 18053.25178 20458.80963 20857.63899 16983.30548

``` r
tpm.out.mixture <- data.frame(tpm.out[,which(colData(dds.ensg.remove.outliers.batch.corrected)$Study == "InHouse")])
tpm.out.mixture.genes <- tibble::rownames_to_column(tpm.out.mixture, "genes")
tpm.out.mixture.genes[1:5,1:5]
```

    ##               genes       X4033       X4034       X4027      X4028
    ## 1 ENSG00000210194.1    78.86842    40.14785    37.68474   153.1497
    ## 2 ENSG00000198695.2  4156.83396  1777.45379  3651.67071  1928.5826
    ## 3 ENSG00000198886.2 16356.39095 11526.73650 10852.35787 11137.5262
    ## 4 ENSG00000212907.2  9048.97845  7373.91771  7037.73644  6651.0578
    ## 5 ENSG00000198938.2 23438.46149 18672.53618  7007.71542  9690.4570

``` r
tpm.out.mixture.genes.outfile <- paste("intermediate_files/cibersort/caf_tpm_mixture_engs_version_batch_corrected_", Sys.Date(), ".txt", sep = "")
#write.table(x = tpm.out.mixture.genes, file = here(tpm.out.mixture.genes.outfile), sep = "\t", quote = F, row.names = F)
```

``` r
# signature is our non-inhouse
tpm.out.signature <- tpm.out[,which(dds.ensg.remove.outliers.batch.corrected$Study != "InHouse")]
colnames(tpm.out.signature) <- dds.ensg.remove.outliers.batch.corrected$Subpopulation[which(dds.ensg.remove.outliers.batch.corrected$Study != "InHouse")]
tpm.out.signature[1:5,1:5]
```

    ##                           S3          S3          S3          S3          S3
    ## ENSG00000210194.1   258.4799    63.58747    59.63943    26.86029    25.18767
    ## ENSG00000198695.2  7802.8100  5704.84946  5191.84982  2531.30774  2620.07417
    ## ENSG00000198886.2 28598.6497 16477.36348 18112.70336 14083.61073  9555.36400
    ## ENSG00000212907.2 21119.0404  7705.67065  7596.00276  7626.78847  6183.90652
    ## ENSG00000198938.2 38516.8024 18053.25178 20458.80963 20857.63899 16983.30548

``` r
s1_cols <- which(colnames(tpm.out.signature) == "S1")
s3_cols <- which(colnames(tpm.out.signature) == "S3")
s4_cols <- which(colnames(tpm.out.signature) == "S4")

s1_cols_length <- length(s1_cols)
s3_cols_length <- length(s3_cols)
s4_cols_length <- length(s4_cols)
```

``` r
tpm.data.no.inhouse.order <- cbind.data.frame(rownames(tpm.out.signature),
                                              tpm.out.signature[,s1_cols],
                                              tpm.out.signature[,s3_cols],
                                              tpm.out.signature[,s4_cols]
                                              )
colnames(tpm.data.no.inhouse.order) <- NULL
rownames(tpm.data.no.inhouse.order) <- NULL
cols_no_inhouse <- c("genes", rep("S1", s1_cols_length), rep("S3", s3_cols_length), rep("S4", s4_cols_length))
tpm.data.no.inhouse.order <- rbind.data.frame(cols_no_inhouse, tpm.data.no.inhouse.order)
tpm.data.no.inhouse.order[1:5,1:5]
```

    ##                                                                       
    ## 1             genes               S1               S1               S1
    ## 2 ENSG00000210194.1 84.6588898736435 52.7477195924449  71.007487992008
    ## 3 ENSG00000198695.2 2968.94528474302 4123.22406999768 2782.02602409809
    ## 4 ENSG00000198886.2 7687.91744528988 12362.2845503602 10542.5483151727
    ## 5 ENSG00000212907.2 4022.41133006327 15549.2803836915 6641.02835167996
    ##                   
    ## 1               S1
    ## 2 22.9989250973712
    ## 3 1789.31114797989
    ## 4 5609.31823399334
    ## 5 3801.30789688231

``` r
tpm.data.no.inhouse.order.outfile <- paste("intermediate_files/cibersort/caf_tpm_for_signature_engs_version_batch_corrected_", Sys.Date(), ".txt", sep = "")
#write.table(tpm.data.no.inhouse.order, file = here(tpm.data.no.inhouse.order.outfile), sep = "\t", quote = F, row.names = F, col.names = F)
```

``` r
S1 <- c(rep(1,s1_cols_length), rep(2, (s3_cols_length+s4_cols_length)))
S3 <- c(rep(2,s1_cols_length), rep(1, s3_cols_length), rep(2, s4_cols_length))
S4 <- c(rep(2, (s1_cols_length+s3_cols_length)), rep(1, s4_cols_length))


if (length(S1) != length(colnames(tpm.out.signature))){
  stop("length s1 row incorrect")
} else if (length(S3) != length(colnames(tpm.out.signature))) {
  stop("length s3 row incorrect")
} else if (length(S4) != length(colnames(tpm.out.signature))) {
  stop("length s4 row incorrect")
}

out <- rbind.data.frame(S1, S3, S4)
colnames(out) <- NULL
row.names(out) <- c("S1", "S3", "S4")
outfile.name <-  paste("intermediate_files/cibersort/phenoclasses_caf_", Sys.Date(), ".txt", sep = "")
#write.table(out, file = here(outfile.name), sep = "\t", quote = F, col.names = F)
```
