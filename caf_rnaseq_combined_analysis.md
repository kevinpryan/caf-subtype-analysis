CAF subtype analysis
================
Kevin Ryan
2022-06-29 22:31:46

-   <a href="#introduction" id="toc-introduction">Introduction</a>
    -   <a href="#preparation" id="toc-preparation">Preparation</a>
        -   <a href="#create-sample-file" id="toc-create-sample-file">Create Sample
            File</a>
        -   <a href="#prepare-transcript-annotations"
            id="toc-prepare-transcript-annotations">Prepare transcript
            annotations</a>
        -   <a href="#read-in-files" id="toc-read-in-files">Read in files</a>
        -   <a href="#pca" id="toc-pca">PCA</a>
    -   <a href="#references" id="toc-references">References</a>

# Introduction

Cancer-associated fibroblasts (CAFs) are a heterogeneous cell type found
in the tumour microenvironment. They have a wide array of functions, and
tend to be immunosuppressive and cancer-promoting. There have been many
attempts to characterise subtypes of CAFs, with much transcriptomic
analysis being carried out in the Mechta-Grigoriou lab in Institut
Curie. They have identified 4 ‘subtypes’ which can be separated based on
the expression of different markers:

-   S1: FAP<sup>High</sup>, CD29<sup>Med-High</sup>,
    α<sup>SMAHigh</sup>, PDPN<sup>High</sup>, PDGFRβ<sup>High</sup>
-   S2: FAP<sup>Neg</sup>, CD29<sup>Low</sup>, αSMANeg-<sup>Low</sup>,
    PDPN<sup>Low</sup>, PDGFRβ<sup>Low</sup>
-   S3: FAP<sup>Neg-Low</sup>, CD29<sup>Med</sup>,
    αSMA<sup>Neg-Low</sup>, PDPN<sup>Low</sup>, PDGFRβ<sup>Low-Med</sup>
-   S4: FAP<sup>Low-Med</sup>, CD29<sup>High</sup>, αSMA<sup>High</sup>,
    PDPN<sup>Low</sup>, PDGFRβ<sup>Med</sup>

(Pelon et al. 2020)

FACS gating strategies can be used to isolate these various subtypes.
The Mechta-Grigoriou group have done this and have generated bulk
RNA-sequencing data for the S1, S3 and S4 subtypes. They generated
scRNA-sequencing data for the S1 subtype. This data was deposited on the
European Genome Phenome Archive, and was accessed via a Data Transfer
Agreement.

The following summarises the data obtained:

<table>
<colgroup>
<col style="width: 18%" />
<col style="width: 18%" />
<col style="width: 27%" />
<col style="width: 36%" />
</colgroup>
<thead>
<tr class="header">
<th>Subtype</th>
<th>Total samples</th>
<th>Studies (Samples)</th>
<th>Notes</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>S1</td>
<td>28</td>
<td><ul>
<li>EGAD00001003808 (16)</li>
<li>EGAD00001005744 (5)</li>
<li>EGAD00001006144 (7)</li>
</ul></td>
<td><ul>
<li>3808 has 12xJuxta-tumor</li>
<li>5744 5 samples from LN</li>
<li>Sorting vs spreading</li>
</ul></td>
</tr>
<tr class="even">
<td>S2</td>
<td>0</td>
<td>N/A</td>
<td>N/A</td>
</tr>
<tr class="odd">
<td>S3</td>
<td>14</td>
<td><ul>
<li>EGAD00001004810 (14)</li>
</ul></td>
<td><ul>
<li>4810 has 11xJuxta-tumor</li>
<li>Ovarian</li>
</ul></td>
</tr>
<tr class="even">
<td>S4</td>
<td>15</td>
<td><ul>
<li>EGAD00001003808 (10)</li>
<li>EGAD00001005744 (5)</li>
</ul></td>
<td><ul>
<li>3808 has 9xJuxta-tumor</li>
<li>5744 5 samples from LN</li>
</ul></td>
</tr>
</tbody>
</table>

With the juxta-tumour data, they got tumour and juxta-tumour data from
the same patient. However, I have not been able to figure out whether
they came from the same patient. Could probably use Optitype to
determine HLA allele - match tumour and juxta tumour.

We also have scRNA-seq data for S1.

The data was processed using nf-core/rnaseq version using the default
parameters.

We would expect our tumour-associated normal to be most like the S3
subtype (usually accumulate in juxta-tumours).

## Preparation

### Create Sample File

Columns will be: Sample, Study, CAF_subtype, Tumor_Juxtatumor

*Here we will be combining data from 5 studies. To begin with, we will
only include the metadata available for all studies (except for our
unknown CAF subtype label). Breast cancer subtype is only available for
certain studies and so is not included at this stage.*

There are also: ovarian cancer samples, EPCAM+ cells, samples prepared
by spreading or spreading and samples from lymph nodes. For the time
being, I will not consider them.

``` r
suppressPackageStartupMessages(library(dplyr))
library(stringr)
library(biomaRt)
suppressPackageStartupMessages(library(tximport))
library(DT)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
library(cowplot)
suppressPackageStartupMessages(library(PCAtools))
library(dplyr)
suppressPackageStartupMessages(library(SummarizedExperiment))
```

``` r
EGAD_4810 <- read.table("/home/kevin/Documents/PhD/CAF_data/EGAD00001004810/delimited_maps/Run_Sample_meta_info.map", 
                        sep = ";")
EGAD_4810_cancers <- str_split_fixed(EGAD_4810$V1, pattern = "=", n = 2)[,2]
EGAD_4810_keep <- which(EGAD_4810_cancers == "BC")
EGAD_4810_filtered <- EGAD_4810[EGAD_4810_keep,]
EGAD_4810_meta <- data.frame(
  Sample = str_split_fixed(EGAD_4810_filtered$V4, pattern = "=", n = 2)[,2],
  Study = "EGAD00001004810",
  Subtype = "S3",
  Tumor_JuxtaTumor = tolower(str_split_fixed(EGAD_4810_filtered$V3, pattern = " ", n = 2)[,2]),
  directory = "/home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon",
  row.names = 1
)
EGAD_3808 <- read.table("/home/kevin/Documents/PhD/CAF_data/EGAD00001003808/meta_CAF-S1_S4_BC_47samples.txt", 
                        header = T, sep = "\t")
EGAD_3808_meta <- data.frame(
  Sample = EGAD_3808$Sample.Name,
  Study = "EGAD00001003808",
  Subtype = EGAD_3808$subset,
  Tumor_JuxtaTumor = EGAD_3808$Type, 
    directory = "/home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon",
  row.names = 1
)
EGAD_6144 <- read.table("/home/kevin/Documents/PhD/CAF_data/EGAD00001006144/meta_7samples.txt", 
                        header = T,sep = "\t")
EGAD_6144_meta <- data.frame(
  Sample = paste("CAF_Culture_", EGAD_6144$Sample.Name, sep = ""),
  Study = "EGAD00001006144",
  Subtype = "S1",
  Tumor_JuxtaTumor = "tumor",
    directory = "/home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon",
  row.names = 1
)
EGAD_5744 <- read.table("/home/kevin/Documents/PhD/CAF_data/EGAD00001005744/metaData_Pelon_et_al.txt", 
                        header =T, check.names = F)
EGAD_5744$Sample.Name <- gsub(pattern = "\\.", replacement = "-", x = EGAD_5744$Sample.Name )
EGAD_5744_filtered <- EGAD_5744[!(EGAD_5744$subset == "EPCAM+") & (EGAD_5744$Type == "T"),]
EGAD_5744_meta <- data.frame(
  Sample = EGAD_5744_filtered$Sample.Name,
  Study = "EGAD00001005744",
  Subtype = EGAD_5744_filtered$subset,
  Tumor_JuxtaTumor = "tumor",
    directory = "/home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon",
  row.names = 1
)
barkley_samples <- read.csv("/home/kevin/Documents/PhD/rna_seq_bc/metadata/reformat_samples.csv", 
                            header = T, row.names = "samples", check.names = F)
barkley_samples_meta <- data.frame(
  Sample = row.names(barkley_samples),
  Study = "In-House",
  Subtype = "Unknown",
  Tumor_JuxtaTumor = ifelse(barkley_samples$Condition == "Tumour", "tumor", "juxta-tumor"),
  directory = "/home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/",
  row.names = 1
)
metadata <- rbind.data.frame(EGAD_4810_meta, EGAD_3808_meta, 
                             EGAD_6144_meta, EGAD_5744_meta, barkley_samples_meta)
rownames(metadata) <- gsub(pattern = "-", replacement = "\\.", x = rownames(metadata))
# drop directory
metadata$directory <- NULL
metadata[1:5,]
```

    ##                  Study Subtype Tumor_JuxtaTumor
    ## B73T39 EGAD00001004810      S3            tumor
    ## B86T3  EGAD00001004810      S3            tumor
    ## B86T7  EGAD00001004810      S3      juxta-tumor
    ## B86T10 EGAD00001004810      S3            tumor
    ## B86T13 EGAD00001004810      S3      juxta-tumor

### Prepare transcript annotations

``` r
#mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="uswest.ensembl.org")
#tx2gene <- getBM(attributes = c("ensembl_transcript_id_version", "hgnc_symbol"), mart = mart, useCache = FALSE)
```

### Read in files

``` r
EGAD00001004810_rds <- readRDS("/home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon/salmon.merged.gene_counts.rds")
EGAD00001004810_tpm <- assays(EGAD00001004810_rds)$abundance
EGAD00001003808_rds <- readRDS("/home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon/salmon.merged.gene_counts.rds")
EGAD00001003808_tpm <- assays(EGAD00001003808_rds)$abundance
EGAD00001006144_rds <- readRDS("/home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon/salmon.merged.gene_counts.rds")
EGAD00001006144_tpm <- assays(EGAD00001006144_rds)$abundance
EGAD00001005744_rds <- readRDS("/home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon/salmon.merged.gene_counts.rds")
EGAD00001005744_tpm <- assays(EGAD00001005744_rds)$abundance
inhouse_rds <- readRDS("/home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/salmon.merged.gene_counts.rds")
inhouse_tpm <- assays(inhouse_rds)$abundance
samps <- c("4810", "3808", "6144", "5744", "inhouse")
dfs <- list(EGAD00001004810_tpm, EGAD00001003808_tpm, EGAD00001006144_tpm, EGAD00001005744_tpm, inhouse_tpm)
```

The inhouse samples were processed with a different version of
nf-core/rnaseq, and so have a different number of genes. Here we test
the pairwise overlap of the gene names for all samples

``` r
mat <- matrix(data = rep(0, length(dfs)*length(dfs)), 
              nrow = length(dfs), ncol = length(dfs), dimnames = list(samps,samps))
mat
```

    ##         4810 3808 6144 5744 inhouse
    ## 4810       0    0    0    0       0
    ## 3808       0    0    0    0       0
    ## 6144       0    0    0    0       0
    ## 5744       0    0    0    0       0
    ## inhouse    0    0    0    0       0

``` r
for (i in 1:length(dfs)){
  for (j in 1:length(dfs)){
    mat[i,j] <- length(intersect(rownames(dfs[[i]]), rownames(dfs[[j]])))
  }
}
mat
```

    ##          4810  3808  6144  5744 inhouse
    ## 4810    60603 60603 60603 60603   56388
    ## 3808    60603 60603 60603 60603   56388
    ## 6144    60603 60603 60603 60603   56388
    ## 5744    60603 60603 60603 60603   56388
    ## inhouse 56388 56388 56388 56388   60728

We want to have the same genes in each matrix

``` r
data_combined <- Reduce(function(x, y) merge(x, y, all=TRUE, by="rn", suffixes=c("", ".2")), 
    lapply(dfs, function(x) data.frame(x, rn = row.names(x))))
# get rid of genes with any NA
data_combined <- na.omit(data_combined)
# get rid of genes not expressed in any samples
thresh <- 0.0001
data_combined <- data_combined[rowSums(data_combined[,-1]) > thresh, ]
data_combined <- data_combined %>% remove_rownames %>% column_to_rownames(var="rn")
# get rid of X in colnames
colnames(data_combined) <- gsub(pattern = "X", replacement = "", x = colnames(data_combined))
# only want to keep samples that we want to analyse - those in the metadata object
data_combined <- data_combined[,intersect(rownames(metadata), colnames(data_combined))]
```

For PCA, we want colnames of our `data_combined` object to be the sample
names and the rownames of our `metadata` to be the sample names

``` r
all(colnames(data_combined) == rownames(metadata))
```

    ## [1] TRUE

### PCA

``` r
p <- pca(data_combined, metadata = metadata, removeVar = 0.1)
```

    ## -- removing the lower 10% of variables based on variance

``` r
biplot(p, colby = 'Study', legendPosition = 'top')
```

    ## Warning: ggrepel: 95 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](caf_rnaseq_combined_analysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
eigencorplot(p,metavars = c('Study', 'Subtype', 'Tumor_JuxtaTumor'))
```

    ## Warning in eigencorplot(p, metavars = c("Study", "Subtype",
    ## "Tumor_JuxtaTumor")): Study is not numeric - please check the source data as
    ## non-numeric variables will be coerced to numeric

    ## Warning in eigencorplot(p, metavars = c("Study", "Subtype",
    ## "Tumor_JuxtaTumor")): Subtype is not numeric - please check the source data as
    ## non-numeric variables will be coerced to numeric

    ## Warning in eigencorplot(p, metavars = c("Study", "Subtype",
    ## "Tumor_JuxtaTumor")): Tumor_JuxtaTumor is not numeric - please check the source
    ## data as non-numeric variables will be coerced to numeric

![](caf_rnaseq_combined_analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

It is not surprising that the first principal component is correlated
with the study of origin. There is also a correlation with subtype, but
this could also be due to the link with study - different studies tend
to have different subtypes.

Next step is to look at batch correction.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Pelon2020" class="csl-entry">

Pelon, Floriane, Brigitte Bourachot, Yann Kieffer, Ilaria Magagna, Fanny
Mermet-Meillon, Isabelle Bonnet, Ana Costa, et al. 2020. “<span
class="nocase">Cancer-associated fibroblast heterogeneity in axillary
lymph nodes drives metastases in breast cancer through complementary
mechanisms</span>.” *Nature Communications 2020 11:1* 11 (1): 1–20.
<https://doi.org/10.1038/s41467-019-14134-w>.

</div>

</div>
