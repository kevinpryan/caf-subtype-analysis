CAF subtype analysis
================
Kevin Ryan
2022-06-28 21:30:45

-   <a href="#introduction" id="toc-introduction">Introduction</a>
    -   <a href="#preparation" id="toc-preparation">Preparation</a>
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
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(stringr)
library(biomaRt)
library(tximport)
```

``` r
EGAD_4810 <- read.table("/home/kevin/Documents/PhD/CAF_data/EGAD00001004810/delimited_maps/Run_Sample_meta_info.map", sep = ";")
EGAD_4810_cancers <- str_split_fixed(EGAD_4810$V1, pattern = "=", n = 2)[,2]
EGAD_4810_keep <- which(EGAD_4810_cancers == "BC")
EGAD_4810_filtered <- EGAD_4810[EGAD_4810_keep,]
EGAD_4810_filtered
```

    ##                 V1            V2                       V3                  V4
    ## 1  disease_site=BC gender=female       phenotype=S3 Tumor   subject_id=B73T39
    ## 2  disease_site=BC gender=female       phenotype=S3 Tumor    subject_id=B86T3
    ## 3  disease_site=BC gender=female phenotype=S3 Juxta-Tumor    subject_id=B86T7
    ## 4  disease_site=BC gender=female       phenotype=S3 Tumor   subject_id=B86T10
    ## 5  disease_site=BC gender=female phenotype=S3 Juxta-Tumor   subject_id=B86T13
    ## 6  disease_site=BC gender=female phenotype=S3 Juxta-Tumor   subject_id=B86T16
    ## 7  disease_site=BC gender=female       phenotype=S3 Tumor   subject_id=B86T22
    ## 8  disease_site=BC gender=female       phenotype=S3 Tumor subject_id=B103T103
    ## 9  disease_site=BC gender=female phenotype=S3 Juxta-Tumor subject_id=B103T107
    ## 10 disease_site=BC gender=female       phenotype=S3 Tumor subject_id=B103T111
    ## 11 disease_site=BC gender=female phenotype=S3 Juxta-Tumor subject_id=B103T115
    ## 12 disease_site=BC gender=female       phenotype=S3 Tumor  subject_id=B103T51
    ## 13 disease_site=BC gender=female phenotype=S3 Juxta-Tumor  subject_id=B103T55
    ## 14 disease_site=BC gender=female       phenotype=S3 Tumor  subject_id=B103T67
    ## 15 disease_site=BC gender=female       phenotype=S3 Tumor  subject_id=B103T71
    ## 16 disease_site=BC gender=female phenotype=S3 Juxta-Tumor  subject_id=B103T75
    ## 17 disease_site=BC gender=female       phenotype=S3 Tumor  subject_id=B103T79
    ## 18 disease_site=BC gender=female phenotype=S3 Juxta-Tumor  subject_id=B103T83
    ## 19 disease_site=BC gender=female       phenotype=S3 Tumor  subject_id=B103T87
    ## 20 disease_site=BC gender=female phenotype=S3 Juxta-Tumor  subject_id=B103T91
    ## 21 disease_site=BC gender=female       phenotype=S3 Tumor  subject_id=B103T95
    ## 22 disease_site=BC gender=female phenotype=S3 Juxta-Tumor  subject_id=B103T99
    ## 25 disease_site=BC gender=female       phenotype=S3 Tumor   subject_id=B123U2
    ## 26 disease_site=BC gender=female       phenotype=S3 Tumor   subject_id=B123U6
    ## 27 disease_site=BC gender=female phenotype=S3 Juxta-Tumor   subject_id=B123U8
    ##                         V5 V6
    ## 1  ENA-CHECKLIST=ERC000026 NA
    ## 2  ENA-CHECKLIST=ERC000026 NA
    ## 3  ENA-CHECKLIST=ERC000026 NA
    ## 4  ENA-CHECKLIST=ERC000026 NA
    ## 5  ENA-CHECKLIST=ERC000026 NA
    ## 6  ENA-CHECKLIST=ERC000026 NA
    ## 7  ENA-CHECKLIST=ERC000026 NA
    ## 8  ENA-CHECKLIST=ERC000026 NA
    ## 9  ENA-CHECKLIST=ERC000026 NA
    ## 10 ENA-CHECKLIST=ERC000026 NA
    ## 11 ENA-CHECKLIST=ERC000026 NA
    ## 12 ENA-CHECKLIST=ERC000026 NA
    ## 13 ENA-CHECKLIST=ERC000026 NA
    ## 14 ENA-CHECKLIST=ERC000026 NA
    ## 15 ENA-CHECKLIST=ERC000026 NA
    ## 16 ENA-CHECKLIST=ERC000026 NA
    ## 17 ENA-CHECKLIST=ERC000026 NA
    ## 18 ENA-CHECKLIST=ERC000026 NA
    ## 19 ENA-CHECKLIST=ERC000026 NA
    ## 20 ENA-CHECKLIST=ERC000026 NA
    ## 21 ENA-CHECKLIST=ERC000026 NA
    ## 22 ENA-CHECKLIST=ERC000026 NA
    ## 25 ENA-CHECKLIST=ERC000026 NA
    ## 26 ENA-CHECKLIST=ERC000026 NA
    ## 27 ENA-CHECKLIST=ERC000026 NA

``` r
EGAD_4810_meta <- data.frame(
  Sample = str_split_fixed(EGAD_4810_filtered$V4, pattern = "=", n = 2)[,2],
  Study = "EGAD00001004810",
  Subtype = "S3",
  Tumor_JuxtaTumor = tolower(str_split_fixed(EGAD_4810_filtered$V3, pattern = " ", n = 2)[,2]),
  directory = "/home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon",
  row.names = 1
)
EGAD_4810_meta
```

    ##                    Study Subtype Tumor_JuxtaTumor
    ## B73T39   EGAD00001004810      S3            tumor
    ## B86T3    EGAD00001004810      S3            tumor
    ## B86T7    EGAD00001004810      S3      juxta-tumor
    ## B86T10   EGAD00001004810      S3            tumor
    ## B86T13   EGAD00001004810      S3      juxta-tumor
    ## B86T16   EGAD00001004810      S3      juxta-tumor
    ## B86T22   EGAD00001004810      S3            tumor
    ## B103T103 EGAD00001004810      S3            tumor
    ## B103T107 EGAD00001004810      S3      juxta-tumor
    ## B103T111 EGAD00001004810      S3            tumor
    ## B103T115 EGAD00001004810      S3      juxta-tumor
    ## B103T51  EGAD00001004810      S3            tumor
    ## B103T55  EGAD00001004810      S3      juxta-tumor
    ## B103T67  EGAD00001004810      S3            tumor
    ## B103T71  EGAD00001004810      S3            tumor
    ## B103T75  EGAD00001004810      S3      juxta-tumor
    ## B103T79  EGAD00001004810      S3            tumor
    ## B103T83  EGAD00001004810      S3      juxta-tumor
    ## B103T87  EGAD00001004810      S3            tumor
    ## B103T91  EGAD00001004810      S3      juxta-tumor
    ## B103T95  EGAD00001004810      S3            tumor
    ## B103T99  EGAD00001004810      S3      juxta-tumor
    ## B123U2   EGAD00001004810      S3            tumor
    ## B123U6   EGAD00001004810      S3            tumor
    ## B123U8   EGAD00001004810      S3      juxta-tumor
    ##                                                                                             directory
    ## B73T39   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B86T3    /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B86T7    /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B86T10   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B86T13   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B86T16   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B86T22   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T103 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T107 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T111 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T115 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T51  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T55  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T67  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T71  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T75  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T79  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T83  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T87  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T91  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T95  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T99  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B123U2   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B123U6   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B123U8   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon

``` r
EGAD_3808 <- read.table("/home/kevin/Documents/PhD/CAF_data/EGAD00001003808/meta_CAF-S1_S4_BC_47samples.txt", header = T, sep = "\t")
EGAD_3808
```

    ##    Sample.Name        Type subset   Cancer.subtype
    ## 1       B73T37       tumor     S1             LumA
    ## 2       B73T40       tumor     S4             LumA
    ## 3        B86T1       tumor     S1               TN
    ## 4        B86T4       tumor     S4               TN
    ## 5        B86T5 juxta-tumor     S1               TN
    ## 6        B86T8 juxta-tumor     S4               TN
    ## 7        B86T9       tumor     S1             LumA
    ## 8       B86T11       tumor     S4             LumA
    ## 9       B86T12 juxta-tumor     S1             LumA
    ## 10      B86T14 juxta-tumor     S4             LumA
    ## 11      B86T15 juxta-tumor     S1             LumA
    ## 12      B86T17 juxta-tumor     S4             LumA
    ## 13      B86T18       tumor     S1             LumA
    ## 14      B86T21 juxta-tumor     S1             LumA
    ## 15      B86T23       tumor     S4             LumA
    ## 16      B86T24       tumor     S1             LumA
    ## 17      B86T26 juxta-tumor     S4             LumA
    ## 18    B103T100 juxta-tumor     S4               TN
    ## 19    B103T101       tumor     S1             LumA
    ## 20    B103T104       tumor     S4             LumA
    ## 21    B103T105 juxta-tumor     S1             LumA
    ## 22    B103T108 juxta-tumor     S4             LumA
    ## 23    B103T109       tumor     S1             LumA
    ## 24    B103T112       tumor     S4             LumA
    ## 25     B103T49       tumor     S1 Mix LumA / Lum B
    ## 26     B103T52       tumor     S4 Mix LumA / Lum B
    ## 27     B103T53 juxta-tumor     S1 Mix LumA / Lum B
    ## 28     B103T56 juxta-tumor     S4 Mix LumA / Lum B
    ## 29     B103T57       tumor     S1             LumA
    ## 30     B103T61 juxta-tumor     S1             LumA
    ## 31     B103T65       tumor     S1             HER2
    ## 32     B103T69       tumor     S1             LumA
    ## 33     B103T72       tumor     S4             LumA
    ## 34     B103T73 juxta-tumor     S1             LumA
    ## 35     B103T76 juxta-tumor     S4             LumA
    ## 36     B103T77       tumor     S1               TN
    ## 37     B103T80       tumor     S4               TN
    ## 38     B103T81 juxta-tumor     S1               TN
    ## 39     B103T84 juxta-tumor     S4               TN
    ## 40     B103T85       tumor     S1               TN
    ## 41     B103T93       tumor     S1               TN
    ## 42     B103T96       tumor     S4               TN
    ## 43     B103T97 juxta-tumor     S1               TN
    ## 44      B123U1       tumor     S1               TN
    ## 45      B123U3 juxta-tumor     S1               TN
    ## 46      B123U4       tumor     S1               TN
    ## 47      B123U7 juxta-tumor     S1               TN

``` r
EGAD_3808_meta <- data.frame(
  Sample = EGAD_3808$Sample.Name,
  Study = "EGAD00001003808",
  Subtype = EGAD_3808$subset,
  Tumor_JuxtaTumor = EGAD_3808$Type, 
    directory = "/home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon",
  row.names = 1
)
EGAD_3808_meta
```

    ##                    Study Subtype Tumor_JuxtaTumor
    ## B73T37   EGAD00001003808      S1            tumor
    ## B73T40   EGAD00001003808      S4            tumor
    ## B86T1    EGAD00001003808      S1            tumor
    ## B86T4    EGAD00001003808      S4            tumor
    ## B86T5    EGAD00001003808      S1      juxta-tumor
    ## B86T8    EGAD00001003808      S4      juxta-tumor
    ## B86T9    EGAD00001003808      S1            tumor
    ## B86T11   EGAD00001003808      S4            tumor
    ## B86T12   EGAD00001003808      S1      juxta-tumor
    ## B86T14   EGAD00001003808      S4      juxta-tumor
    ## B86T15   EGAD00001003808      S1      juxta-tumor
    ## B86T17   EGAD00001003808      S4      juxta-tumor
    ## B86T18   EGAD00001003808      S1            tumor
    ## B86T21   EGAD00001003808      S1      juxta-tumor
    ## B86T23   EGAD00001003808      S4            tumor
    ## B86T24   EGAD00001003808      S1            tumor
    ## B86T26   EGAD00001003808      S4      juxta-tumor
    ## B103T100 EGAD00001003808      S4      juxta-tumor
    ## B103T101 EGAD00001003808      S1            tumor
    ## B103T104 EGAD00001003808      S4            tumor
    ## B103T105 EGAD00001003808      S1      juxta-tumor
    ## B103T108 EGAD00001003808      S4      juxta-tumor
    ## B103T109 EGAD00001003808      S1            tumor
    ## B103T112 EGAD00001003808      S4            tumor
    ## B103T49  EGAD00001003808      S1            tumor
    ## B103T52  EGAD00001003808      S4            tumor
    ## B103T53  EGAD00001003808      S1      juxta-tumor
    ## B103T56  EGAD00001003808      S4      juxta-tumor
    ## B103T57  EGAD00001003808      S1            tumor
    ## B103T61  EGAD00001003808      S1      juxta-tumor
    ## B103T65  EGAD00001003808      S1            tumor
    ## B103T69  EGAD00001003808      S1            tumor
    ## B103T72  EGAD00001003808      S4            tumor
    ## B103T73  EGAD00001003808      S1      juxta-tumor
    ## B103T76  EGAD00001003808      S4      juxta-tumor
    ## B103T77  EGAD00001003808      S1            tumor
    ## B103T80  EGAD00001003808      S4            tumor
    ## B103T81  EGAD00001003808      S1      juxta-tumor
    ## B103T84  EGAD00001003808      S4      juxta-tumor
    ## B103T85  EGAD00001003808      S1            tumor
    ## B103T93  EGAD00001003808      S1            tumor
    ## B103T96  EGAD00001003808      S4            tumor
    ## B103T97  EGAD00001003808      S1      juxta-tumor
    ## B123U1   EGAD00001003808      S1            tumor
    ## B123U3   EGAD00001003808      S1      juxta-tumor
    ## B123U4   EGAD00001003808      S1            tumor
    ## B123U7   EGAD00001003808      S1      juxta-tumor
    ##                                                                                             directory
    ## B73T37   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B73T40   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T1    /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T4    /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T5    /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T8    /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T9    /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T11   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T12   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T14   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T15   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T17   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T18   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T21   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T23   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T24   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T26   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T100 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T101 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T104 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T105 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T108 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T109 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T112 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T49  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T52  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T53  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T56  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T57  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T61  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T65  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T69  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T72  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T73  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T76  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T77  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T80  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T81  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T84  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T85  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T93  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T96  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T97  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B123U1   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B123U3   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B123U4   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B123U7   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon

``` r
EGAD_6144 <- read.table("/home/kevin/Documents/PhD/CAF_data/EGAD00001006144/meta_7samples.txt", header = T,sep = "\t")
EGAD_6144_meta <- data.frame(
  Sample = paste("CAF_Culture_", EGAD_6144$Sample.Name, sep = ""),
  Study = "EGAD00001006144",
  Subtype = "S1",
  Tumor_JuxtaTumor = "tumor",
    directory = "/home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon",
  row.names = 1
)
EGAD_6144_meta
```

    ##                               Study Subtype Tumor_JuxtaTumor
    ## CAF_Culture_D220T13 EGAD00001006144      S1            tumor
    ## CAF_Culture_D220T17 EGAD00001006144      S1            tumor
    ## CAF_Culture_D220T21 EGAD00001006144      S1            tumor
    ## CAF_Culture_D220T25 EGAD00001006144      S1            tumor
    ## CAF_Culture_D220T29 EGAD00001006144      S1            tumor
    ## CAF_Culture_D220T33 EGAD00001006144      S1            tumor
    ## CAF_Culture_D220T37 EGAD00001006144      S1            tumor
    ##                                                                                                        directory
    ## CAF_Culture_D220T13 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon
    ## CAF_Culture_D220T17 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon
    ## CAF_Culture_D220T21 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon
    ## CAF_Culture_D220T25 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon
    ## CAF_Culture_D220T29 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon
    ## CAF_Culture_D220T33 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon
    ## CAF_Culture_D220T37 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon

``` r
EGAD_5744 <- read.table("/home/kevin/Documents/PhD/CAF_data/EGAD00001005744/metaData_Pelon_et_al.txt", header =T, check.names = F)
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
EGAD_5744_meta
```

    ##                             Study Subtype Tumor_JuxtaTumor
    ## A461-A462-A465U15 EGAD00001005744      S1            tumor
    ## A461-A462-A465U17 EGAD00001005744      S4            tumor
    ## A461-A462-A465U25 EGAD00001005744      S1            tumor
    ## A461-A462-A465U27 EGAD00001005744      S4            tumor
    ## A461-A462-A465U3  EGAD00001005744      S1            tumor
    ## A461-A462-A465U32 EGAD00001005744      S1            tumor
    ## A461-A462-A465U35 EGAD00001005744      S4            tumor
    ## A461-A462-A465U5  EGAD00001005744      S4            tumor
    ## A461-A462-A465U7  EGAD00001005744      S1            tumor
    ## A461-A462-A465U9  EGAD00001005744      S4            tumor
    ##                                                                                                      directory
    ## A461-A462-A465U15 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U17 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U25 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U27 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U3  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U32 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U35 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U5  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U7  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U9  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon

``` r
barkley_samples <- read.csv("/home/kevin/Documents/PhD/rna_seq_bc/metadata/reformat_samples.csv", header = T, row.names = "samples", check.names = F)
barkley_samples_meta <- data.frame(
  Sample = row.names(barkley_samples),
  Study = "In-House",
  Subtype = "Unknown",
  Tumor_JuxtaTumor = ifelse(barkley_samples$Condition == "Tumour", "tumor", "juxta-tumor"),
  directory = "/home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/",
  row.names = 1
)
barkley_samples_meta
```

    ##         Study Subtype Tumor_JuxtaTumor
    ## 4033 In-House Unknown            tumor
    ## 4034 In-House Unknown      juxta-tumor
    ## 4027 In-House Unknown            tumor
    ## 4028 In-House Unknown      juxta-tumor
    ## 4112 In-House Unknown            tumor
    ## 4113 In-House Unknown      juxta-tumor
    ## 4116 In-House Unknown            tumor
    ## 4117 In-House Unknown      juxta-tumor
    ## 4214 In-House Unknown            tumor
    ## 4215 In-House Unknown      juxta-tumor
    ## 4315 In-House Unknown            tumor
    ## 4316 In-House Unknown      juxta-tumor
    ## 4340 In-House Unknown            tumor
    ## 4341 In-House Unknown      juxta-tumor
    ## 4344 In-House Unknown            tumor
    ## 4345 In-House Unknown      juxta-tumor
    ## 3532 In-House Unknown            tumor
    ## 3533 In-House Unknown      juxta-tumor
    ## 3536 In-House Unknown            tumor
    ## 3537 In-House Unknown      juxta-tumor
    ## 4299 In-House Unknown            tumor
    ## 4300 In-House Unknown      juxta-tumor
    ## 4722 In-House Unknown            tumor
    ## 4723 In-House Unknown      juxta-tumor
    ##                                                                                             directory
    ## 4033 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4034 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4027 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4028 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4112 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4113 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4116 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4117 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4214 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4215 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4315 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4316 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4340 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4341 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4344 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4345 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 3532 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 3533 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 3536 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 3537 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4299 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4300 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4722 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4723 /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/

``` r
metadata <- rbind.data.frame(EGAD_4810_meta, EGAD_3808_meta, EGAD_6144_meta, EGAD_5744_meta, barkley_samples_meta)
metadata
```

    ##                               Study Subtype Tumor_JuxtaTumor
    ## B73T39              EGAD00001004810      S3            tumor
    ## B86T3               EGAD00001004810      S3            tumor
    ## B86T7               EGAD00001004810      S3      juxta-tumor
    ## B86T10              EGAD00001004810      S3            tumor
    ## B86T13              EGAD00001004810      S3      juxta-tumor
    ## B86T16              EGAD00001004810      S3      juxta-tumor
    ## B86T22              EGAD00001004810      S3            tumor
    ## B103T103            EGAD00001004810      S3            tumor
    ## B103T107            EGAD00001004810      S3      juxta-tumor
    ## B103T111            EGAD00001004810      S3            tumor
    ## B103T115            EGAD00001004810      S3      juxta-tumor
    ## B103T51             EGAD00001004810      S3            tumor
    ## B103T55             EGAD00001004810      S3      juxta-tumor
    ## B103T67             EGAD00001004810      S3            tumor
    ## B103T71             EGAD00001004810      S3            tumor
    ## B103T75             EGAD00001004810      S3      juxta-tumor
    ## B103T79             EGAD00001004810      S3            tumor
    ## B103T83             EGAD00001004810      S3      juxta-tumor
    ## B103T87             EGAD00001004810      S3            tumor
    ## B103T91             EGAD00001004810      S3      juxta-tumor
    ## B103T95             EGAD00001004810      S3            tumor
    ## B103T99             EGAD00001004810      S3      juxta-tumor
    ## B123U2              EGAD00001004810      S3            tumor
    ## B123U6              EGAD00001004810      S3            tumor
    ## B123U8              EGAD00001004810      S3      juxta-tumor
    ## B73T37              EGAD00001003808      S1            tumor
    ## B73T40              EGAD00001003808      S4            tumor
    ## B86T1               EGAD00001003808      S1            tumor
    ## B86T4               EGAD00001003808      S4            tumor
    ## B86T5               EGAD00001003808      S1      juxta-tumor
    ## B86T8               EGAD00001003808      S4      juxta-tumor
    ## B86T9               EGAD00001003808      S1            tumor
    ## B86T11              EGAD00001003808      S4            tumor
    ## B86T12              EGAD00001003808      S1      juxta-tumor
    ## B86T14              EGAD00001003808      S4      juxta-tumor
    ## B86T15              EGAD00001003808      S1      juxta-tumor
    ## B86T17              EGAD00001003808      S4      juxta-tumor
    ## B86T18              EGAD00001003808      S1            tumor
    ## B86T21              EGAD00001003808      S1      juxta-tumor
    ## B86T23              EGAD00001003808      S4            tumor
    ## B86T24              EGAD00001003808      S1            tumor
    ## B86T26              EGAD00001003808      S4      juxta-tumor
    ## B103T100            EGAD00001003808      S4      juxta-tumor
    ## B103T101            EGAD00001003808      S1            tumor
    ## B103T104            EGAD00001003808      S4            tumor
    ## B103T105            EGAD00001003808      S1      juxta-tumor
    ## B103T108            EGAD00001003808      S4      juxta-tumor
    ## B103T109            EGAD00001003808      S1            tumor
    ## B103T112            EGAD00001003808      S4            tumor
    ## B103T49             EGAD00001003808      S1            tumor
    ## B103T52             EGAD00001003808      S4            tumor
    ## B103T53             EGAD00001003808      S1      juxta-tumor
    ## B103T56             EGAD00001003808      S4      juxta-tumor
    ## B103T57             EGAD00001003808      S1            tumor
    ## B103T61             EGAD00001003808      S1      juxta-tumor
    ## B103T65             EGAD00001003808      S1            tumor
    ## B103T69             EGAD00001003808      S1            tumor
    ## B103T72             EGAD00001003808      S4            tumor
    ## B103T73             EGAD00001003808      S1      juxta-tumor
    ## B103T76             EGAD00001003808      S4      juxta-tumor
    ## B103T77             EGAD00001003808      S1            tumor
    ## B103T80             EGAD00001003808      S4            tumor
    ## B103T81             EGAD00001003808      S1      juxta-tumor
    ## B103T84             EGAD00001003808      S4      juxta-tumor
    ## B103T85             EGAD00001003808      S1            tumor
    ## B103T93             EGAD00001003808      S1            tumor
    ## B103T96             EGAD00001003808      S4            tumor
    ## B103T97             EGAD00001003808      S1      juxta-tumor
    ## B123U1              EGAD00001003808      S1            tumor
    ## B123U3              EGAD00001003808      S1      juxta-tumor
    ## B123U4              EGAD00001003808      S1            tumor
    ## B123U7              EGAD00001003808      S1      juxta-tumor
    ## CAF_Culture_D220T13 EGAD00001006144      S1            tumor
    ## CAF_Culture_D220T17 EGAD00001006144      S1            tumor
    ## CAF_Culture_D220T21 EGAD00001006144      S1            tumor
    ## CAF_Culture_D220T25 EGAD00001006144      S1            tumor
    ## CAF_Culture_D220T29 EGAD00001006144      S1            tumor
    ## CAF_Culture_D220T33 EGAD00001006144      S1            tumor
    ## CAF_Culture_D220T37 EGAD00001006144      S1            tumor
    ## A461-A462-A465U15   EGAD00001005744      S1            tumor
    ## A461-A462-A465U17   EGAD00001005744      S4            tumor
    ## A461-A462-A465U25   EGAD00001005744      S1            tumor
    ## A461-A462-A465U27   EGAD00001005744      S4            tumor
    ## A461-A462-A465U3    EGAD00001005744      S1            tumor
    ## A461-A462-A465U32   EGAD00001005744      S1            tumor
    ## A461-A462-A465U35   EGAD00001005744      S4            tumor
    ## A461-A462-A465U5    EGAD00001005744      S4            tumor
    ## A461-A462-A465U7    EGAD00001005744      S1            tumor
    ## A461-A462-A465U9    EGAD00001005744      S4            tumor
    ## 4033                       In-House Unknown            tumor
    ## 4034                       In-House Unknown      juxta-tumor
    ## 4027                       In-House Unknown            tumor
    ## 4028                       In-House Unknown      juxta-tumor
    ## 4112                       In-House Unknown            tumor
    ## 4113                       In-House Unknown      juxta-tumor
    ## 4116                       In-House Unknown            tumor
    ## 4117                       In-House Unknown      juxta-tumor
    ## 4214                       In-House Unknown            tumor
    ## 4215                       In-House Unknown      juxta-tumor
    ## 4315                       In-House Unknown            tumor
    ## 4316                       In-House Unknown      juxta-tumor
    ## 4340                       In-House Unknown            tumor
    ## 4341                       In-House Unknown      juxta-tumor
    ## 4344                       In-House Unknown            tumor
    ## 4345                       In-House Unknown      juxta-tumor
    ## 3532                       In-House Unknown            tumor
    ## 3533                       In-House Unknown      juxta-tumor
    ## 3536                       In-House Unknown            tumor
    ## 3537                       In-House Unknown      juxta-tumor
    ## 4299                       In-House Unknown            tumor
    ## 4300                       In-House Unknown      juxta-tumor
    ## 4722                       In-House Unknown            tumor
    ## 4723                       In-House Unknown      juxta-tumor
    ##                                                                                                            directory
    ## B73T39                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B86T3                   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B86T7                   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B86T10                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B86T13                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B86T16                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B86T22                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T103                /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T107                /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T111                /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T115                /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T51                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T55                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T67                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T71                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T75                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T79                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T83                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T87                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T91                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T95                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B103T99                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B123U2                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B123U6                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B123U8                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon
    ## B73T37                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B73T40                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T1                   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T4                   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T5                   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T8                   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T9                   /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T11                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T12                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T14                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T15                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T17                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T18                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T21                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T23                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T24                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B86T26                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T100                /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T101                /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T104                /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T105                /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T108                /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T109                /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T112                /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T49                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T52                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T53                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T56                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T57                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T61                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T65                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T69                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T72                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T73                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T76                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T77                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T80                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T81                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T84                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T85                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T93                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T96                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B103T97                 /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B123U1                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B123U3                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B123U4                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## B123U7                  /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon
    ## CAF_Culture_D220T13     /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon
    ## CAF_Culture_D220T17     /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon
    ## CAF_Culture_D220T21     /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon
    ## CAF_Culture_D220T25     /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon
    ## CAF_Culture_D220T29     /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon
    ## CAF_Culture_D220T33     /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon
    ## CAF_Culture_D220T37     /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon
    ## A461-A462-A465U15       /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U17       /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U25       /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U27       /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U3        /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U32       /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U35       /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U5        /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U7        /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## A461-A462-A465U9        /home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon
    ## 4033                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4034                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4027                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4028                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4112                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4113                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4116                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4117                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4214                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4215                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4315                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4316                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4340                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4341                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4344                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4345                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 3532                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 3533                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 3536                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 3537                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4299                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4300                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4722                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/
    ## 4723                /home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon/

### Prepare transcript annotations

``` r
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="uswest.ensembl.org")
```

    ## Warning: Ensembl will soon enforce the use of https.
    ## Ensure the 'host' argument includes "https://"

``` r
tx2gene <- getBM(attributes = c("ensembl_transcript_id_version", "hgnc_symbol"), mart = mart, useCache = FALSE)
```

### Read in files

``` r
#basedirs <- c("/home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results/star_salmon", "/home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon", "/home/kevin/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon", "/home/kevin/Documents/PhD/CAF_data/nfcore_results/inhouse_caf_nfcore_rnaseq_results/star_salmon")
files <- file.path(metadata$directory, rownames(metadata), "quant.sf")
names(files) <- rownames(metadata)
file.exists(files)
```

    ##   [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [31] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [46] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [61] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [76] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [91] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [106] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

``` r
#txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
```

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
