CAF Subpopulation Analysis Edited
================
Kevin Ryan
2022-09-02 10:43:22

- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#preparation" id="toc-preparation">Preparation</a>

# Introduction

Cancer-associated fibroblasts (CAFs) are a heterogeneous cell type found
in the tumour microenvironment. They have a wide array of functions, and
tend to be immunosuppressive and cancer-promoting. There have been many
attempts to characterise subpopulations of CAFs, with much
transcriptomic analysis being carried out in the Mechta-Grigoriou lab in
Institut Curie. They have identified 4 ‘subpopulations’ which can be
separated based on the expression of different markers:

- S1: FAP<sup>High</sup>, CD29<sup>Med-High</sup>, α<sup>SMAHigh</sup>,
  PDPN<sup>High</sup>, PDGFRβ<sup>High</sup>
- S2: FAP<sup>Neg</sup>, CD29<sup>Low</sup>, αSMANeg-<sup>Low</sup>,
  PDPN<sup>Low</sup>, PDGFRβ<sup>Low</sup>
- S3: FAP<sup>Neg-Low</sup>, CD29<sup>Med</sup>, αSMA<sup>Neg-Low</sup>,
  PDPN<sup>Low</sup>, PDGFRβ<sup>Low-Med</sup>
- S4: FAP<sup>Low-Med</sup>, CD29<sup>High</sup>, αSMA<sup>High</sup>,
  PDPN<sup>Low</sup>, PDGFRβ<sup>Med</sup>

(Pelon et al. 2020)

FACS gating strategies can be used to isolate these various
subpopulations. The Mechta-Grigoriou group have done this and have
generated bulk RNA-sequencing data for the S1, S3 and S4 subpopulations.
They generated scRNA-sequencing data for the S1 subpopulation. This data
was deposited on the European Genome Phenome Archive, and was accessed
via a Data Transfer Agreement.

The following summarises the data obtained:

<table style="width:98%;">
<colgroup>
<col style="width: 17%" />
<col style="width: 17%" />
<col style="width: 26%" />
<col style="width: 35%" />
</colgroup>
<thead>
<tr class="header">
<th>Subpopulation</th>
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

With the juxta-tumour data, they got tumour and juxta-tumour samples
from the same patient. However, I have not been able to figure out
whether they came from the same patient. We could possibly use Optitype
to determine HLA allele - match tumour and juxta tumour.

We also have scRNA-seq data for S1, labelled with 8 subpopulations of S1
CAFs. It may be possible to

It is likely that sorting the cells using FACS alters the
transcriptional properties of the cells compared to if they are
separated using spreading approaches, as is seen in study
`EGAD00001006144`. This is something that we will have to keep in mind.

The data was processed using nf-core/rnaseq version `3.8.1` using the
default parameters. STAR/Salmon were used for alignment/quantification.

We would expect our tumour-associated normal to be most like the S3
subpopulation (usually accumulate in juxta-tumours). The S2
subpopulation has been found to accumulate more in luminal A breast
cancer, whereas the S4 subpopulation tends to be present in Her2+ breast
cancers. Unfortunately, data is not available for the S2 subpopulation
and 11 of the 12 cancers encountered in our samples are Luminal A.

Combining RNA-sequencing datasets from different studies can be very
challenging. We can expect batch effects to be present, so it might not
be possible to determine whether differences we observe are due to
actual biological effects or technical artifacts. In addition, a recent
study suggests that DESeq2 and edgeR (the most popular differential
expression tools) experience large rates of false positives when used
with large sample sizes (Li et al. 2022). However, this assertion has
been refuted, and it has been suggested that the Li 2022 study did not
apply appropriate batch correction and quality control ([Twitter
thread](https://threadreaderapp.com/thread/1513468597288452097.html)
from Mike Love and associated [code on
GitHub](https://github.com/mikelove/preNivolumabOnNivolumab/blob/main/preNivolumabOnNivolumab.knit.md)).
One of the datasets (`EGAD00001006144`) was produced using stranded
RNA-seq, whereas the other datasets were unstranded. This can lead to a
lack of comparability of the datasets (Zhao, Ye, and Stanton 2020). It
may be necessary to drop this dataset from the analysis. All samples
were prepared by poly(A) selection (use of oligo-dT).

# Preparation

Columns will be: Sample, Study, Subpopulation, Tumor_Juxtatumor

*Here we will be combining data from 5 studies. To begin with, we will
only include the metadata available for all studies (except for our
unknown CAF Subpopulation label). Breast cancer subtype is only
available for certain studies and so is not included at this stage.*

There are also: ovarian cancer samples, EPCAM+ cells (an epithelial
marker) and samples from lymph nodes. For the time being, I will not
consider them.

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Li2022" class="csl-entry">

Li, Yumei, Xinzhou Ge, Fanglue Peng, Wei Li, and Jingyi Jessica Li.
2022. “<span class="nocase">Exaggerated false positives by popular
differential expression methods when analyzing human population
samples</span>.” *Genome Biology* 23 (1): 1–13.
<https://doi.org/10.1186/S13059-022-02648-4/FIGURES/2>.

</div>

<div id="ref-Pelon2020" class="csl-entry">

Pelon, Floriane, Brigitte Bourachot, Yann Kieffer, Ilaria Magagna, Fanny
Mermet-Meillon, Isabelle Bonnet, Ana Costa, et al. 2020. “<span
class="nocase">Cancer-associated fibroblast heterogeneity in axillary
lymph nodes drives metastases in breast cancer through complementary
mechanisms</span>.” *Nature Communications 2020 11:1* 11 (1): 1–20.
<https://doi.org/10.1038/s41467-019-14134-w>.

</div>

<div id="ref-Zhao2020" class="csl-entry">

Zhao, Shanrong, Zhan Ye, and Robert Stanton. 2020. “<span
class="nocase">Misuse of RPKM or TPM normalization when comparing across
samples and sequencing protocols</span>.” *RNA* 26 (8): 903.
<https://doi.org/10.1261/RNA.074922.120>.

</div>

</div>
