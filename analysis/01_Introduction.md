CAF Subpopulation Analysis
================

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

FACS gating strategies can be used to isolate these subpopulations. The
Mechta-Grigoriou group have done this and have generated bulk
RNA-sequencing data for the S1, S3 and S4 subpopulations. They generated
scRNA-sequencing data for the S1 subpopulation. This data was deposited
on the European Genome Phenome Archive, and was accessed via a Data
Transfer Agreement.

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

With the juxta-tumour data, tumour and juxta-tumour samples came from
the same patient. However, the metadata gives no indication of these
pairings. We could possibly use Optitype (Szolek et al. 2014) to
determine HLA allele and match the tumour and juxta-tumour samples.

We also have scRNA-seq data for S1, labelled with 8 subpopulations of S1
CAFs. It may be possible to use CIBERSORT (Newman et al. 2015) and
BayesPrism (Chu et al. 2022) to deconvolve the bulk S1 RNA-sequencing
data to further confirm the presence of these subpopulations.

It is likely that sorting the cells using FACS alters the
transcriptional properties of the cells compared to if they are
separated using spreading approaches, as is seen in study
`EGAD00001006144` and described in (Kieffer et al. 2020). This is
something that we will have to keep in mind.

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
been refuted, and it has been implied that the Li 2022 study did not
apply appropriate batch correction and quality control ([Twitter
thread](https://threadreaderapp.com/thread/1513468597288452097.html)
from Mike Love and associated [code on
GitHub](https://github.com/mikelove/preNivolumabOnNivolumab/blob/main/preNivolumabOnNivolumab.knit.md)).
One of the datasets (`EGAD00001006144`) was produced using stranded
RNA-seq, whereas the other datasets were unstranded. This can lead to a
lack of comparability of the datasets (Zhao, Ye, and Stanton 2020). It
may be necessary to drop this dataset from the analysis. All samples
were prepared by poly(A) selection (use of oligo-dT).

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Chu2022" class="csl-entry">

Chu, Tinyi, Zhong Wang, Dana Pe’er, and Charles G. Danko. 2022. “<span
class="nocase">Cell type and gene expression deconvolution with
BayesPrism enables Bayesian integrative analysis across bulk and
single-cell RNA sequencing in oncology</span>.” *Nature Cancer 2022 3:4*
3 (4): 505–17. <https://doi.org/10.1038/s43018-022-00356-3>.

</div>

<div id="ref-Kieffer2020" class="csl-entry">

Kieffer, Yann, Hocine R. Hocine, Géraldine Gentric, Floriane Pelon,
Charles Bernard, Brigitte Bourachot, Sonia Lameiras, et al. 2020. “<span
class="nocase">Single-cell analysis reveals fibroblast clusters linked
to immunotherapy resistance in cancer</span>.” *Cancer Discovery* 10
(9): 1330–51.
<https://doi.org/10.1158/2159-8290.CD-19-1384/333435/AM/SINGLE-CELL-ANALYSIS-REVEALS-FIBROBLAST-CLUSTERS>.

</div>

<div id="ref-Li2022" class="csl-entry">

Li, Yumei, Xinzhou Ge, Fanglue Peng, Wei Li, and Jingyi Jessica Li.
2022. “<span class="nocase">Exaggerated false positives by popular
differential expression methods when analyzing human population
samples</span>.” *Genome Biology* 23 (1): 1–13.
<https://doi.org/10.1186/S13059-022-02648-4/FIGURES/2>.

</div>

<div id="ref-Newman2015" class="csl-entry">

Newman, Aaron M., Chih Long Liu, Michael R. Green, Andrew J. Gentles,
Weiguo Feng, Yue Xu, Chuong D. Hoang, Maximilian Diehn, and Ash A.
Alizadeh. 2015. “<span class="nocase">Robust enumeration of cell subsets
from tissue expression profiles</span>.” *Nature Methods 2015 12:5* 12
(5): 453–57. <https://doi.org/10.1038/nmeth.3337>.

</div>

<div id="ref-Pelon2020" class="csl-entry">

Pelon, Floriane, Brigitte Bourachot, Yann Kieffer, Ilaria Magagna, Fanny
Mermet-Meillon, Isabelle Bonnet, Ana Costa, et al. 2020. “<span
class="nocase">Cancer-associated fibroblast heterogeneity in axillary
lymph nodes drives metastases in breast cancer through complementary
mechanisms</span>.” *Nature Communications 2020 11:1* 11 (1): 1–20.
<https://doi.org/10.1038/s41467-019-14134-w>.

</div>

<div id="ref-Szolek2014" class="csl-entry">

Szolek, András, Benjamin Schubert, Christopher Mohr, Marc Sturm,
Magdalena Feldhahn, and Oliver Kohlbacher. 2014. “<span
class="nocase">OptiType: precision HLA typing from next-generation
sequencing data</span>.” *Bioinformatics* 30 (23): 3310–16.
<https://doi.org/10.1093/BIOINFORMATICS/BTU548>.

</div>

<div id="ref-Zhao2020" class="csl-entry">

Zhao, Shanrong, Zhan Ye, and Robert Stanton. 2020. “<span
class="nocase">Misuse of RPKM or TPM normalization when comparing across
samples and sequencing protocols</span>.” *RNA* 26 (8): 903.
<https://doi.org/10.1261/RNA.074922.120>.

</div>

</div>
