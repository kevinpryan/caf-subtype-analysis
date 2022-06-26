CAF subtype analysis
================
Kevin Ryan
6/26/2022

-   [Introduction](#introduction)

# Introduction

Cancer-associated fibroblasts (CAFs) are a heterogeneous cell type found
in the tumour microenvironment. They have a wide array of functions, and
tend to be immunosuppressive and cancer-promoting. There have been many
attempts to characterise subtypes of CAFs, with much transcriptomic
analysis being carried out in the Mechta-Grigoriou lab in Institut
Curie. They have identified 4 ‘subtypes’ which can be separated based on
the expression of different markers:

-   S1: FAPHigh CD29Med-High αSMAHigh PDPNHigh PDGFRβHigh
-   S2: FAPNeg CD29Low αSMANeg-Low PDPNLow PDGFRβLow
-   S3: FAPNeg-Low CD29Med αSMANeg-Low PDPNLow PDGFRβLow-Med
-   S4: FAPLow-Med CD29High αSMAHigh PDPNLow PDGFRβMed

(Kieffer et al, Nature Communications, 2020)

FACS gating strategies can be used to isolate these various subtypes.
The Mechta-Grigoriou group have done this and have generated bulk
RNA-sequencing data for the S1, S3 and S4 subtypes. They generated
scRNA-sequencing data for the S1 subtype. This data was deposited on the
European Genome Phenome Archive, and was accessed via a Data Transfer
Agreement.

The following summarises the data obtained:

The data was processed using nf-core/rnaseq version using the default
parameters.
