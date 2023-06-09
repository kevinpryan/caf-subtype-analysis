library(gplots)
library(RColorBrewer)
library(dendextend)
library(ggplot2)
library(biomaRt)
library(reshape2)
library(limma)
library(VennDiagram)
library(grid)
library(gridExtra)
library(goseq)
library(dplyr)
library(colorspace)

metadata <- read.csv(file = "~/Documents/PhD/subtypes/caf-subtype-analysis/intermediate_files/metadata/reformat_samples_extra_info.csv")

colnames(metadata)[1] <- "Mixture"
metadata$Mixture <- as.character(metadata$Mixture)
metadata$Condition <- ifelse(metadata$Condition == "Tumour", "CAF", "TAN")

cibersort_output <- read.csv("~/Documents/PhD/subtypes/caf-subtype-analysis/intermediate_files/cibersort/CIBERSORTx_Job14_Results.csv")
cibersort_output$Mixture <- gsub(pattern = "X", replacement = "", x = cibersort_output$Mixture)

cibersort_output_absolute <- read.csv("~/Documents/PhD/subtypes/caf-subtype-analysis/intermediate_files/cibersort/CIBERSORTx_Job17_Results_20230518_absolutemode.csv")
cibersort_output_absolute$Mixture <- gsub(pattern = "X", replacement = "", x = cibersort_output_absolute$Mixture)

cibersort_output_with_metadata <- full_join(cibersort_output, metadata)
hclust_data <- cibersort_output_with_metadata[,2:4]
rownames(hclust_data) <- cibersort_output_with_metadata$Mixture
hclust_test <- hclust(dist(hclust_data))

## cluster and label by grade

cibersort_output_with_metadata$Grade <- as.factor(cibersort_output_with_metadata$Grade)
grade <- levels(cibersort_output_with_metadata$Grade)

dend <- as.dendrogram(hclust_test)
labels_colors(dend) <-
  colorspace::rainbow_hcl(2)[sort_levels_values(
    as.numeric(cibersort_output_with_metadata$Grade)[order.dendrogram(dend)]
  )]
# labels(dend) <- paste(as.character(cibersort_output_with_metadata$Grade)[order.dendrogram(dend)],
#                       "(",labels(dend),")", 
#                       sep = "")

labels(dend) <- as.character(cibersort_output_with_metadata$Grade)[order.dendrogram(dend)]
plot(dend)


## cluster and label by patient, colour by condition

cibersort_output_with_metadata$Patient <- as.factor(cibersort_output_with_metadata$Patient)
cibersort_output_with_metadata$Condition <- as.factor(cibersort_output_with_metadata$Condition)
patient <- levels(cibersort_output_with_metadata$Patient)
condition <- levels(cibersort_output_with_metadata$Condition)
dend2 <- as.dendrogram(hclust_test)
labels_colors(dend2) <-
  colorspace::rainbow_hcl(2)[sort_levels_values(
    as.numeric(cibersort_output_with_metadata$Condition)[order.dendrogram(dend2)]
  )]
labels(dend2) <- as.character(cibersort_output_with_metadata$Patient)[order.dendrogram(dend2)]
plot(dend2, main = "Clustered CIBERSORT results", xlab = "Patient")
legend("topleft", legend = condition, fill = rainbow_hcl(2))

## cluster and label by patient, colour by histology

dend3 <- as.dendrogram(hclust_test)
cibersort_output_with_metadata$Histology <- as.factor(cibersort_output_with_metadata$Histology)
histology <- levels(cibersort_output_with_metadata$Histology)
labels(dend3) <- as.character(cibersort_output_with_metadata$Patient)[order.dendrogram(dend3)]
labels_colors(dend3) <-
  colorspace::rainbow_hcl(2)[sort_levels_values(
    as.numeric(cibersort_output_with_metadata$Histology)[order.dendrogram(dend3)]
  )]
plot(dend3, main = "Clustered CIBERSORT results", xlab = "Patient")
legend("topleft", legend = histology, fill = rainbow_hcl(2))


