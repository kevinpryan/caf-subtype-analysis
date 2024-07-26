dat <- data.frame(
  "Kazakova_no" = c(20000, 154),
  "Kazakova_yes" = c(343, 25),
  row.names = c("InHouse_no", "InHouse_yes"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("Kazakova_no", "Kazakova_yes")

dat

mosaicplot(dat,
           main = "Mosaic plot",
           color = TRUE
)

test <- fisher.test(dat)
test


### edit tpm for online version of epic

star_salmon_tpm <- read.table("/home/rstudio/Documents/PhD/CAF_data/nfcore_results/inhouse_data_nfcore_results_version_3_8_1/star_salmon/salmon.merged.gene_tpm.tsv", header = TRUE)
star_salmon_tpm <- star_salmon_tpm %>% dplyr::select(-gene_id)
write.table(star_salmon_tpm, file = "/home/rstudio/Desktop/star_salmon_tpm_inhouse_no_gene_id.txt", quote = F, sep = "\t", row.names = F)

star_salmon_tpm_3808 <- read.table("/home/rstudio/Documents/PhD/CAF_data/nfcore_results/EGAD00001003808_nfcore_results//star_salmon/salmon.merged.gene_tpm.tsv", header = TRUE)
star_salmon_tpm_3808 <- star_salmon_tpm_3808 %>% dplyr::select(-gene_id)
write.table(star_salmon_tpm_3808, file = "/home/rstudio/Desktop/star_salmon_tpm_3808_no_gene_id.txt", quote = F, sep = "\t", row.names = F)

star_salmon_tpm_6144 <- read.table("/home/rstudio/Documents/PhD/CAF_data/nfcore_results/EGAD00001006144_nfcore_results/star_salmon/salmon.merged.gene_tpm.tsv", header = TRUE)
star_salmon_tpm_6144 <- star_salmon_tpm_6144 %>% dplyr::select(-gene_id)
write.table(star_salmon_tpm_3808, file = "/home/rstudio/Desktop/star_salmon_tpm_6144_no_gene_id.txt", quote = F, sep = "\t", row.names = F)

star_salmon_tpm_4810 <- read.table("/home/rstudio/Documents/PhD/CAF_data/nfcore_results/EGAD00001004810_nfcore_results/star_salmon/salmon.merged.gene_tpm.tsv", header = TRUE)
star_salmon_tpm_4810 <- star_salmon_tpm_4810 %>% dplyr::select(-gene_id)
write.table(star_salmon_tpm_4810, file = "/home/rstudio/Desktop/star_salmon_tpm_4810_no_gene_id.txt", quote = F, sep = "\t", row.names = F)

star_salmon_tpm_5744 <- read.table("/home/rstudio/Documents/PhD/CAF_data/nfcore_results/EGAD00001005744_nfcore_results/star_salmon/salmon.merged.gene_tpm.tsv", header = TRUE)
star_salmon_tpm_5744 <- star_salmon_tpm_5744 %>% dplyr::select(-gene_id)
write.table(star_salmon_tpm_5744, file = "/home/rstudio/Desktop/star_salmon_tpm_5744_no_gene_id.txt", quote = F, sep = "\t", row.names = F)
