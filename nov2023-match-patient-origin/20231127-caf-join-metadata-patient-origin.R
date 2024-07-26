library(dplyr)
matched_patients_data_path <- "~/Documents/PhD/CAF_data/caf_subtype_data_patient_of_origin_20230311.csv"
metadata_full_path <- "~/Documents/PhD/subtypes/caf-subtype-analysis/nf-subpop/outdir/metadata/metadata_full.txt"
matched_patients <- read.csv(matched_patients_data_path)
metadata <- read.table(metadata_full_path)
metadata <- metadata %>% rownames_to_column(var = "Sample")
metadata$Samplenames_old <- metadata$Sample
metadata$Sample <- gsub("CAF_Culture_", replacement = "", x = metadata$Sample)
metadata_inhouse <- metadata %>% dplyr::filter(Study == "InHouse")
metadata_no_inhouse <- metadata %>% dplyr::filter(Study != "InHouse")
inhouse_with_patient_in <- read.csv("~/Documents/PhD/CAF_data/InHouse/reformat_samples.csv")
colnames(inhouse_with_patient_in)[1] <- "Sample"
inhouse_with_patient_in <- inhouse_with_patient_in[,1:2]
inhouse_with_patient_in$Sample <- as.character(inhouse_with_patient_in$Sample)
metadata_inhouse <- full_join(inhouse_with_patient_in, metadata_inhouse, by = "Sample")
write.csv(metadata_inhouse, "~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/metadata_inhouse_with_patient.csv", row.names = F, quote = F)
matched_patients_cols_interest <- matched_patients[,1:2]
metadata_no_inhouse_with_patient <- left_join(metadata_no_inhouse, matched_patients_cols_interest, by = "Sample")
write.csv(metadata_no_inhouse_with_patient, "~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/metadata_no_inhouse_with_patient_20231127.csv", row.names = F, quote = F)
metadata_full_with_patient <- rbind.data.frame(metadata_inhouse, metadata_no_inhouse_with_patient)
write.csv(metadata_full_with_patient, "~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/metadata_full_with_patient_20231127.csv", row.names = F, quote = F)
